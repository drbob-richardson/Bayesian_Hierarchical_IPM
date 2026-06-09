source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
replicates <- as.integer(Sys.getenv("BENCHMARK_REPLICATES", "50"))
posterior_draws <- as.integer(Sys.getenv("BENCHMARK_DRAWS", "500"))
populations <- 3L
transitions <- 4L
initial_n <- 500L
profile_detection <- 0.25
capture_detection <- 0.35

task_directory <- "results/paired_benchmark/tasks"
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)

settings <- data.frame(
  replicates = replicates,
  posterior_draws = posterior_draws,
  populations = populations,
  transitions = transitions,
  initial_n = initial_n,
  profile_detection = profile_detection,
  capture_detection = capture_detection
)
write.csv(settings, "results/paired_benchmark/settings.csv", row.names = FALSE)

fit_method <- function(
  analysis,
  log_posterior,
  prior,
  replicate,
  method_index,
  valid_parameter
) {
  tryCatch({
    fit <- fit_laplace_custom(
      log_posterior,
      prior,
      analysis = analysis,
      maxit = 1200L
    )
    draws <- sample_custom_laplace_draws(
      fit,
      log_posterior,
      valid_parameter = valid_parameter,
      draws = posterior_draws,
      seed = 20274100 + 100L * replicate + method_index
    )
    summary <- summarize_inverse_ipm_draws(draws, mesh)
    summary$parameter$prior_sd <- prior$sd[summary$parameter$quantity]
    summary$parameter$contraction <- 1 -
      summary$parameter$posterior_sd / summary$parameter$prior_sd
    summary$parameter$group <- vapply(
      summary$parameter$quantity,
      parameter_group,
      character(1)
    )
    summary$derived$group <- vapply(
      summary$derived$quantity,
      parameter_group,
      character(1)
    )

    list(
      analysis = analysis,
      fit = fit,
      draws = draws,
      summary = summary,
      error = NULL
    )
  }, error = function(error) {
    list(
      analysis = analysis,
      fit = NULL,
      draws = NULL,
      summary = NULL,
      error = conditionMessage(error)
    )
  })
}

run_replicate <- function(replicate) {
  output_path <- file.path(
    task_directory,
    sprintf("replicate_%03d.rds", replicate)
  )
  if (file.exists(output_path)) {
    existing <- readRDS(output_path)
    required_analyses <- c(
      "profile_only",
      "profile_gamma_poisson",
      "capture_recapture",
      "integrated",
      "integrated_gamma_poisson",
      "oracle"
    )
    complete <- all(required_analyses %in% names(existing$analyses)) &&
      all(vapply(existing$analyses[required_analyses], function(analysis) {
      is.null(analysis$error) &&
        !is.null(analysis$draws) &&
        nrow(analysis$draws) >= posterior_draws
      }, logical(1)))
    if (complete) {
      return(output_path)
    }
  }

  set.seed(20274000 + replicate)
  simulations <- lapply(seq_len(populations), function(population_index) {
    simulate_population(
      years = transitions + 1L,
      initial_population = simulate_initial_population(
        n = initial_n,
        latent_quality_sd = 0
      ),
      vital_rates = default_fish_vital_rates(),
      environment = rep(0, transitions)
    )
  })

  profiles <- lapply(simulations, function(simulation) {
    sample_size_profiles(
      simulation$census,
      detection = function(size, year) {
        rep(profile_detection, length(size))
      },
      measurement_sd = 0
    )
  })
  counts <- lapply(profiles, make_profile_count_matrix, mesh = mesh)

  captures <- lapply(simulations, function(simulation) {
    sample_mark_recapture(
      simulation$census,
      detection = function(size, year) {
        rep(capture_detection, length(size))
      },
      measurement_sd = 0
    )
  })
  capture_events <- do.call(rbind, lapply(captures, function(capture) {
    prepare_capture_recapture_events(
      capture,
      final_year = transitions,
      recapture_detection = capture_detection
    )
  }))
  histories <- lapply(simulations, prepare_complete_history)

  prior <- default_inverse_ipm_priors("weak", "poisson")
  gamma_prior <- default_inverse_ipm_priors("weak", "gamma_poisson")
  valid_parameter <- function(parameter) {
    isTRUE(tryCatch({
      build_centered_ipm(parameter, mesh)
      TRUE
    }, error = function(error) FALSE))
  }
  profile_log_posterior <- function(parameter) {
    inverse_ipm_log_posterior(
      parameter,
      counts,
      mesh,
      prior,
      process_model = "poisson"
    )
  }
  capture_log_posterior <- function(parameter) {
    capture_recapture_log_likelihood(parameter, capture_events) +
      inverse_ipm_log_prior(parameter, prior)
  }
  integrated_log_posterior <- function(parameter) {
    inverse_ipm_log_posterior(
      parameter,
      counts,
      mesh,
      prior,
      process_model = "poisson"
    ) + capture_recapture_log_likelihood(parameter, capture_events)
  }
  oracle_log_posterior <- function(parameter) {
    complete_history_log_likelihood(parameter, histories) +
      inverse_ipm_log_prior(parameter, prior)
  }
  gamma_profile_log_posterior <- function(parameter) {
    inverse_ipm_log_posterior(
      parameter,
      counts,
      mesh,
      gamma_prior,
      process_model = "gamma_poisson"
    )
  }
  gamma_integrated_log_posterior <- function(parameter) {
    inverse_ipm_log_posterior(
      parameter,
      counts,
      mesh,
      gamma_prior,
      process_model = "gamma_poisson"
    ) + capture_recapture_log_likelihood(parameter, capture_events)
  }

  analyses <- list(
    profile_only = fit_method(
      "profile_only",
      profile_log_posterior,
      prior,
      replicate,
      1L,
      valid_parameter
    ),
    profile_gamma_poisson = fit_method(
      "profile_gamma_poisson",
      gamma_profile_log_posterior,
      gamma_prior,
      replicate,
      5L,
      valid_parameter
    ),
    capture_recapture = fit_method(
      "capture_recapture",
      capture_log_posterior,
      prior,
      replicate,
      2L,
      valid_parameter
    ),
    integrated = fit_method(
      "integrated",
      integrated_log_posterior,
      prior,
      replicate,
      3L,
      valid_parameter
    ),
    integrated_gamma_poisson = fit_method(
      "integrated_gamma_poisson",
      gamma_integrated_log_posterior,
      gamma_prior,
      replicate,
      6L,
      valid_parameter
    ),
    oracle = fit_method(
      "oracle",
      oracle_log_posterior,
      prior,
      replicate,
      4L,
      valid_parameter
    )
  )

  result <- list(
    replicate = replicate,
    observed = list(
      profile_counts = vapply(counts, sum, numeric(1)),
      captures = vapply(captures, nrow, integer(1)),
      capture_events = nrow(capture_events),
      consecutive_recaptures = sum(capture_events$recaptured)
    ),
    analyses = analyses
  )
  saveRDS(result, output_path, compress = "xz")
  output_path
}

available_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN",
  intern = TRUE
)))
if (!is.finite(available_cores)) {
  available_cores <- detectCores(logical = TRUE)
}
cores <- max(1L, min(6L, available_cores - 1L))

cat("Running", replicates, "paired benchmark replicates on", cores, "cores.\n")
paths <- mclapply(
  seq_len(replicates),
  run_replicate,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/paired_benchmark/task_paths.rds")
cat("Completed paired data benchmark.\n")

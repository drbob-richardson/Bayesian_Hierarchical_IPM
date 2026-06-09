source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")
source("R/profile_process_alternatives.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
replicates <- as.integer(Sys.getenv("PROFILE_ALT_REPLICATES", "50"))
posterior_draws <- as.integer(Sys.getenv("PROFILE_ALT_DRAWS", "500"))
experiment <- Sys.getenv("PROFILE_ALT_EXPERIMENT", "main")

scenarios <- data.frame(
  scenario = c(
    "standard",
    "dense_profiles",
    "longer_followup",
    "diverse_initial",
    "rich_design",
    "parent_excitation_sparse",
    "parent_excitation_dense"
  ),
  populations = c(3L, 3L, 3L, 3L, 6L, 8L, 8L),
  transitions = c(4L, 4L, 8L, 4L, 6L, 2L, 2L),
  initial_n = c(500L, 500L, 500L, 500L, 500L, 500L, 500L),
  profile_detection = c(0.25, 1.00, 0.25, 0.25, 0.50, 0.25, 1.00),
  initial_design = c(
    "baseline",
    "baseline",
    "baseline",
    "diverse",
    "baseline",
    "parent_excitation",
    "parent_excitation"
  ),
  stringsAsFactors = FALSE
)
methods <- c(
  "recursive_poisson",
  "recursive_gamma_poisson",
  "one_step_poisson",
  "one_step_gamma_poisson",
  "one_step_shared_gamma_poisson",
  "one_step_finite_population",
  "one_step_finite_population_overdispersed",
  "one_step_gamma_poisson_shape_informed",
  "one_step_gamma_poisson_fixed_recruitment_slope"
)
requested_scenarios <- Sys.getenv("PROFILE_ALT_SCENARIOS", "")
if (nzchar(requested_scenarios)) {
  requested_scenarios <- strsplit(requested_scenarios, ",", fixed = TRUE)[[1L]]
  scenarios <- scenarios[scenarios$scenario %in% requested_scenarios, ]
}
requested_methods <- Sys.getenv("PROFILE_ALT_METHODS", "")
if (nzchar(requested_methods)) {
  requested_methods <- strsplit(requested_methods, ",", fixed = TRUE)[[1L]]
  methods <- methods[methods %in% requested_methods]
}
tasks <- merge(
  scenarios,
  data.frame(replicate = seq_len(replicates)),
  all = TRUE
)
tasks$task_id <- seq_len(nrow(tasks))

results_directory <- file.path("results/profile_process_alternatives", experiment)
task_directory <- file.path(results_directory, "tasks")
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(tasks, file.path(results_directory, "task_grid.csv"), row.names = FALSE)
write.csv(
  data.frame(
    experiment = experiment,
    replicates = replicates,
    posterior_draws = posterior_draws,
    scenarios = paste(scenarios$scenario, collapse = ","),
    methods = paste(methods, collapse = ",")
  ),
  file.path(results_directory, "settings.csv"),
  row.names = FALSE
)

make_initial_population <- function(n, design, population_index) {
  initial <- simulate_initial_population(n = n, latent_quality_sd = 0)
  if (design == "baseline") {
    return(initial)
  }
  if (design == "parent_excitation") {
    center <- seq(60, 140, length.out = 8L)[population_index]
    initial$size <- pmin(145, pmax(35, rnorm(n, center, 4)))
    initial$age <- pmax(0L, round((initial$size - 35) / 18))
    return(initial)
  }

  component_probabilities <- list(
    c(0.90, 0.09, 0.01),
    c(0.62, 0.28, 0.10),
    c(0.20, 0.35, 0.45)
  )
  component <- sample.int(
    3L,
    n,
    replace = TRUE,
    prob = component_probabilities[[population_index]]
  )
  initial$size <- rnorm(
    n,
    mean = c(45, 80, 115)[component],
    sd = c(4, 8, 9)[component]
  )
  initial$size <- pmin(145, pmax(35, initial$size))
  initial$age <- pmax(0L, round((initial$size - 35) / 18))
  initial
}

fit_method <- function(
  counts,
  method,
  task_id,
  method_index,
  detection_probability
) {
  full_prior <- profile_process_prior(method, "weak")
  fixed <- profile_process_fixed_parameters(method)
  prior <- remove_fixed_from_prior(full_prior, fixed)
  log_posterior <- function(parameter) {
    full_parameter <- expand_profile_process_parameter(parameter, fixed)
    profile_process_log_posterior(
      full_parameter,
      counts,
      mesh,
      full_prior,
      method,
      detection_probability = detection_probability
    )
  }
  valid_parameter <- function(parameter) {
    full_parameter <- expand_profile_process_parameter(parameter, fixed)
    isTRUE(tryCatch({
      build_centered_ipm(full_parameter, mesh)
      TRUE
    }, error = function(error) FALSE))
  }

  tryCatch({
    fit <- fit_laplace_custom(
      log_posterior,
      prior,
      analysis = method,
      maxit = 1200L
    )
    draws <- sample_custom_laplace_draws(
      fit,
      log_posterior,
      valid_parameter = valid_parameter,
      draws = posterior_draws,
      seed = 20278100 + 10L * task_id + method_index
    )
    if (length(fixed) > 0L) {
      draws <- t(apply(draws, 1L, function(parameter) {
        expand_profile_process_parameter(parameter, fixed)
      }))
    }
    summary <- summarize_inverse_ipm_draws(draws, mesh)
    summary$parameter$prior_sd <- full_prior$sd[summary$parameter$quantity]
    summary$parameter$prior_sd[
      summary$parameter$quantity %in% names(fixed)
    ] <- NA_real_
    summary$parameter$contraction <- 1 -
      summary$parameter$posterior_sd / summary$parameter$prior_sd
    list(fit = fit, draws = draws, summary = summary, error = NULL)
  }, error = function(error) {
    list(
      fit = NULL,
      draws = NULL,
      summary = NULL,
      error = conditionMessage(error)
    )
  })
}

run_task <- function(task) {
  output_path <- file.path(
    task_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(output_path)) {
    existing <- readRDS(output_path)
    complete <- all(methods %in% names(existing$analyses)) &&
      all(vapply(existing$analyses[methods], function(analysis) {
        is.null(analysis$error) &&
          !is.null(analysis$draws) &&
          nrow(analysis$draws) >= posterior_draws
      }, logical(1)))
    if (complete) {
      return(output_path)
    }
  }

  set.seed(20278000 + task$task_id)
  counts <- lapply(seq_len(task$populations), function(population_index) {
    initial <- make_initial_population(
      task$initial_n,
      task$initial_design,
      population_index
    )
    simulation <- simulate_population(
      years = task$transitions + 1L,
      initial_population = initial,
      vital_rates = default_fish_vital_rates(),
      environment = rep(0, task$transitions)
    )
    profiles <- sample_size_profiles(
      simulation$census,
      detection = function(size, year) {
        rep(task$profile_detection, length(size))
      },
      measurement_sd = 0
    )
    make_profile_count_matrix(profiles, mesh)
  })

  analyses <- lapply(seq_along(methods), function(method_index) {
    fit_method(
      counts,
      methods[method_index],
      task$task_id,
      method_index,
      task$profile_detection
    )
  })
  names(analyses) <- methods

  result <- list(
    task = task,
    profile_counts = vapply(counts, sum, numeric(1)),
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

cat("Running", nrow(tasks), "profile-process tasks on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, file.path(results_directory, "task_paths.rds"))
cat("Completed profile-process alternatives.\n")

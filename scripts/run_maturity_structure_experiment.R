source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")
source("R/profile_process_alternatives.R")
source("R/maturity_structured_ipm.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
truth <- inverse_ipm_truth()
truth_derived <- maturity_inverse_ipm_quantities(truth, mesh)
replicates <- as.integer(Sys.getenv("MATURITY_REPLICATES", "50"))
posterior_draws <- as.integer(Sys.getenv("MATURITY_DRAWS", "500"))

scenarios <- data.frame(
  scenario = c("standard", "dense_profiles", "rich_design"),
  populations = c(3L, 3L, 6L),
  transitions = c(4L, 4L, 6L),
  initial_n = c(500L, 500L, 500L),
  profile_detection = c(0.25, 1.00, 0.50),
  stringsAsFactors = FALSE
)
methods <- c(
  "ungated_weak",
  "maturity_weak",
  "maturity_shape_informed",
  "maturity_fixed_recruitment_slope"
)
tasks <- merge(
  scenarios,
  data.frame(replicate = seq_len(replicates)),
  all = TRUE
)
tasks$task_id <- seq_len(nrow(tasks))

task_directory <- "results/maturity_structure/tasks"
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(tasks, "results/maturity_structure/task_grid.csv", row.names = FALSE)
write.csv(
  data.frame(replicates = replicates, posterior_draws = posterior_draws),
  "results/maturity_structure/settings.csv",
  row.names = FALSE
)

make_model_specification <- function(method) {
  full_prior <- default_inverse_ipm_priors("weak", "gamma_poisson")
  fixed <- numeric()
  kernel_builder <- build_maturity_centered_ipm
  draw_derived <- maturity_inverse_ipm_quantities

  if (method == "ungated_weak") {
    kernel_builder <- build_centered_ipm
    draw_derived <- derived_inverse_ipm_quantities
  } else if (method == "maturity_shape_informed") {
    full_prior$mean["recruitment_slope_20"] <-
      truth["recruitment_slope_20"]
    full_prior$sd["recruitment_slope_20"] <- 0.20
  } else if (method == "maturity_fixed_recruitment_slope") {
    fixed <- truth["recruitment_slope_20"]
  }

  list(
    full_prior = full_prior,
    prior = remove_fixed_from_prior(full_prior, fixed),
    fixed = fixed,
    kernel_builder = kernel_builder,
    draw_derived = draw_derived
  )
}

fit_method <- function(counts, method, task_id, method_index) {
  specification <- make_model_specification(method)
  log_posterior <- function(parameter) {
    full_parameter <- expand_profile_process_parameter(
      parameter,
      specification$fixed
    )
    one_step_profile_log_posterior(
      full_parameter,
      counts,
      mesh,
      specification$full_prior,
      "gamma_poisson",
      specification$kernel_builder
    )
  }
  valid_parameter <- function(parameter) {
    full_parameter <- expand_profile_process_parameter(
      parameter,
      specification$fixed
    )
    isTRUE(tryCatch({
      specification$kernel_builder(full_parameter, mesh)
      TRUE
    }, error = function(error) FALSE))
  }

  tryCatch({
    fit <- fit_laplace_custom(
      log_posterior,
      specification$prior,
      analysis = method,
      maxit = 1200L
    )
    draws <- sample_custom_laplace_draws(
      fit,
      log_posterior,
      valid_parameter = valid_parameter,
      draws = posterior_draws,
      seed = 20279100 + 10L * task_id + method_index
    )
    if (length(specification$fixed) > 0L) {
      draws <- t(apply(draws, 1L, function(parameter) {
        expand_profile_process_parameter(parameter, specification$fixed)
      }))
    }
    summary <- summarize_draws_with_derived_target(
      draws,
      mesh,
      specification$draw_derived,
      truth_derived,
      truth
    )
    summary$parameter$prior_sd <-
      specification$full_prior$sd[summary$parameter$quantity]
    summary$parameter$prior_sd[
      summary$parameter$quantity %in% names(specification$fixed)
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

  set.seed(20279000 + task$task_id)
  counts <- lapply(seq_len(task$populations), function(population_index) {
    simulation <- simulate_population(
      years = task$transitions + 1L,
      initial_population = simulate_initial_population(
        n = task$initial_n,
        latent_quality_sd = 0
      ),
      vital_rates = fish_vital_rates_from_maturity_ipm(truth),
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
    fit_method(counts, methods[method_index], task$task_id, method_index)
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

cat("Running", nrow(tasks), "maturity-structure tasks on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/maturity_structure/task_paths.rds")
cat("Completed maturity-structure experiment.\n")

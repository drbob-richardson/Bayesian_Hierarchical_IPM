# Observation-robustness experiment. The recovery model assumes the size
# profile is a constant multiple of latent abundance (size-independent,
# time-constant detection), which the linear projection absorbs. This run
# stresses that assumption with (i) size-dependent detection and (ii) unknown
# time-varying survey effort, and measures vital-rate coverage / RMSE
# degradation relative to the absorbed-constant baseline. Data are generated
# from the individual-based simulator, matching the main experiments.

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

library(parallel)

replicates <- as.integer(Sys.getenv("DETECTION_REPLICATES", "20"))
core_cap <- as.integer(Sys.getenv("DETECTION_CORES", "8"))
mesh <- seq(25, 160, by = 5)
results_directory <- "results/detection_robustness"
dir.create(file.path(results_directory, "tasks"), recursive = TRUE, showWarnings = FALSE)

truth <- inverse_ipm_truth()
vital_rates <- fish_vital_rates_from_inverse_ipm(truth)

detection_regimes <- c("constant", "size_dependent", "time_varying_effort")
process_models <- c("poisson", "gamma_poisson")

grid <- expand.grid(
  regime = detection_regimes,
  process_model = process_models,
  replicate = seq_len(replicates),
  stringsAsFactors = FALSE
)
grid$task_id <- seq_len(nrow(grid))
write.csv(grid, file.path(results_directory, "task_grid.csv"), row.names = FALSE)

observe_with_regime <- function(census, regime, seed) {
  set.seed(seed)
  years <- sort(unique(census$year))
  if (regime == "constant") {
    sample_size_profiles(
      census,
      detection = function(size, year) rep(0.5, length(size))
    )
  } else if (regime == "size_dependent") {
    sample_size_profiles(
      census,
      detection = function(size, year) 0.2 + 0.7 * inv_logit((size - 70) / 12)
    )
  } else {
    effort <- runif(length(years), 0.3, 0.9)
    sample_size_profiles(
      census,
      detection = function(size, year) rep(0.6, length(size)),
      effort = effort
    )
  }
}

run_task <- function(task) {
  output_path <- file.path(
    results_directory, "tasks", sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(output_path)) {
    return(output_path)
  }
  set.seed(20300000 + task$task_id)
  counts <- lapply(seq_len(3L), function(population_index) {
    initial <- simulate_initial_population(n = 500, latent_quality_sd = 0)
    simulation <- simulate_population(
      years = 5L,
      initial_population = initial,
      vital_rates = vital_rates,
      environment = rep(0, 4L)
    )
    profiles <- observe_with_regime(
      simulation$census, task$regime,
      seed = 20301000 + task$task_id * 10L + population_index
    )
    make_profile_count_matrix(profiles, mesh)
  })
  prior <- default_inverse_ipm_priors("weak", task$process_model)
  result <- tryCatch({
    fit <- fit_laplace_inverse_ipm(counts, mesh, prior, task$process_model)
    summary <- monte_carlo_fit_summary(
      fit, draws = 1000L, seed = 20302000 + task$task_id, truth = truth
    )
    list(task = task, summary = summary, error = NULL)
  }, error = function(error) {
    list(task = task, summary = NULL, error = conditionMessage(error))
  })
  saveRDS(result, output_path, compress = "xz")
  output_path
}

available_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN", intern = TRUE
)))
if (!is.finite(available_cores)) available_cores <- detectCores(logical = TRUE)
cores <- max(1L, min(core_cap, available_cores - 1L))

cat(sprintf("Detection robustness: %d tasks on %d cores.\n", nrow(grid), cores))
paths <- mclapply(
  split(grid, seq_len(nrow(grid))), run_task,
  mc.cores = cores, mc.preschedule = FALSE
)
saveRDS(paths, file.path(results_directory, "task_paths.rds"))
cat("Completed detection-robustness experiment.\n")

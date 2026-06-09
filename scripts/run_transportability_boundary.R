source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/transportability.R")
source("R/monte_carlo.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
truths <- transportability_truths()
replicates <- 10L
bias_multipliers <- c(0, 0.25, 0.50, 0.75, 1.00, 1.50)
data_directory <- "results/transportability_boundary/data"
task_directory <- "results/transportability_boundary/tasks"
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)

data_grid <- expand.grid(
  truth_name = names(truths),
  replicate = seq_len(replicates),
  stringsAsFactors = FALSE
)
data_grid$data_id <- seq_len(nrow(data_grid))
write.csv(
  data_grid,
  "results/transportability_boundary/data_grid.csv",
  row.names = FALSE
)

make_dataset <- function(row) {
  path <- file.path(
    data_directory,
    sprintf("data_%04d.rds", row$data_id)
  )
  if (file.exists(path)) {
    return(path)
  }

  truth <- truths[[row$truth_name]]
  vital_rates <- fish_vital_rates_from_inverse_ipm(truth)
  set.seed(20273000 + row$data_id)
  counts <- lapply(seq_len(3L), function(population_index) {
    initial <- simulate_initial_population(n = 500, latent_quality_sd = 0)
    simulation <- simulate_population(
      years = 5L,
      initial_population = initial,
      vital_rates = vital_rates,
      environment = rep(0, 4L)
    )
    profiles <- sample_size_profiles(
      simulation$census,
      detection = function(size, year) rep(0.25, length(size)),
      measurement_sd = 0
    )
    make_profile_count_matrix(profiles, mesh)
  })
  saveRDS(
    list(counts = counts, truth = truth, truth_name = row$truth_name),
    path,
    compress = "xz"
  )
  path
}

borrow_tasks <- expand.grid(
  data_id = data_grid$data_id,
  borrowing = c("full", "robust"),
  bias_multiplier = bias_multipliers,
  stringsAsFactors = FALSE
)
weak_tasks <- data.frame(
  data_id = data_grid$data_id,
  borrowing = "none",
  bias_multiplier = NA_real_,
  stringsAsFactors = FALSE
)
tasks <- rbind(weak_tasks, borrow_tasks)
tasks <- merge(tasks, data_grid, by = "data_id", all.x = TRUE, sort = FALSE)
tasks$task_id <- seq_len(nrow(tasks))
write.csv(
  tasks,
  "results/transportability_boundary/task_grid.csv",
  row.names = FALSE
)

available_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN",
  intern = TRUE
)))
if (!is.finite(available_cores)) {
  available_cores <- detectCores(logical = TRUE)
}
cores <- max(1L, min(6L, available_cores - 1L))

cat("Generating", nrow(data_grid), "multi-truth datasets on", cores, "cores.\n")
data_paths <- mclapply(
  split(data_grid, seq_len(nrow(data_grid))),
  make_dataset,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(data_paths, "results/transportability_boundary/data_paths.rds")

run_task <- function(task) {
  output_path <- file.path(
    task_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(output_path)) {
    return(output_path)
  }

  data <- readRDS(file.path(
    data_directory,
    sprintf("data_%04d.rds", task$data_id)
  ))
  truth <- data$truth

  if (task$borrowing == "none") {
    prior <- default_inverse_ipm_priors("weak", "gamma_poisson")
  } else {
    prior <- make_fecundity_study_prior(
      sample_size = 800L,
      accuracy = "biased",
      borrowing = task$borrowing,
      process_model = "gamma_poisson",
      robust_weight = 0.8,
      truth = truth,
      bias = fecundity_transport_bias(task$bias_multiplier),
      seed = 20274000 + task$data_id +
        round(1000 * task$bias_multiplier)
    )
  }

  result <- tryCatch({
    if (task$borrowing == "robust") {
      fit <- fit_robust_laplace_inverse_ipm(
        data$counts,
        mesh,
        prior,
        "gamma_poisson"
      )
      draws <- sample_robust_laplace_draws(
        fit,
        draws = 1000L,
        seed = 20275000 + task$task_id
      )
      summary <- summarize_inverse_ipm_draws(
        draws,
        mesh,
        truth = truth
      )
      robust_weight <- fit$posterior_external_weight
    } else {
      fit <- fit_laplace_inverse_ipm(
        data$counts,
        mesh,
        prior,
        "gamma_poisson"
      )
      summary <- monte_carlo_fit_summary(
        fit,
        draws = 1000L,
        seed = 20275000 + task$task_id,
        truth = truth
      )
      robust_weight <- NA_real_
    }
    list(
      task = task,
      truth = truth,
      fit = fit,
      summary = summary,
      robust_weight = robust_weight,
      error = NULL
    )
  }, error = function(error) {
    list(
      task = task,
      truth = truth,
      fit = NULL,
      summary = NULL,
      robust_weight = NA_real_,
      error = conditionMessage(error)
    )
  })
  saveRDS(result, output_path, compress = "xz")
  output_path
}

cat("Running", nrow(tasks), "transport-boundary fits on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/transportability_boundary/task_paths.rds")
cat("Completed transportability-boundary experiment.\n")

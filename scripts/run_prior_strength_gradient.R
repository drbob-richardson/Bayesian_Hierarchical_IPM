source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
scenario_grid <- read.csv("results/monte_carlo/scenario_grid.csv")
selected <- scenario_grid[
  scenario_grid$design == "three_by_four" &
    scenario_grid$initial_n == 500 &
    scenario_grid$profile_detection == 0.25 &
    scenario_grid$initial_profile == "baseline" &
    scenario_grid$replicate <= 5L,
]

gradient <- expand.grid(
  accuracy = c("correct", "biased"),
  uncertainty_multiplier = c(0.5, 1, 2, 4),
  stringsAsFactors = FALSE
)
tasks <- merge(selected, gradient, all = TRUE)
weak_tasks <- selected
weak_tasks$accuracy <- "weak"
weak_tasks$uncertainty_multiplier <- NA_real_
tasks <- rbind(weak_tasks, tasks)
tasks$task_id <- seq_len(nrow(tasks))

output_directory <- "results/prior_strength/tasks"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(tasks, "results/prior_strength/task_grid.csv", row.names = FALSE)

run_task <- function(task) {
  output_path <- file.path(
    output_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(output_path)) {
    return(output_path)
  }

  source_result <- readRDS(file.path(
    "results/monte_carlo/tasks",
    sprintf("scenario_%04d.rds", task$scenario_id)
  ))
  counts <- source_result$models$gamma_poisson$fit$counts
  if (task$accuracy == "weak") {
    prior <- make_external_information_prior("weak", "gamma_poisson")
  } else {
    prior <- make_external_information_prior(
      if (task$accuracy == "correct") {
        "all_vital_rates"
      } else {
        "all_vital_rates_biased"
      },
      "gamma_poisson",
      uncertainty_multiplier = task$uncertainty_multiplier
    )
  }

  result <- tryCatch({
    fit <- fit_laplace_inverse_ipm(counts, mesh, prior, "gamma_poisson")
    summary <- monte_carlo_fit_summary(
      fit,
      draws = 1000L,
      seed = 20264000 + task$task_id
    )
    list(task = task, fit = fit, summary = summary, error = NULL)
  }, error = function(error) {
    list(task = task, fit = NULL, summary = NULL, error = conditionMessage(error))
  })
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
cat("Running", nrow(tasks), "prior-strength fits on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/prior_strength/task_paths.rds")
cat("Completed prior-strength gradient.\n")

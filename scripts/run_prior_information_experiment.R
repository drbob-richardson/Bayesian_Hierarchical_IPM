source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
regimes <- c(
  "weak",
  "survival_growth",
  "recruit_size",
  "fecundity",
  "all_vital_rates",
  "all_vital_rates_biased"
)

scenario_grid <- read.csv("results/monte_carlo/scenario_grid.csv")
selected <- scenario_grid[
  scenario_grid$initial_n == 500 &
    scenario_grid$profile_detection == 0.25 &
    scenario_grid$initial_profile == "baseline",
]
tasks <- merge(
  selected,
  data.frame(prior_regime = regimes, stringsAsFactors = FALSE),
  all = TRUE
)
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[tasks$replicate <= 5L, ]

output_directory <- "results/prior_information/tasks"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(tasks, "results/prior_information/task_grid.csv", row.names = FALSE)

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
  prior <- make_external_information_prior(
    task$prior_regime,
    "gamma_poisson",
    uncertainty_multiplier = 2
  )

  result <- tryCatch({
    fit <- fit_laplace_inverse_ipm(
      counts = counts,
      mesh = mesh,
      prior = prior,
      process_model = "gamma_poisson"
    )
    summary <- monte_carlo_fit_summary(
      fit,
      draws = 800L,
      seed = 20263000 + task$task_id
    )
    list(task = task, fit = fit, summary = summary, error = NULL)
  }, error = function(error) {
    list(task = task, fit = NULL, summary = NULL, error = conditionMessage(error))
  })

  saveRDS(result, output_path, compress = "xz")
  output_path
}

available_cores <- detectCores(logical = FALSE)
if (!is.finite(available_cores)) {
  available_cores <- detectCores(logical = TRUE)
}
system_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN",
  intern = TRUE
)))
if (length(system_cores) == 1L && is.finite(system_cores)) {
  available_cores <- max(available_cores, system_cores, na.rm = TRUE)
}
if (!is.finite(available_cores)) {
  available_cores <- 2L
}
cores <- max(1L, min(6L, available_cores - 1L))
cat("Running", nrow(tasks), "prior-information fits on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/prior_information/task_paths.rds")
cat("Completed prior-information experiment.\n")

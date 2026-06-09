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
    scenario_grid$initial_profile == "baseline",
]

external_design <- expand.grid(
  sample_size = c(50L, 200L, 800L),
  accuracy = c("correct", "biased"),
  borrowing = c("full", "robust"),
  stringsAsFactors = FALSE
)
tasks <- merge(selected, external_design, all = TRUE)
weak_tasks <- selected
weak_tasks$sample_size <- 0L
weak_tasks$accuracy <- "weak"
weak_tasks$borrowing <- "none"
tasks <- rbind(weak_tasks, tasks)
tasks$task_id <- seq_len(nrow(tasks))

output_directory <- "results/auxiliary_fecundity/tasks"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  tasks,
  "results/auxiliary_fecundity/task_grid.csv",
  row.names = FALSE
)

run_task <- function(task) {
  output_path <- file.path(
    output_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(output_path)) {
    existing <- readRDS(output_path)
    if (
      task$borrowing != "robust" ||
        !is.null(existing$robust_mixture_weight)
    ) {
      return(output_path)
    }
  }

  source_result <- readRDS(file.path(
    "results/monte_carlo/tasks",
    sprintf("scenario_%04d.rds", task$scenario_id)
  ))
  counts <- source_result$models$gamma_poisson$fit$counts

  prior <- if (task$accuracy == "weak") {
    default_inverse_ipm_priors("weak", "gamma_poisson")
  } else {
    make_fecundity_study_prior(
      sample_size = task$sample_size,
      accuracy = task$accuracy,
      borrowing = task$borrowing,
      process_model = "gamma_poisson",
      robust_weight = 0.8,
      seed = 20265000 + task$scenario_id + 10L * task$sample_size +
        ifelse(task$accuracy == "biased", 100000L, 0L)
    )
  }

  result <- tryCatch({
    if (task$borrowing == "robust") {
      fit <- fit_robust_laplace_inverse_ipm(
        counts,
        mesh,
        prior,
        "gamma_poisson"
      )
      draws <- sample_robust_laplace_draws(
        fit,
        draws = 1000L,
        seed = 20266000 + task$task_id
      )
      summary <- summarize_inverse_ipm_draws(draws, mesh)
      summary$diagnostics <- data.frame(
        process_model = "gamma_poisson",
        convergence = fit$convergence,
        positive_definite = fit$positive_definite,
        condition_number = fit$condition_number,
        laplace_draws = nrow(draws),
        parameter_coverage = mean(summary$parameter$covers_truth),
        derived_coverage = mean(summary$derived$covers_truth),
        parameter_mean_abs_z = mean(
          abs(summary$parameter$standardized_bias)
        ),
        derived_mean_abs_z = mean(abs(summary$derived$standardized_bias))
      )
      robust_responsibility <- fit$posterior_external_weight
    } else {
      fit <- fit_laplace_inverse_ipm(counts, mesh, prior, "gamma_poisson")
      summary <- monte_carlo_fit_summary(
        fit,
        draws = 1000L,
        seed = 20266000 + task$task_id
      )
      robust_responsibility <- NA_real_
    }

    list(
      task = task,
      fit = fit,
      summary = summary,
      robust_responsibility = robust_responsibility,
      robust_mixture_weight = robust_responsibility,
      error = NULL
    )
  }, error = function(error) {
    list(
      task = task,
      fit = NULL,
      summary = NULL,
      robust_responsibility = NA_real_,
      error = conditionMessage(error)
    )
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
cat("Running", nrow(tasks), "auxiliary-fecundity fits on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/auxiliary_fecundity/task_paths.rds")
cat("Completed auxiliary-fecundity experiment.\n")

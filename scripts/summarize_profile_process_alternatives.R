source("R/paired_benchmark.R")
source("R/profile_process_alternatives.R")

experiment <- Sys.getenv("PROFILE_ALT_EXPERIMENT", "main")
results_directory <- file.path("results/profile_process_alternatives", experiment)
task_paths <- list.files(
  file.path(results_directory, "tasks"),
  pattern = "\\.rds$",
  full.names = TRUE
)
results <- lapply(task_paths, readRDS)

parameter_rows <- list()
derived_rows <- list()
diagnostic_rows <- list()
observed_rows <- list()

for (result in results) {
  task <- result$task
  observed_rows[[length(observed_rows) + 1L]] <- data.frame(
    scenario = task$scenario,
    replicate = task$replicate,
    total_profile_counts = sum(result$profile_counts)
  )

  for (method in names(result$analyses)) {
    analysis <- result$analyses[[method]]
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      method = method,
      success = is.null(analysis$error),
      error = if (is.null(analysis$error)) "" else analysis$error,
      positive_definite = if (is.null(analysis$fit)) {
        NA
      } else {
        analysis$fit$positive_definite
      },
      condition_number = if (is.null(analysis$fit)) {
        NA_real_
      } else {
        analysis$fit$condition_number
      }
    )

    if (
      !is.null(analysis$error) ||
        !isTRUE(analysis$fit$positive_definite)
    ) {
      next
    }
    parameter <- analysis$summary$parameter
    parameter$scenario <- task$scenario
    parameter$replicate <- task$replicate
    parameter$method <- method
    parameter$group <- vapply(
      parameter$quantity,
      profile_process_parameter_group,
      character(1)
    )
    parameter_rows[[length(parameter_rows) + 1L]] <- parameter

    derived <- analysis$summary$derived
    derived$scenario <- task$scenario
    derived$replicate <- task$replicate
    derived$method <- method
    derived$group <- vapply(
      derived$quantity,
      profile_process_parameter_group,
      character(1)
    )
    derived_rows[[length(derived_rows) + 1L]] <- derived
  }
}

parameters <- do.call(rbind, parameter_rows)
derived <- do.call(rbind, derived_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
observed <- do.call(rbind, observed_rows)

summarize_split <- function(data, grouping, include_contraction = FALSE) {
  groups <- split(data, interaction(data[grouping], drop = TRUE))
  output <- lapply(groups, function(group) {
    result <- group[1L, grouping, drop = FALSE]
    result$n <- nrow(group)
    result$coverage <- mean(group$covers_truth, na.rm = TRUE)
    result$bias <- mean(group$posterior_mean - group$truth, na.rm = TRUE)
    result$rmse <- sqrt(mean(
      (group$posterior_mean - group$truth)^2,
      na.rm = TRUE
    ))
    result$mean_posterior_sd <- mean(group$posterior_sd, na.rm = TRUE)
    result$mean_abs_standardized_bias <- mean(
      abs(group$standardized_bias),
      na.rm = TRUE
    )
    if (include_contraction) {
      result$mean_contraction <- mean(group$contraction, na.rm = TRUE)
    }
    result
  })
  rownames(output) <- NULL
  do.call(rbind, output)
}

parameter_by_group <- summarize_split(
  parameters,
  c("scenario", "method", "group"),
  include_contraction = TRUE
)
derived_by_group <- summarize_split(
  derived,
  c("scenario", "method", "group")
)
derived_by_quantity <- summarize_split(
  derived,
  c("scenario", "method", "quantity")
)

write.csv(
  parameters,
  file.path(results_directory, "parameter_results.csv"),
  row.names = FALSE
)
write.csv(
  derived,
  file.path(results_directory, "derived_results.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(results_directory, "diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  observed,
  file.path(results_directory, "observed_data_summary.csv"),
  row.names = FALSE
)
write.csv(
  parameter_by_group,
  file.path(results_directory, "parameter_by_group.csv"),
  row.names = FALSE
)
write.csv(
  derived_by_group,
  file.path(results_directory, "derived_by_group.csv"),
  row.names = FALSE
)
write.csv(
  derived_by_quantity,
  file.path(results_directory, "derived_by_quantity.csv"),
  row.names = FALSE
)

cat("Fits:", nrow(diagnostics), "\n")
cat("Success rate:", mean(diagnostics$success), "\n")
cat("Positive-definite rate:", mean(
  diagnostics$positive_definite,
  na.rm = TRUE
), "\n\n")
cat("Parameter recovery by scenario, method, and group:\n")
print(parameter_by_group)

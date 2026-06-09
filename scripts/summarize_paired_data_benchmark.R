source("R/paired_benchmark.R")

task_paths <- list.files(
  "results/paired_benchmark/tasks",
  pattern = "\\.rds$",
  full.names = TRUE
)
results <- lapply(task_paths, readRDS)
methods <- c(
  "profile_only",
  "profile_gamma_poisson",
  "capture_recapture",
  "integrated",
  "integrated_gamma_poisson",
  "oracle"
)

parameter_rows <- list()
derived_rows <- list()
diagnostic_rows <- list()
observed_rows <- list()
row_index <- 1L

for (result in results) {
  observed_rows[[length(observed_rows) + 1L]] <- data.frame(
    replicate = result$replicate,
    profiles_observed = sum(result$observed$profile_counts),
    captures = sum(result$observed$captures),
    capture_events = result$observed$capture_events,
    consecutive_recaptures = result$observed$consecutive_recaptures
  )

  for (method in methods) {
    analysis <- result$analyses[[method]]
    diagnostic_rows[[row_index]] <- data.frame(
      replicate = result$replicate,
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
      },
      laplace_draws = if (is.null(analysis$draws)) {
        0L
      } else {
        nrow(analysis$draws)
      }
    )
    row_index <- row_index + 1L

    if (
      !is.null(analysis$error) ||
        !isTRUE(analysis$fit$positive_definite)
    ) {
      next
    }
    parameter <- analysis$summary$parameter
    parameter$replicate <- result$replicate
    parameter$method <- method
    parameter_rows[[length(parameter_rows) + 1L]] <- parameter

    derived <- analysis$summary$derived
    derived$replicate <- result$replicate
    derived$method <- method
    derived_rows[[length(derived_rows) + 1L]] <- derived
  }
}

parameters <- do.call(rbind, parameter_rows)
derived <- do.call(rbind, derived_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
observed <- do.call(rbind, observed_rows)
parameters$group <- vapply(parameters$quantity, parameter_group, character(1))
derived$group <- vapply(derived$quantity, parameter_group, character(1))

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
  c("method", "group"),
  include_contraction = TRUE
)
parameter_by_quantity <- summarize_split(
  parameters,
  c("method", "quantity"),
  include_contraction = TRUE
)
derived_by_group <- summarize_split(derived, c("method", "group"))
derived_by_quantity <- summarize_split(derived, c("method", "quantity"))

dir.create("results/paired_benchmark", recursive = TRUE, showWarnings = FALSE)
write.csv(
  parameters,
  "results/paired_benchmark/parameter_results.csv",
  row.names = FALSE
)
write.csv(
  derived,
  "results/paired_benchmark/derived_results.csv",
  row.names = FALSE
)
write.csv(
  diagnostics,
  "results/paired_benchmark/diagnostics.csv",
  row.names = FALSE
)
write.csv(
  observed,
  "results/paired_benchmark/observed_data_summary.csv",
  row.names = FALSE
)
write.csv(
  parameter_by_group,
  "results/paired_benchmark/parameter_by_group.csv",
  row.names = FALSE
)
write.csv(
  parameter_by_quantity,
  "results/paired_benchmark/parameter_by_quantity.csv",
  row.names = FALSE
)
write.csv(
  derived_by_group,
  "results/paired_benchmark/derived_by_group.csv",
  row.names = FALSE
)
write.csv(
  derived_by_quantity,
  "results/paired_benchmark/derived_by_quantity.csv",
  row.names = FALSE
)

cat("Fits:", nrow(diagnostics), "\n")
cat("Success rate:", mean(diagnostics$success), "\n")
cat("Positive-definite rate:", mean(
  diagnostics$positive_definite,
  na.rm = TRUE
), "\n\n")
cat("Parameter recovery by method and group:\n")
print(parameter_by_group)
cat("\nDerived recovery by method and group:\n")
print(derived_by_group)

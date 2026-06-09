task_grid <- read.csv("results/auxiliary_fecundity/task_grid.csv")
task_files <- file.path(
  "results/auxiliary_fecundity/tasks",
  sprintf("task_%04d.rds", task_grid$task_id)
)
stopifnot(all(file.exists(task_files)))
results <- lapply(task_files, readRDS)

diagnostics <- list()
parameters <- list()
derived <- list()

for (index in seq_along(results)) {
  result <- results[[index]]
  task <- result$task
  base <- data.frame(
    task_id = task$task_id,
    scenario_id = task$scenario_id,
    replicate = task$replicate,
    sample_size = task$sample_size,
    accuracy = task$accuracy,
    borrowing = task$borrowing,
    robust_responsibility = result$robust_responsibility,
    failed = !is.null(result$error),
    error = ifelse(is.null(result$error), NA, result$error),
    stringsAsFactors = FALSE
  )

  if (is.null(result$error)) {
    diagnostics[[index]] <- cbind(base, result$summary$diagnostics)
    parameters[[index]] <- cbind(
      base,
      positive_definite = result$fit$positive_definite,
      result$summary$parameter
    )
    derived[[index]] <- cbind(
      base,
      positive_definite = result$fit$positive_definite,
      result$summary$derived
    )
  } else {
    diagnostics[[index]] <- base
  }
}

diagnostics <- do.call(rbind, diagnostics)
parameters <- do.call(rbind, parameters)
derived <- do.call(rbind, derived)

parameters$error <- parameters$posterior_mean - parameters$truth
parameters$interval_width <- parameters$q975 - parameters$q025
derived$error <- derived$posterior_mean - derived$truth
derived$interval_width <- derived$q975 - derived$q025

weak <- derived[
  derived$accuracy == "weak",
  c("scenario_id", "quantity", "error", "interval_width")
]
names(weak)[3:4] <- c("weak_error", "weak_interval_width")
derived <- merge(
  derived,
  weak,
  by = c("scenario_id", "quantity"),
  all.x = TRUE,
  sort = FALSE
)
derived$absolute_error_ratio_to_weak <- abs(derived$error) /
  pmax(abs(derived$weak_error), 1e-8)
derived$width_ratio_to_weak <- derived$interval_width /
  derived$weak_interval_width

dir.create("results/auxiliary_fecundity", recursive = TRUE, showWarnings = FALSE)
write.csv(
  diagnostics,
  "results/auxiliary_fecundity/diagnostics.csv",
  row.names = FALSE
)
write.csv(
  parameters,
  "results/auxiliary_fecundity/parameter_summary.csv",
  row.names = FALSE
)
write.csv(
  derived,
  "results/auxiliary_fecundity/derived_summary.csv",
  row.names = FALSE
)

target <- derived[
  derived$positive_definite &
    derived$quantity %in% c("recruitment_50", "recruitment_80", "recruitment_110"),
]
target$regime <- ifelse(
  target$accuracy == "weak",
  "weak",
  paste(target$accuracy, target$borrowing, sep = "_")
)

aggregate_table <- do.call(rbind, lapply(
  split(target, list(target$regime, target$sample_size), drop = TRUE),
  function(group) {
    data.frame(
      regime = group$regime[1L],
      sample_size = group$sample_size[1L],
      coverage = mean(group$covers_truth),
      bias = mean(group$error),
      rmse = sqrt(mean(group$error^2)),
      mean_interval_width = mean(group$interval_width),
      median_absolute_error_ratio_to_weak = median(
        group$absolute_error_ratio_to_weak
      ),
      mean_robust_responsibility = mean(
        group$robust_responsibility,
        na.rm = TRUE
      )
    )
  }
))
write.csv(
  aggregate_table,
  "results/auxiliary_fecundity/aggregate.csv",
  row.names = FALSE
)

cat("Fits:", nrow(diagnostics), "\n")
cat("Failures:", sum(diagnostics$failed), "\n")
cat("Positive-definite Hessians:", mean(
  diagnostics$positive_definite,
  na.rm = TRUE
), "\n\n")
print(aggregate_table, row.names = FALSE)

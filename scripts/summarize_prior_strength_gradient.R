task_grid <- read.csv("results/prior_strength/task_grid.csv")
task_files <- file.path(
  "results/prior_strength/tasks",
  sprintf("task_%04d.rds", task_grid$task_id)
)
stopifnot(all(file.exists(task_files)))
results <- lapply(task_files, readRDS)

rows <- list()
diagnostics <- list()
for (index in seq_along(results)) {
  result <- results[[index]]
  task <- result$task
  if (!is.null(result$error)) {
    next
  }
  base <- data.frame(
    scenario_id = task$scenario_id,
    replicate = task$replicate,
    accuracy = task$accuracy,
    uncertainty_multiplier = task$uncertainty_multiplier,
    positive_definite = result$fit$positive_definite,
    stringsAsFactors = FALSE
  )
  rows[[index]] <- cbind(base, result$summary$derived)
  diagnostics[[index]] <- cbind(base, result$summary$diagnostics)
}
derived <- do.call(rbind, rows)
diagnostics <- do.call(rbind, diagnostics)
derived$interval_width <- derived$q975 - derived$q025
derived$relative_absolute_error <- abs(derived$posterior_mean - derived$truth) /
  pmax(abs(derived$truth), 0.1)

weak <- derived[
  derived$accuracy == "weak",
  c("scenario_id", "quantity", "interval_width", "relative_absolute_error")
]
names(weak)[3:4] <- c("weak_width", "weak_error")
derived <- merge(
  derived,
  weak,
  by = c("scenario_id", "quantity"),
  all.x = TRUE,
  sort = FALSE
)
derived$width_ratio_to_weak <- derived$interval_width / derived$weak_width
derived$error_ratio_to_weak <- derived$relative_absolute_error /
  pmax(derived$weak_error, 1e-8)

dir.create("results/prior_strength", recursive = TRUE, showWarnings = FALSE)
write.csv(derived, "results/prior_strength/derived_summary.csv", row.names = FALSE)
write.csv(
  diagnostics,
  "results/prior_strength/diagnostics.csv",
  row.names = FALSE
)

target <- derived[derived$quantity %in% c("recruitment_80", "lambda"), ]
aggregate_table <- aggregate(
  cbind(
    covers_truth,
    relative_absolute_error,
    width_ratio_to_weak,
    error_ratio_to_weak
  ) ~ accuracy + uncertainty_multiplier + quantity,
  target,
  mean
)
write.csv(
  aggregate_table,
  "results/prior_strength/aggregate.csv",
  row.names = FALSE
)
cat("Positive-definite Hessians:", mean(diagnostics$positive_definite), "\n\n")
print(aggregate_table, row.names = FALSE)

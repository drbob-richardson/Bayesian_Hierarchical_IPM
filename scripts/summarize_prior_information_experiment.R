task_grid <- read.csv("results/prior_information/task_grid.csv")
task_files <- file.path(
  "results/prior_information/tasks",
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
    design = task$design,
    replicate = task$replicate,
    prior_regime = task$prior_regime,
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

parameters$interval_width <- parameters$q975 - parameters$q025
parameters$absolute_error <- abs(parameters$posterior_mean - parameters$truth)
derived$interval_width <- derived$q975 - derived$q025
derived$absolute_error <- abs(derived$posterior_mean - derived$truth)
derived$relative_absolute_error <- derived$absolute_error /
  pmax(abs(derived$truth), 0.1)

weak_width <- derived[
  derived$prior_regime == "weak",
  c("scenario_id", "quantity", "interval_width", "absolute_error")
]
names(weak_width)[3:4] <- c("weak_interval_width", "weak_absolute_error")
derived <- merge(
  derived,
  weak_width,
  by = c("scenario_id", "quantity"),
  all.x = TRUE,
  sort = FALSE
)
derived$width_ratio_to_weak <- derived$interval_width /
  derived$weak_interval_width
derived$error_ratio_to_weak <- derived$absolute_error /
  pmax(derived$weak_absolute_error, 1e-8)

dir.create("results/prior_information", recursive = TRUE, showWarnings = FALSE)
write.csv(
  diagnostics,
  "results/prior_information/diagnostics.csv",
  row.names = FALSE
)
write.csv(
  parameters,
  "results/prior_information/parameter_summary.csv",
  row.names = FALSE
)
write.csv(
  derived,
  "results/prior_information/derived_summary.csv",
  row.names = FALSE
)

valid <- derived[derived$positive_definite, ]
aggregate_table <- do.call(rbind, lapply(
  split(valid, list(valid$prior_regime, valid$design), drop = TRUE),
  function(group) {
    data.frame(
      prior_regime = group$prior_regime[1L],
      design = group$design[1L],
      coverage = mean(group$covers_truth),
      mean_relative_absolute_error = mean(group$relative_absolute_error),
      median_width_ratio_to_weak = median(group$width_ratio_to_weak),
      median_error_ratio_to_weak = median(group$error_ratio_to_weak)
    )
  }
))
write.csv(
  aggregate_table,
  "results/prior_information/aggregate.csv",
  row.names = FALSE
)

cat("Fits:", nrow(diagnostics), "\n")
cat("Failures:", sum(diagnostics$failed), "\n")
cat("Positive-definite Hessians:", mean(
  diagnostics$positive_definite,
  na.rm = TRUE
), "\n\n")
print(aggregate(
  cbind(covers_truth, relative_absolute_error, width_ratio_to_weak) ~
    prior_regime,
  valid,
  mean
))

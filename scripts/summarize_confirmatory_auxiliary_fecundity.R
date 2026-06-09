task_grid <- read.csv("results/auxiliary_confirmatory/task_grid.csv")
task_files <- file.path(
  "results/auxiliary_confirmatory/tasks",
  sprintf("task_%04d.rds", task_grid$task_id)
)
stopifnot(all(file.exists(task_files)))
results <- lapply(task_files, readRDS)

derived <- list()
diagnostics <- list()
for (index in seq_along(results)) {
  result <- results[[index]]
  task <- result$task
  base <- data.frame(
    replicate = task$replicate,
    regime = task$regime,
    sample_size = task$sample_size,
    robust_weight = result$robust_weight,
    positive_definite = ifelse(
      is.null(result$error),
      result$fit$positive_definite,
      NA
    ),
    failed = !is.null(result$error),
    stringsAsFactors = FALSE
  )
  if (is.null(result$error)) {
    derived[[index]] <- cbind(base, result$summary$derived)
    diagnostics[[index]] <- base
  } else {
    diagnostics[[index]] <- base
  }
}
derived <- do.call(rbind, derived)
diagnostics <- do.call(rbind, diagnostics)
derived$error <- derived$posterior_mean - derived$truth
derived$interval_width <- derived$q975 - derived$q025

recruitment <- derived[
  derived$positive_definite &
    derived$quantity %in% c(
      "recruitment_50",
      "recruitment_80",
      "recruitment_110"
    ),
]
per_dataset <- do.call(rbind, lapply(
  split(recruitment, list(recruitment$regime, recruitment$replicate), drop = TRUE),
  function(group) {
    at_80 <- group[group$quantity == "recruitment_80", ]
    data.frame(
      regime = group$regime[1L],
      replicate = group$replicate[1L],
      curve_rmse = sqrt(mean(group$error^2)),
      coverage_at_80 = at_80$covers_truth,
      error_at_80 = at_80$error,
      width_at_80 = at_80$interval_width,
      robust_weight = group$robust_weight[1L]
    )
  }
))

bootstrap_mean_interval <- function(value, iterations = 5000L, seed = 1L) {
  set.seed(seed)
  means <- replicate(
    iterations,
    mean(sample(value, length(value), replace = TRUE))
  )
  c(mean = mean(value), quantile(means, c(0.025, 0.975)))
}

summary_table <- do.call(rbind, lapply(
  split(per_dataset, per_dataset$regime),
  function(group) {
    rmse_interval <- bootstrap_mean_interval(
      group$curve_rmse,
      seed = 20271000 + group$replicate[1L]
    )
    coverage_interval <- binom.test(
      sum(group$coverage_at_80),
      nrow(group)
    )$conf.int
    data.frame(
      regime = group$regime[1L],
      datasets = nrow(group),
      mean_curve_rmse = rmse_interval["mean"],
      curve_rmse_q025 = rmse_interval["2.5%"],
      curve_rmse_q975 = rmse_interval["97.5%"],
      coverage_at_80 = mean(group$coverage_at_80),
      coverage_q025 = coverage_interval[1L],
      coverage_q975 = coverage_interval[2L],
      mean_width_at_80 = mean(group$width_at_80),
      mean_robust_weight = if (all(is.na(group$robust_weight))) {
        NA_real_
      } else {
        mean(group$robust_weight, na.rm = TRUE)
      }
    )
  }
))

weak <- per_dataset[
  per_dataset$regime == "weak",
  c("replicate", "curve_rmse")
]
names(weak)[2L] <- "weak_curve_rmse"
paired <- merge(per_dataset, weak, by = "replicate")
paired_summary <- do.call(rbind, lapply(
  split(paired, paired$regime),
  function(group) {
    set.seed(20272000 + group$replicate[1L])
    ratio <- mean(group$curve_rmse) / mean(group$weak_curve_rmse)
    bootstrap_ratio <- replicate(
      5000L,
      {
        index <- sample.int(nrow(group), nrow(group), replace = TRUE)
        mean(group$curve_rmse[index]) /
          mean(group$weak_curve_rmse[index])
      }
    )
    data.frame(
      regime = group$regime[1L],
      ratio_of_mean_curve_rmse_to_weak = ratio,
      ratio_q025 = quantile(bootstrap_ratio, 0.025),
      ratio_q975 = quantile(bootstrap_ratio, 0.975)
    )
  }
))

dir.create(
  "results/auxiliary_confirmatory",
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(
  diagnostics,
  "results/auxiliary_confirmatory/diagnostics.csv",
  row.names = FALSE
)
write.csv(
  per_dataset,
  "results/auxiliary_confirmatory/per_dataset.csv",
  row.names = FALSE
)
write.csv(
  summary_table,
  "results/auxiliary_confirmatory/summary.csv",
  row.names = FALSE
)
write.csv(
  paired_summary,
  "results/auxiliary_confirmatory/paired_summary.csv",
  row.names = FALSE
)

cat("Fits:", nrow(diagnostics), "\n")
cat("Failures:", sum(diagnostics$failed), "\n")
cat("Positive-definite Hessians:", mean(
  diagnostics$positive_definite,
  na.rm = TRUE
), "\n\n")
print(summary_table, row.names = FALSE)
cat("\nPaired RMSE ratios:\n")
print(paired_summary, row.names = FALSE)

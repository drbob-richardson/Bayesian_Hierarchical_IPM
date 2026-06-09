task_grid <- read.csv("results/transportability_boundary/task_grid.csv")
task_files <- file.path(
  "results/transportability_boundary/tasks",
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
    data_id = task$data_id,
    truth_name = task$truth_name,
    replicate = task$replicate,
    borrowing = task$borrowing,
    bias_multiplier = task$bias_multiplier,
    robust_weight = unname(result$robust_weight),
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
  }
  diagnostics[[index]] <- base
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
recruitment$group_id <- paste(
  recruitment$truth_name,
  recruitment$borrowing,
  ifelse(
    is.na(recruitment$bias_multiplier),
    "weak",
    recruitment$bias_multiplier
  ),
  recruitment$replicate,
  sep = ":"
)
per_dataset <- do.call(rbind, lapply(
  split(recruitment, recruitment$group_id),
  function(group) {
    at_80 <- group[group$quantity == "recruitment_80", ]
    data.frame(
      truth_name = group$truth_name[1L],
      borrowing = group$borrowing[1L],
      bias_multiplier = group$bias_multiplier[1L],
      replicate = group$replicate[1L],
      curve_rmse = sqrt(mean(group$error^2)),
      coverage_at_80 = at_80$covers_truth,
      width_at_80 = at_80$interval_width,
      robust_weight = group$robust_weight[1L]
    )
  }
))

weak <- per_dataset[
  per_dataset$borrowing == "none",
  c("truth_name", "replicate", "curve_rmse", "coverage_at_80")
]
names(weak)[3:4] <- c("weak_curve_rmse", "weak_coverage_at_80")
paired <- merge(
  per_dataset[per_dataset$borrowing != "none", ],
  weak,
  by = c("truth_name", "replicate"),
  all.x = TRUE
)

bootstrap_boundary_summary <- function(group, iterations = 3000L) {
  group <- group[
    is.finite(group$curve_rmse) &
      is.finite(group$weak_curve_rmse) &
      !is.na(group$coverage_at_80) &
      !is.na(group$weak_coverage_at_80),
  ]
  stopifnot(nrow(group) >= 2L)
  set.seed(
    20276000 +
      round(1000 * group$bias_multiplier[1L]) +
      sum(utf8ToInt(group$truth_name[1L])) +
      ifelse(group$borrowing[1L] == "robust", 10000L, 0L)
  )
  ratio <- mean(group$curve_rmse) / mean(group$weak_curve_rmse)
  bootstrap_ratio <- replicate(
    iterations,
    {
      index <- sample.int(nrow(group), nrow(group), replace = TRUE)
      mean(group$curve_rmse[index]) /
        mean(group$weak_curve_rmse[index])
    }
  )
  coverage <- mean(group$coverage_at_80)
  coverage_interval <- binom.test(
    sum(group$coverage_at_80),
    nrow(group)
  )$conf.int
  data.frame(
    truth_name = group$truth_name[1L],
    borrowing = group$borrowing[1L],
    bias_multiplier = group$bias_multiplier[1L],
    datasets = nrow(group),
    rmse_ratio_to_profiles_only = ratio,
    rmse_ratio_q025 = quantile(bootstrap_ratio, 0.025),
    rmse_ratio_q975 = quantile(bootstrap_ratio, 0.975),
    coverage_at_80 = coverage,
    coverage_q025 = coverage_interval[1L],
    coverage_q975 = coverage_interval[2L],
    weak_coverage_at_80 = mean(group$weak_coverage_at_80),
    mean_robust_weight = if (all(is.na(group$robust_weight))) {
      NA_real_
    } else {
      mean(group$robust_weight, na.rm = TRUE)
    }
  )
}

boundary_summary <- do.call(rbind, lapply(
  split(
    paired,
    list(paired$truth_name, paired$borrowing, paired$bias_multiplier),
    drop = TRUE
  ),
  bootstrap_boundary_summary
))

safety_boundary <- do.call(rbind, lapply(
  split(
    boundary_summary,
    list(boundary_summary$truth_name, boundary_summary$borrowing),
    drop = TRUE
  ),
  function(group) {
    group <- group[order(group$bias_multiplier), ]
    safe <- group$rmse_ratio_q975 < 1 &
      group$coverage_at_80 >= group$weak_coverage_at_80
    data.frame(
      truth_name = group$truth_name[1L],
      borrowing = group$borrowing[1L],
      largest_tested_safe_bias = if (any(safe)) {
        max(group$bias_multiplier[safe])
      } else {
        NA_real_
      },
      first_bias_with_worse_coverage = if (any(
        group$coverage_at_80 < group$weak_coverage_at_80
      )) {
        min(group$bias_multiplier[
          group$coverage_at_80 < group$weak_coverage_at_80
        ])
      } else {
        NA_real_
      }
    )
  }
))

pooled_summary <- do.call(rbind, lapply(
  split(
    paired,
    list(paired$borrowing, paired$bias_multiplier),
    drop = TRUE
  ),
  function(group) {
    group <- group[
      is.finite(group$curve_rmse) &
        is.finite(group$weak_curve_rmse) &
        !is.na(group$coverage_at_80),
    ]
    data.frame(
      borrowing = group$borrowing[1L],
      bias_multiplier = group$bias_multiplier[1L],
      datasets = nrow(group),
      rmse_ratio_to_profiles_only =
        mean(group$curve_rmse) / mean(group$weak_curve_rmse),
      coverage_at_80 = mean(group$coverage_at_80),
      mean_robust_weight = if (all(is.na(group$robust_weight))) {
        NA_real_
      } else {
        mean(group$robust_weight, na.rm = TRUE)
      }
    )
  }
))

dir.create(
  "results/transportability_boundary",
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(
  diagnostics,
  "results/transportability_boundary/diagnostics.csv",
  row.names = FALSE
)
write.csv(
  per_dataset,
  "results/transportability_boundary/per_dataset.csv",
  row.names = FALSE
)
write.csv(
  boundary_summary,
  "results/transportability_boundary/boundary_summary.csv",
  row.names = FALSE
)
write.csv(
  safety_boundary,
  "results/transportability_boundary/safety_boundary.csv",
  row.names = FALSE
)
write.csv(
  pooled_summary,
  "results/transportability_boundary/pooled_summary.csv",
  row.names = FALSE
)

cat("Fits:", nrow(diagnostics), "\n")
cat("Failures:", sum(diagnostics$failed), "\n")
cat("Positive-definite Hessians:", mean(
  diagnostics$positive_definite,
  na.rm = TRUE
), "\n\n")
print(boundary_summary, row.names = FALSE)
cat("\nSafety boundary:\n")
print(safety_boundary, row.names = FALSE)
cat("\nPooled across biological truths:\n")
print(pooled_summary, row.names = FALSE)

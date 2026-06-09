derived <- read.csv("results/paired_benchmark/derived_results.csv")
parameters <- read.csv("results/paired_benchmark/parameter_results.csv")

artifact_directory <- "artifacts/paired_benchmark"
dir.create(artifact_directory, recursive = TRUE, showWarnings = FALSE)

method_order <- c(
  "profile_only",
  "profile_gamma_poisson",
  "capture_recapture",
  "integrated",
  "integrated_gamma_poisson",
  "oracle"
)
method_labels <- c(
  "Profiles: Poisson",
  "Profiles: Gamma-Poisson",
  "Capture-recapture",
  "Integrated: Poisson",
  "Integrated: Gamma-Poisson",
  "Oracle"
)
method_colors <- c(
  "#4C78A8", "#72B7B2", "#F58518", "#54A24B", "#ECA82C", "#B279A2"
)
group_order <- c("survival", "growth", "recruitment", "recruit_size", "lambda")
group_labels <- c(
  "Survival", "Growth", "Recruitment", "Recruit size", "Population growth"
)

coverage <- aggregate(
  covers_truth ~ method + group,
  derived,
  mean,
  na.rm = TRUE
)
coverage_matrix <- matrix(
  NA_real_,
  nrow = length(method_order),
  ncol = length(group_order),
  dimnames = list(method_labels, group_labels)
)
for (index in seq_len(nrow(coverage))) {
  method_index <- match(coverage$method[index], method_order)
  group_index <- match(coverage$group[index], group_order)
  if (!is.na(method_index) && !is.na(group_index)) {
    coverage_matrix[method_index, group_index] <- coverage$covers_truth[index]
  }
}

png(
  file.path(artifact_directory, "coverage_by_data_source.png"),
  width = 1800,
  height = 1050,
  res = 180
)
par(mar = c(7, 5, 2, 1))
barplot(
  coverage_matrix,
  beside = TRUE,
  col = method_colors,
  ylim = c(0, 1),
  ylab = "Empirical 95% interval coverage",
  las = 2,
  border = NA
)
abline(h = 0.95, lty = 2, col = "gray35")
legend(
  "topright",
  legend = method_labels,
  fill = method_colors,
  border = NA,
  bty = "n"
)
dev.off()

contraction <- aggregate(
  contraction ~ method + group,
  parameters,
  mean,
  na.rm = TRUE
)
parameter_groups <- c("survival", "growth", "recruitment", "recruit_size")
contraction_matrix <- matrix(
  NA_real_,
  nrow = length(method_order),
  ncol = length(parameter_groups),
  dimnames = list(
    method_labels,
    c("Survival", "Growth", "Recruitment", "Recruit size")
  )
)
for (index in seq_len(nrow(contraction))) {
  method_index <- match(contraction$method[index], method_order)
  group_index <- match(contraction$group[index], parameter_groups)
  if (!is.na(method_index) && !is.na(group_index)) {
    contraction_matrix[method_index, group_index] <-
      contraction$contraction[index]
  }
}

png(
  file.path(artifact_directory, "prior_to_posterior_contraction.png"),
  width = 1800,
  height = 1050,
  res = 180
)
par(mar = c(6, 5, 2, 1))
barplot(
  contraction_matrix,
  beside = TRUE,
  col = method_colors,
  ylim = range(c(0, contraction_matrix), na.rm = TRUE),
  ylab = "Mean prior-to-posterior SD contraction",
  las = 2,
  border = NA
)
abline(h = 0, col = "gray35")
legend(
  "bottomright",
  legend = method_labels,
  fill = method_colors,
  border = NA,
  bty = "n"
)
dev.off()

selected <- c(
  "survival_80",
  "growth_increment_80",
  "recruitment_80",
  "lambda"
)
selected_labels <- c(
  "Survival at 80 mm",
  "Growth increment at 80 mm",
  "Recruitment at 80 mm",
  "Population growth rate"
)
derived$squared_error <- (derived$posterior_mean - derived$truth)^2
rmse <- aggregate(
  squared_error ~ method + quantity,
  derived[derived$quantity %in% selected, ],
  mean
)
rmse$rmse <- sqrt(rmse$squared_error)

png(
  file.path(artifact_directory, "rmse_ratio_selected_quantities.png"),
  width = 1800,
  height = 1200,
  res = 180
)
par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))
for (quantity_index in seq_along(selected)) {
  quantity <- selected[quantity_index]
  values <- rmse$rmse[match(
    paste(method_order, quantity),
    paste(rmse$method, rmse$quantity)
  )]
  ratio <- values / values[1L]
  barplot(
    ratio,
    names.arg = method_labels,
    col = method_colors,
    border = NA,
    las = 2,
    ylab = "RMSE ratio vs profiles",
    main = selected_labels[quantity_index],
  ylim = c(0, max(1.1, ratio, na.rm = TRUE))
  )
  abline(h = 1, lty = 2, col = "gray35")
}
dev.off()

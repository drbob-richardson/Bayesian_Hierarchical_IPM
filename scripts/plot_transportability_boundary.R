boundary <- read.csv("results/transportability_boundary/boundary_summary.csv")
truths <- c("baseline", "low_survival", "flatter_fecundity")
truth_labels <- c("Baseline", "Lower survival", "Flatter fecundity")
colors <- c(full = "#E15759", robust = "#F28E2B")
symbols <- c(full = 16, robust = 17)

dir.create(
  "artifacts/transportability_boundary",
  recursive = TRUE,
  showWarnings = FALSE
)
png(
  "artifacts/transportability_boundary/borrowing_safety_boundary.png",
  width = 2100,
  height = 1800,
  res = 190
)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 3.5, 1))
for (truth_index in seq_along(truths)) {
  truth <- truths[truth_index]
  subset <- boundary[boundary$truth_name == truth, ]

  plot(
    NA,
    xlim = range(subset$bias_multiplier),
    ylim = range(c(0, 1.5, subset$rmse_ratio_q025, subset$rmse_ratio_q975)),
    xlab = "Transport-bias multiplier",
    ylab = "Recruitment-curve RMSE / profiles-only RMSE",
    main = paste(truth_labels[truth_index], "- accuracy")
  )
  abline(h = 1, lty = 3, col = "gray30")
  for (borrowing in names(colors)) {
    group <- subset[subset$borrowing == borrowing, ]
    group <- group[order(group$bias_multiplier), ]
    polygon(
      c(group$bias_multiplier, rev(group$bias_multiplier)),
      c(group$rmse_ratio_q025, rev(group$rmse_ratio_q975)),
      border = NA,
      col = adjustcolor(colors[[borrowing]], alpha.f = 0.12)
    )
    lines(
      group$bias_multiplier,
      group$rmse_ratio_to_profiles_only,
      type = "b",
      lwd = 2,
      pch = symbols[[borrowing]],
      col = colors[[borrowing]]
    )
  }

  plot(
    NA,
    xlim = range(subset$bias_multiplier),
    ylim = c(0, 1),
    xlab = "Transport-bias multiplier",
    ylab = "Recruitment-at-80 95% coverage",
    main = paste(truth_labels[truth_index], "- calibration")
  )
  abline(h = unique(subset$weak_coverage_at_80), lty = 3, col = "gray30")
  abline(h = 0.95, lty = 2, col = "gray60")
  for (borrowing in names(colors)) {
    group <- subset[subset$borrowing == borrowing, ]
    group <- group[order(group$bias_multiplier), ]
    polygon(
      c(group$bias_multiplier, rev(group$bias_multiplier)),
      c(group$coverage_q025, rev(group$coverage_q975)),
      border = NA,
      col = adjustcolor(colors[[borrowing]], alpha.f = 0.12)
    )
    lines(
      group$bias_multiplier,
      group$coverage_at_80,
      type = "b",
      lwd = 2,
      pch = symbols[[borrowing]],
      col = colors[[borrowing]]
    )
  }
}
legend(
  "bottom",
  inset = -0.08,
  xpd = NA,
  horiz = TRUE,
  bty = "n",
  lwd = 2,
  pch = symbols,
  col = colors,
  legend = c("Full borrowing", "Robust borrowing")
)
dev.off()

robust <- boundary[boundary$borrowing == "robust", ]
png(
  "artifacts/transportability_boundary/robust_external_weight.png",
  width = 1500,
  height = 950,
  res = 180
)
par(mar = c(5, 4.5, 4, 1))
plot(
  NA,
  xlim = range(robust$bias_multiplier),
  ylim = c(0, 1),
  xlab = "Transport-bias multiplier",
  ylab = "Posterior external-study component probability",
  main = "Can profiles detect non-transportability?"
)
truth_colors <- c("#4E79A7", "#59A14F", "#9467BD")
for (index in seq_along(truths)) {
  group <- robust[robust$truth_name == truths[index], ]
  group <- group[order(group$bias_multiplier), ]
  lines(
    group$bias_multiplier,
    group$mean_robust_weight,
    type = "b",
    lwd = 2,
    pch = 16,
    col = truth_colors[index]
  )
}
legend(
  "bottomleft",
  bty = "n",
  lwd = 2,
  pch = 16,
  col = truth_colors,
  legend = truth_labels
)
dev.off()

cat("Created transportability-boundary plots.\n")

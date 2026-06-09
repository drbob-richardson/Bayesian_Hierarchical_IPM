summary_table <- read.csv("results/auxiliary_confirmatory/summary.csv")
paired <- read.csv("results/auxiliary_confirmatory/paired_summary.csv")

regimes <- c(
  "weak",
  "correct_full_50",
  "correct_full_200",
  "correct_full_800",
  "biased_full_800",
  "biased_robust_800"
)
labels <- c(
  "Profiles only",
  "Correct study, n=50",
  "Correct study, n=200",
  "Correct study, n=800",
  "Biased study, full",
  "Biased study, robust"
)
colors <- c(
  "#BDBDBD",
  "#9ECAE1",
  "#4E79A7",
  "#225EA8",
  "#E15759",
  "#F28E2B"
)
summary_table <- summary_table[match(regimes, summary_table$regime), ]
paired <- paired[match(regimes, paired$regime), ]

dir.create(
  "artifacts/auxiliary_fecundity",
  recursive = TRUE,
  showWarnings = FALSE
)
png(
  "artifacts/auxiliary_fecundity/confirmatory_results.png",
  width = 2000,
  height = 1050,
  res = 180
)
par(mfrow = c(1, 2), mar = c(10, 4.5, 4, 1))
mid <- barplot(
  summary_table$mean_curve_rmse,
  names.arg = labels,
  las = 2,
  col = colors,
  ylim = c(0, max(summary_table$curve_rmse_q975) * 1.1),
  ylab = "Mean recruitment-curve RMSE",
  main = "Accuracy across 30 fresh datasets"
)
arrows(
  mid,
  summary_table$curve_rmse_q025,
  mid,
  summary_table$curve_rmse_q975,
  angle = 90,
  code = 3,
  length = 0.04
)

mid <- barplot(
  summary_table$coverage_at_80,
  names.arg = labels,
  las = 2,
  col = colors,
  ylim = c(0, 1),
  ylab = "Recruitment-at-80 95% coverage",
  main = "Calibration across 30 fresh datasets"
)
arrows(
  mid,
  summary_table$coverage_q025,
  mid,
  summary_table$coverage_q975,
  angle = 90,
  code = 3,
  length = 0.04
)
abline(h = 0.95, lty = 3, col = "gray30")
dev.off()

cat("Created confirmatory auxiliary-fecundity plot.\n")

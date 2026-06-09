derived <- read.csv("results/prior_information/derived_summary.csv")
derived <- derived[derived$positive_definite, ]

regimes <- c(
  "weak",
  "survival_growth",
  "recruit_size",
  "fecundity",
  "all_vital_rates",
  "all_vital_rates_biased"
)
regime_labels <- c(
  "Diffuse",
  "Survival +\ngrowth",
  "Recruit-size\ndistribution",
  "Fecundity",
  "All rates",
  "All rates,\nbiased"
)
colors <- c(
  "#BDBDBD",
  "#4E79A7",
  "#59A14F",
  "#F28E2B",
  "#9467BD",
  "#E15759"
)

dir.create("artifacts/prior_information", recursive = TRUE, showWarnings = FALSE)

coverage <- sapply(regimes, function(regime) {
  mean(derived$covers_truth[derived$prior_regime == regime])
})
width <- sapply(regimes, function(regime) {
  median(derived$width_ratio_to_weak[derived$prior_regime == regime])
})

png(
  "artifacts/prior_information/overall_prior_value.png",
  width = 1900,
  height = 950,
  res = 180
)
par(mfrow = c(1, 2), mar = c(9, 4, 4, 1))
barplot(
  coverage,
  names.arg = regime_labels,
  col = colors,
  las = 2,
  ylim = c(0, 1),
  ylab = "Mean 95% interval coverage",
  main = "Calibration"
)
abline(h = 0.95, lty = 3, col = "gray30")
barplot(
  width,
  names.arg = regime_labels,
  col = colors,
  las = 2,
  ylim = c(0, max(1, width)),
  ylab = "Median interval width / diffuse-prior width",
  main = "Precision gained from prior information"
)
abline(h = 1, lty = 3, col = "gray30")
dev.off()

quantities <- c(
  "survival_80",
  "growth_increment_80",
  "growth_sd.log_growth_sd",
  "recruitment_80",
  "recruit_mean_80",
  "recruit_sd.log_recruit_sd",
  "lambda"
)
quantity_labels <- c(
  "Survival",
  "Growth increment",
  "Growth SD",
  "Recruitment",
  "Recruit mean",
  "Recruit SD",
  "Lambda"
)
matrix_value <- sapply(regimes, function(regime) {
  sapply(quantities, function(quantity) {
    subset <- derived[
      derived$prior_regime == regime & derived$quantity == quantity,
    ]
    median(subset$error_ratio_to_weak)
  })
})
matrix_value <- pmin(matrix_value, 2)

png(
  "artifacts/prior_information/error_ratio_heatmap.png",
  width = 1750,
  height = 1200,
  res = 180
)
par(mar = c(8, 9, 4, 2))
image(
  x = seq_along(regimes),
  y = seq_along(quantities),
  z = t(matrix_value),
  col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  zlim = c(0, 2),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Posterior error relative to diffuse-prior analysis"
)
axis(1, at = seq_along(regimes), labels = regime_labels, las = 2)
axis(2, at = seq_along(quantities), labels = quantity_labels, las = 2)
for (i in seq_along(regimes)) {
  for (j in seq_along(quantities)) {
    text(i, j, sprintf("%.2f", matrix_value[j, i]), cex = 0.85)
  }
}
box()
dev.off()

cat("Created prior-information plots in artifacts/prior_information/.\n")

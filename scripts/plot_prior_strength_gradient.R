aggregate_table <- read.csv("results/prior_strength/aggregate.csv")
dir.create("artifacts/prior_information", recursive = TRUE, showWarnings = FALSE)

colors <- c(correct = "#4E79A7", biased = "#E15759")
quantities <- c(recruitment_80 = "Recruitment at 80 mm", lambda = "Lambda")

png(
  "artifacts/prior_information/prior_strength_gradient.png",
  width = 1900,
  height = 950,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1))
for (quantity in names(quantities)) {
  subset <- aggregate_table[aggregate_table$quantity == quantity, ]
  plot(
    NA,
    xlim = rev(range(subset$uncertainty_multiplier)),
    ylim = range(c(0.5, 1, subset$error_ratio_to_weak)),
    log = "x",
    xlab = "External-prior uncertainty multiplier\n(smaller is stronger)",
    ylab = "Mean posterior error / diffuse-prior error",
    main = quantities[[quantity]]
  )
  abline(h = 1, lty = 3, col = "gray30")
  for (accuracy in names(colors)) {
    group <- subset[subset$accuracy == accuracy, ]
    group <- group[order(group$uncertainty_multiplier), ]
    lines(
      group$uncertainty_multiplier,
      group$error_ratio_to_weak,
      type = "b",
      pch = 16,
      lwd = 2,
      col = colors[[accuracy]]
    )
  }
  if (quantity == "recruitment_80") {
    legend(
      "topleft",
      bty = "n",
      lwd = 2,
      pch = 16,
      col = colors,
      legend = c("Correctly centered", "Systematically biased")
    )
  }
}
dev.off()

cat("Created prior-strength plot.\n")

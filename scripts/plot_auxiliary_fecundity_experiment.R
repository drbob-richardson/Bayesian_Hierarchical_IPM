aggregate_table <- read.csv("results/auxiliary_fecundity/aggregate.csv")
diagnostics <- read.csv("results/auxiliary_fecundity/diagnostics.csv")

regimes <- c("correct_full", "correct_robust", "biased_full", "biased_robust")
labels <- c(
  correct_full = "Correct, full borrowing",
  correct_robust = "Correct, robust borrowing",
  biased_full = "Biased, full borrowing",
  biased_robust = "Biased, robust borrowing"
)
colors <- c(
  correct_full = "#4E79A7",
  correct_robust = "#59A14F",
  biased_full = "#E15759",
  biased_robust = "#F28E2B"
)
symbols <- c(correct_full = 16, correct_robust = 17, biased_full = 16, biased_robust = 17)

weak_rmse <- aggregate_table$rmse[aggregate_table$regime == "weak"]
weak_coverage <- aggregate_table$coverage[aggregate_table$regime == "weak"]
external <- aggregate_table[aggregate_table$regime != "weak", ]

dir.create("artifacts/auxiliary_fecundity", recursive = TRUE, showWarnings = FALSE)
png(
  "artifacts/auxiliary_fecundity/auxiliary_study_performance.png",
  width = 1900,
  height = 950,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1))
plot(
  NA,
  xlim = range(external$sample_size),
  ylim = range(c(external$rmse, weak_rmse)),
  log = "x",
  xlab = "Fish in independent fecundity study",
  ylab = "Recruitment-curve RMSE",
  main = "Accuracy"
)
abline(h = weak_rmse, lty = 3, col = "gray30")
for (regime in regimes) {
  group <- external[external$regime == regime, ]
  group <- group[order(group$sample_size), ]
  lines(
    group$sample_size,
    group$rmse,
    type = "b",
    lwd = 2,
    pch = symbols[[regime]],
    col = colors[[regime]]
  )
}
legend(
  "topright",
  bty = "n",
  lwd = 2,
  pch = symbols,
  col = colors,
  legend = labels
)

plot(
  NA,
  xlim = range(external$sample_size),
  ylim = c(0, 1),
  log = "x",
  xlab = "Fish in independent fecundity study",
  ylab = "Recruitment-curve 95% coverage",
  main = "Calibration"
)
abline(h = weak_coverage, lty = 3, col = "gray30")
abline(h = 0.95, lty = 2, col = "gray60")
for (regime in regimes) {
  group <- external[external$regime == regime, ]
  group <- group[order(group$sample_size), ]
  lines(
    group$sample_size,
    group$coverage,
    type = "b",
    lwd = 2,
    pch = symbols[[regime]],
    col = colors[[regime]]
  )
}
dev.off()

robust <- diagnostics[
  diagnostics$borrowing == "robust" & diagnostics$positive_definite,
]
responsibility <- aggregate(
  robust_responsibility ~ accuracy + sample_size,
  robust,
  mean
)

png(
  "artifacts/auxiliary_fecundity/robust_borrowing_weight.png",
  width = 1300,
  height = 900,
  res = 180
)
par(mar = c(5, 4.5, 4, 1))
plot(
  NA,
  xlim = range(responsibility$sample_size),
  ylim = c(0, 1),
  log = "x",
  xlab = "Fish in independent fecundity study",
  ylab = "Posterior external-study component probability",
  main = "Does robust borrowing detect non-transportability?"
)
for (accuracy in c("correct", "biased")) {
  group <- responsibility[responsibility$accuracy == accuracy, ]
  group <- group[order(group$sample_size), ]
  lines(
    group$sample_size,
    group$robust_responsibility,
    type = "b",
    lwd = 2,
    pch = 16,
    col = if (accuracy == "correct") "#59A14F" else "#E15759"
  )
}
legend(
  "bottomleft",
  bty = "n",
  lwd = 2,
  pch = 16,
  col = c("#59A14F", "#E15759"),
  legend = c("Correctly transportable", "Biased transport")
)
dev.off()

cat("Created auxiliary-fecundity plots.\n")

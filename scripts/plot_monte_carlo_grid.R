diagnostics <- read.csv("results/monte_carlo/diagnostics.csv")
parameters <- read.csv("results/monte_carlo/parameter_summary.csv")
derived <- read.csv("results/monte_carlo/derived_summary.csv")
aggregate_table <- read.csv("results/monte_carlo/aggregate.csv")

dir.create("artifacts/monte_carlo", recursive = TRUE, showWarnings = FALSE)

design_levels <- c("one_by_twelve", "three_by_four", "six_by_two")
design_labels <- c("1 x 12", "3 x 4", "6 x 2")
model_colors <- c(poisson = "#2C7FB8", gamma_poisson = "#F28E2B")

successful <- diagnostics[
  !diagnostics$failed & diagnostics$positive_definite,
]
parameters <- parameters[parameters$positive_definite, ]
derived <- derived[derived$positive_definite, ]

png(
  "artifacts/monte_carlo/coverage_by_design.png",
  width = 1900,
  height = 1000,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 1), oma = c(2, 0, 3, 0))
for (metric in c("parameter_coverage", "derived_coverage")) {
  values <- sapply(design_levels, function(design) {
    sapply(names(model_colors), function(model) {
      mean(successful[
        successful$design == design & successful$process_model == model,
        metric
      ])
    })
  })
  barplot(
    values,
    beside = TRUE,
    col = model_colors[rownames(values)],
    names.arg = design_labels,
    ylim = c(0, 1),
    ylab = "Mean 95% interval coverage",
    main = if (metric == "parameter_coverage") {
      "Vital-rate parameters"
    } else {
      "Derived demographic quantities"
    }
  )
  abline(h = 0.95, lty = 3, col = "gray40")
}
legend(
  "bottom",
  inset = -0.15,
  xpd = NA,
  horiz = TRUE,
  bty = "n",
  fill = model_colors,
  legend = c("Poisson", "Gamma-Poisson")
)
mtext(
  "Monte Carlo recovery across all population sizes, sampling rates, and initial profiles",
  outer = TRUE,
  line = 1,
  cex = 1.15,
  font = 2
)
dev.off()

png(
  "artifacts/monte_carlo/coverage_by_sample_size.png",
  width = 1900,
  height = 1000,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 1), oma = c(2, 0, 3, 0))
for (metric in c("parameter_coverage", "derived_coverage")) {
  combinations <- expand.grid(
    initial_n = sort(unique(successful$initial_n)),
    detection = sort(unique(successful$profile_detection))
  )
  values <- sapply(seq_len(nrow(combinations)), function(index) {
    mean(successful[
      successful$initial_n == combinations$initial_n[index] &
        successful$profile_detection == combinations$detection[index] &
        successful$process_model == "poisson",
      metric
    ])
  })
  names(values) <- paste0(
    combinations$initial_n,
    " / ",
    round(100 * combinations$detection),
    "%"
  )
  barplot(
    values,
    col = "#2C7FB8",
    ylim = c(0, 1),
    ylab = "Mean 95% interval coverage",
    main = if (metric == "parameter_coverage") {
      "Vital-rate parameters"
    } else {
      "Derived demographic quantities"
    }
  )
  abline(h = 0.95, lty = 3, col = "gray40")
}
mtext(
  "Poisson recovery by population size and profile sampling rate",
  outer = TRUE,
  line = 1,
  cex = 1.15,
  font = 2
)
dev.off()

selected_quantities <- c(
  "survival_80",
  "growth_increment_80",
  "growth_sd.log_growth_sd",
  "recruitment_80",
  "recruit_mean_80",
  "recruit_sd.log_recruit_sd",
  "lambda"
)
quantity_labels <- c(
  survival_80 = "Survival\nat 80 mm",
  growth_increment_80 = "Growth\nat 80 mm",
  growth_sd.log_growth_sd = "Growth SD",
  recruitment_80 = "Recruitment\nat 80 mm",
  recruit_mean_80 = "Recruit\nmean",
  recruit_sd.log_recruit_sd = "Recruit SD",
  lambda = "Lambda"
)
quantity_coverage <- sapply(selected_quantities, function(quantity) {
  sapply(names(model_colors), function(model) {
    mean(derived$covers_truth[
      derived$quantity == quantity & derived$process_model == model
    ])
  })
})

png(
  "artifacts/monte_carlo/coverage_by_quantity.png",
  width = 2100,
  height = 1050,
  res = 190
)
par(mar = c(6, 4, 4, 1))
barplot(
  quantity_coverage,
  beside = TRUE,
  col = model_colors[rownames(quantity_coverage)],
  names.arg = quantity_labels[colnames(quantity_coverage)],
  ylim = c(0, 1),
  ylab = "Mean 95% interval coverage",
  main = "Which demographic quantities are recovered?"
)
abline(h = 0.95, lty = 3, col = "gray40")
legend(
  "bottomleft",
  bty = "n",
  horiz = TRUE,
  fill = model_colors,
  legend = c("Poisson", "Gamma-Poisson")
)
dev.off()

profile_levels <- c("baseline", "cohort_pulse")
profile_colors <- c(baseline = "#4E79A7", cohort_pulse = "#59A14F")
profile_values <- sapply(design_levels, function(design) {
  sapply(profile_levels, function(profile) {
    mean(successful$derived_coverage[
      successful$design == design &
        successful$initial_profile == profile &
        successful$process_model == "poisson"
    ])
  })
})

png(
  "artifacts/monte_carlo/initial_profile_comparison.png",
  width = 1500,
  height = 950,
  res = 180
)
par(mar = c(5, 4, 4, 1))
barplot(
  profile_values,
  beside = TRUE,
  col = profile_colors[rownames(profile_values)],
  names.arg = design_labels,
  ylim = c(0, 1),
  ylab = "Mean derived-quantity coverage",
  main = "A common cohort pulse does not uniformly improve recovery"
)
abline(h = 0.95, lty = 3, col = "gray40")
legend(
  "bottomleft",
  bty = "n",
  horiz = TRUE,
  fill = profile_colors,
  legend = c("Baseline profile", "Small-fish cohort pulse")
)
dev.off()

target_quantities <- c(
  "survival_80",
  "growth_increment_80",
  "recruitment_80",
  "lambda"
)
target <- derived[
  derived$quantity %in% target_quantities &
    derived$process_model == "poisson",
]

png(
  "artifacts/monte_carlo/bias_distributions.png",
  width = 2100,
  height = 1600,
  res = 190
)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
for (quantity in target_quantities) {
  subset <- target[target$quantity == quantity, ]
  errors <- lapply(design_levels, function(design) {
    data <- subset[subset$design == design, ]
    (data$posterior_mean - data$truth) / pmax(abs(data$truth), 0.1)
  })
  boxplot(
    errors,
    names = design_labels,
    col = "#9ECAE1",
    ylab = "Relative posterior-mean error",
    main = quantity,
    outline = FALSE
  )
  abline(h = 0, lty = 2, col = "gray30")
}
dev.off()

audit_path <- "results/monte_carlo/mcmc_audit_comparison.csv"
if (file.exists(audit_path)) {
  audit <- read.csv(audit_path)
  png(
    "artifacts/monte_carlo/laplace_mcmc_audit.png",
    width = 1900,
    height = 900,
    res = 180
  )
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))
  plot(
    audit$posterior_sd_mcmc,
    audit$posterior_sd_laplace,
    pch = 16,
    col = adjustcolor("#2C7FB8", alpha.f = 0.55),
    xlab = "Full-MCMC posterior SD",
    ylab = "Laplace posterior SD",
    main = "Posterior uncertainty"
  )
  abline(0, 1, lty = 2, lwd = 2)
  audit_errors <- lapply(design_levels, function(design) {
    abs(audit$mean_difference_in_mcmc_sd[audit$design == design])
  })
  boxplot(
    audit_errors,
    names = design_labels,
    col = "#9ECAE1",
    ylab = "Absolute mean difference\n(in MCMC posterior SD units)",
    xlab = "",
    main = "Posterior location"
  )
  abline(h = 0.25, lty = 3, col = "gray40")
  dev.off()
}

cat("Created Monte Carlo plots in artifacts/monte_carlo/.\n")

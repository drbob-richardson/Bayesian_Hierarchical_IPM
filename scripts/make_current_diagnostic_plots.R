source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/correlated_process_ipm.R")

dir.create("artifacts", showWarnings = FALSE)

colors <- c(
  truth = "#1F1F1F",
  poisson = "#2C7FB8",
  gamma_poisson = "#F28E2B",
  correlated_log_gaussian = "#7A5195"
)

design_labels <- c(
  one_by_twelve = "1 population x 12 transitions",
  three_by_four = "3 populations x 4 transitions",
  six_by_two = "6 populations x 2 transitions"
)

model_labels <- c(
  poisson = "Poisson",
  gamma_poisson = "Gamma-Poisson",
  correlated_log_gaussian = "Correlated process"
)

draw_subset <- function(fit, n = 600L) {
  draws <- do.call(rbind, fit$chains)
  if (nrow(draws) > n) {
    draws <- draws[unique(round(seq(1, nrow(draws), length.out = n))), ]
  }
  draws
}

curve_summary <- function(draws, function_name, sizes) {
  values <- t(apply(draws, 1, function(parameter) {
    z <- (sizes - 80) / 20
    if (function_name == "survival") {
      inv_logit(
        parameter["survival_at_80"] +
          parameter["survival_slope_20"] * z
      )
    } else if (function_name == "growth") {
      pmax(
        0,
        parameter["growth_increment_80"] +
          parameter["growth_slope_20"] * z
      )
    } else if (function_name == "recruitment") {
      exp(
        parameter["log_recruitment_at_80"] +
          parameter["recruitment_slope_20"] * z
      )
    } else {
      parameter["recruit_mean_80"] +
        parameter["recruit_mean_slope_20"] * z
    }
  }))
  list(
    mean = colMeans(values),
    lower = apply(values, 2, quantile, 0.025),
    upper = apply(values, 2, quantile, 0.975)
  )
}

truth_curve <- function(function_name, sizes) {
  parameter <- inverse_ipm_truth()
  z <- (sizes - 80) / 20
  if (function_name == "survival") {
    inv_logit(
      parameter["survival_at_80"] +
        parameter["survival_slope_20"] * z
    )
  } else if (function_name == "growth") {
    pmax(
      0,
      parameter["growth_increment_80"] +
        parameter["growth_slope_20"] * z
    )
  } else if (function_name == "recruitment") {
    exp(
      parameter["log_recruitment_at_80"] +
        parameter["recruitment_slope_20"] * z
    )
  } else {
    parameter["recruit_mean_80"] +
      parameter["recruit_mean_slope_20"] * z
  }
}

draw_ribbon <- function(x, summary, color, alpha = 0.20) {
  polygon(
    c(x, rev(x)),
    c(summary$lower, rev(summary$upper)),
    border = NA,
    col = adjustcolor(color, alpha.f = alpha)
  )
  lines(x, summary$mean, col = color, lwd = 2)
}

# Figure 1: annual profile fit for the converged one-long-series Poisson model.
counts_by_design <- readRDS("results/replicate_design_counts.rds")
fit <- readRDS("results/one_by_twelve_poisson_bin_integrated_fit.rds")
counts <- counts_by_design$one_by_twelve[[1L]]
draws <- draw_subset(fit, 500L)
mesh <- fit$mesh
years_to_plot <- c(1L, 4L, 8L, 12L)

profile_predictions <- lapply(seq_len(nrow(draws)), function(index) {
  expected_profile_counts(
    setNames(draws[index, ], colnames(draws)),
    counts[1L, ],
    mesh,
    nrow(counts) - 1L
  )
})

png(
  "artifacts/profile_fit_one_by_twelve.png",
  width = 1800,
  height = 1400,
  res = 180
)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))
for (year in years_to_plot) {
  predicted <- t(vapply(
    profile_predictions,
    function(item) item[year + 1L, ],
    numeric(length(mesh))
  ))
  prediction_mean <- colMeans(predicted)
  prediction_lower <- apply(predicted, 2, quantile, 0.025)
  prediction_upper <- apply(predicted, 2, quantile, 0.975)

  plot(
    mesh,
    counts[year + 1L, ],
    type = "h",
    lwd = 3,
    col = adjustcolor(colors["truth"], alpha.f = 0.65),
    xlab = "Fish size (mm)",
    ylab = "Count",
    main = paste("Year", year),
    ylim = range(c(counts[year + 1L, ], prediction_upper))
  )
  polygon(
    c(mesh, rev(mesh)),
    c(prediction_lower, rev(prediction_upper)),
    border = NA,
    col = adjustcolor(colors["poisson"], alpha.f = 0.22)
  )
  lines(mesh, prediction_mean, col = colors["poisson"], lwd = 2)
  points(mesh, counts[year + 1L, ], pch = 16, cex = 0.6)
}
mtext(
  "Observed annual profiles and posterior IPM projections",
  outer = TRUE,
  line = 1.5,
  cex = 1.25,
  font = 2
)
dev.off()

# Figure 2: vital-rate curve recovery across the three equal-budget designs.
fits <- readRDS("results/replicate_identifiability_fits.rds")
sizes <- seq(35, 145, length.out = 111)
functions_to_plot <- c("survival", "growth", "recruitment", "recruit_mean")
function_labels <- c(
  survival = "Survival probability",
  growth = "Expected growth increment (mm)",
  recruitment = "Expected realized recruits",
  recruit_mean = "Mean recruit size (mm)"
)
designs <- c("one_by_twelve", "three_by_four", "six_by_two")

png(
  "artifacts/vital_rate_recovery_by_design.png",
  width = 2400,
  height = 2100,
  res = 200
)
par(
  mfrow = c(3, 4),
  mar = c(3.2, 2.2, 4.2, 1),
  oma = c(4, 1, 3, 1)
)
for (design in designs) {
  for (function_name in functions_to_plot) {
    fit <- fits[[paste(design, "poisson", sep = "__")]]
    summary <- curve_summary(draw_subset(fit), function_name, sizes)
    truth <- truth_curve(function_name, sizes)
    y_range <- range(c(summary$lower, summary$upper, truth))

    plot(
      sizes,
      truth,
      type = "n",
      xlab = "Size (mm)",
      ylab = "",
      main = paste(
        design_labels[design],
        function_labels[function_name],
        sep = "\n"
      ),
      cex.main = 0.76,
      ylim = y_range
    )
    draw_ribbon(sizes, summary, colors["poisson"])
    lines(sizes, truth, col = colors["truth"], lwd = 2.5, lty = 2)
  }
}
mtext(
  "Vital-rate recovery from profiles alone: Poisson inverse IPM",
  outer = TRUE,
  line = 1,
  cex = 1.35,
  font = 2
)
mtext(
  "Blue: posterior mean and 95% interval. Dashed black: individual-level truth.",
  outer = TRUE,
  side = 1,
  line = 1.5,
  cex = 0.85
)
dev.off()

# Figure 3: posterior confounding ridge and induced curve uncertainty.
fit_short <- fits[["six_by_two__poisson"]]
draws_short <- draw_subset(fit_short, 2000L)
truth <- inverse_ipm_truth()

png(
  "artifacts/recruitment_confounding.png",
  width = 1900,
  height = 900,
  res = 180
)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))
plot(
  draws_short[, "log_recruitment_at_80"],
  draws_short[, "recruitment_slope_20"],
  pch = 16,
  cex = 0.55,
  col = adjustcolor(colors["poisson"], alpha.f = 0.30),
  xlab = "Log recruitment at 80 mm",
  ylab = "Recruitment slope per 20 mm",
  main = "Coefficient ridge: 6 populations x 2 transitions"
)
points(
  truth["log_recruitment_at_80"],
  truth["recruitment_slope_20"],
  pch = 4,
  lwd = 3,
  cex = 1.5,
  col = colors["truth"]
)
legend(
  "topright",
  bty = "n",
  legend = c("Posterior draws", "Truth"),
  pch = c(16, 4),
  col = c(colors["poisson"], colors["truth"])
)

summary <- curve_summary(draws_short, "recruitment", sizes)
truth_recruitment <- truth_curve("recruitment", sizes)
plot(
  sizes,
  truth_recruitment,
  type = "n",
  xlab = "Size (mm)",
  ylab = "Expected realized recruits",
  main = "The ridge's biological consequence",
  ylim = range(c(summary$lower, summary$upper, truth_recruitment))
)
draw_ribbon(sizes, summary, colors["poisson"], alpha = 0.25)
lines(sizes, truth_recruitment, col = colors["truth"], lwd = 2.5, lty = 2)
dev.off()

# Figure 4: exploratory design and process-model comparison.
aggregate_results <- read.csv("results/replicate_identifiability_aggregate.csv")
aggregate_results$design <- factor(
  aggregate_results$design,
  levels = c("one_by_twelve", "three_by_four", "six_by_two")
)
aggregate_results$process_model <- factor(
  aggregate_results$process_model,
  levels = c("poisson", "gamma_poisson", "correlated_log_gaussian")
)
metric_matrix <- rbind(
  parameter_coverage = aggregate_results$coverage,
  derived_coverage = aggregate_results$derived_coverage
)

png(
  "artifacts/exploratory_design_comparison.png",
  width = 2000,
  height = 1100,
  res = 180
)
layout(
  matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
  heights = c(9, 1)
)
par(mar = c(5, 4, 4, 1), oma = c(0, 0, 4, 0))
for (metric in c("coverage", "derived_coverage")) {
  value <- matrix(
    NA_real_,
    nrow = 3,
    ncol = 3,
    dimnames = list(levels(aggregate_results$process_model), levels(
      aggregate_results$design
    ))
  )
  warning_matrix <- value
  for (index in seq_len(nrow(aggregate_results))) {
    value[
      as.character(aggregate_results$process_model[index]),
      as.character(aggregate_results$design[index])
    ] <- aggregate_results[[metric]][index]
    warning_matrix[
      as.character(aggregate_results$process_model[index]),
      as.character(aggregate_results$design[index])
    ] <- aggregate_results$max_rhat[index] > 1.05
  }
  positions <- barplot(
    value,
    beside = TRUE,
    col = colors[rownames(value)],
    ylim = c(0, 1.08),
    ylab = "Coverage proportion",
    names.arg = c("1 x 12", "3 x 4", "6 x 2"),
    las = 1,
    main = if (metric == "coverage") {
      "Vital-rate parameter coverage"
    } else {
      "Derived-quantity coverage"
    }
  )
  abline(h = 0.95, lty = 3, col = "gray45")
  for (row in seq_len(nrow(value))) {
    for (column in seq_len(ncol(value))) {
      if (warning_matrix[row, column]) {
        text(
          positions[row, column],
          value[row, column] + 0.045,
          labels = "*",
          cex = 1.6,
          font = 2
        )
      }
    }
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  horiz = TRUE,
  bty = "n",
  fill = colors[levels(aggregate_results$process_model)],
  legend = c("Poisson", "Gamma-Poisson", "Correlated process")
)
mtext(
  "Exploratory single-simulation comparison; * indicates maximum R-hat > 1.05",
  outer = TRUE,
  line = 1.5,
  cex = 1.15,
  font = 2
)
dev.off()

cat("Created diagnostic plots in artifacts/.\n")

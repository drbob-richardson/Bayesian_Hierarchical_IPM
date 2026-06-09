source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")

sizes <- seq(35, 145, by = 2)
z <- (sizes - 80) / 20
regime_specs <- list(
  diffuse = list(regime = "weak", multiplier = 2),
  moderate = list(regime = "all_vital_rates", multiplier = 2),
  strong = list(regime = "all_vital_rates", multiplier = 0.5),
  strong_biased = list(regime = "all_vital_rates_biased", multiplier = 0.5)
)
labels <- c(
  diffuse = "Diffuse",
  moderate = "Moderate external",
  strong = "Strong external",
  strong_biased = "Strong biased external"
)
colors <- c(
  diffuse = "#BDBDBD",
  moderate = "#4E79A7",
  strong = "#59A14F",
  strong_biased = "#E15759"
)

set.seed(20260608)
curves <- lapply(regime_specs, function(spec) {
  prior <- make_external_information_prior(
    spec$regime,
    "poisson",
    spec$multiplier
  )
  parameter_names <- names(inverse_ipm_truth())
  draws <- sapply(parameter_names, function(parameter) {
    rnorm(3000, prior$mean[parameter], prior$sd[parameter])
  })

  linear_curve <- function(intercept, slope) {
    outer(draws[, intercept], rep(1, length(sizes))) +
      outer(draws[, slope], z)
  }
  growth <- linear_curve("growth_increment_80", "growth_slope_20")
  growth[growth < 0] <- 0
  list(
    survival = inv_logit(linear_curve("survival_at_80", "survival_slope_20")),
    growth = growth,
    recruitment = exp(linear_curve(
      "log_recruitment_at_80",
      "recruitment_slope_20"
    )),
    recruit_mean = linear_curve("recruit_mean_80", "recruit_mean_slope_20")
  )
})

summarize_curve <- function(curve) {
  apply(curve, 2, quantile, probs = c(0.05, 0.5, 0.95))
}

dir.create("artifacts/prior_information", recursive = TRUE, showWarnings = FALSE)
png(
  "artifacts/prior_information/prior_predictive_vital_rates.png",
  width = 1900,
  height = 1500,
  res = 180
)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1))
settings <- list(
  survival = list(title = "Survival", ylim = c(0, 1), log = ""),
  growth = list(title = "Growth increment", ylim = c(0, 30), log = ""),
  recruitment = list(
    title = "Recruitment intensity",
    ylim = c(0.005, 20),
    log = "y"
  ),
  recruit_mean = list(
    title = "Recruit mean size",
    ylim = c(15, 80),
    log = ""
  )
)
for (quantity in names(settings)) {
  setting <- settings[[quantity]]
  plot(
    NA,
    xlim = range(sizes),
    ylim = setting$ylim,
    log = setting$log,
    xlab = "Parent size (mm)",
    ylab = setting$title,
    main = paste("Prior predictive", tolower(setting$title))
  )
  for (regime in names(curves)) {
    summary <- summarize_curve(curves[[regime]][[quantity]])
    if (regime != "strong_biased") {
      polygon(
        c(sizes, rev(sizes)),
        c(summary[1, ], rev(summary[3, ])),
        border = NA,
        col = adjustcolor(colors[[regime]], alpha.f = 0.12)
      )
    }
    lines(
      sizes,
      summary[2, ],
      lwd = 2,
      lty = ifelse(regime == "strong_biased", 2, 1),
      col = colors[[regime]]
    )
  }
  if (quantity == "recruit_mean") {
    legend(
      "bottomleft",
      bty = "n",
      lwd = 2,
      lty = c(1, 1, 1, 2),
      col = colors,
      legend = labels
    )
  }
}
dev.off()

cat("Created prior-predictive vital-rate plot.\n")

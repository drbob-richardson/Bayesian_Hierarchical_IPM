source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")

results_directory <- Sys.getenv(
  "STATE_SPACE_PLOT_RESULTS",
  "results/state_space_ipm_pilot"
)
artifact_directory <- Sys.getenv(
  "STATE_SPACE_PLOT_ARTIFACTS",
  "artifacts/state_space_ipm_pilot"
)
dir.create(artifact_directory, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path(results_directory, "multichain_fit.rds"))
draws <- do.call(rbind, lapply(fit$chains, `[[`, "draws"))
sizes <- seq(35, 145, by = 1)
z <- (sizes - 80) / 20

curve_draws <- list(
  survival = t(apply(draws, 1L, function(parameter) {
    plogis(
      parameter["survival_at_80"] +
        parameter["survival_slope_20"] * z
    )
  })),
  growth = t(apply(draws, 1L, function(parameter) {
    pmax(
      0,
      parameter["growth_increment_80"] +
        parameter["growth_slope_20"] * z
    )
  })),
  recruitment = t(apply(draws, 1L, function(parameter) {
    exp(
      parameter["log_recruitment_at_80"] +
        parameter["recruitment_slope_20"] * z
    )
  })),
  recruit_mean = t(apply(draws, 1L, function(parameter) {
    parameter["recruit_mean_80"] +
      parameter["recruit_mean_slope_20"] * z
  }))
)

truth <- inverse_ipm_truth()
truth_curves <- list(
  survival = plogis(
    truth["survival_at_80"] + truth["survival_slope_20"] * z
  ),
  growth = pmax(
    0,
    truth["growth_increment_80"] + truth["growth_slope_20"] * z
  ),
  recruitment = exp(
    truth["log_recruitment_at_80"] + truth["recruitment_slope_20"] * z
  ),
  recruit_mean = truth["recruit_mean_80"] +
    truth["recruit_mean_slope_20"] * z
)
titles <- c(
  survival = "Survival probability",
  growth = "Expected growth increment",
  recruitment = "Expected recruitment",
  recruit_mean = "Expected recruit size"
)
ylabels <- c(
  survival = "Probability",
  growth = "Increment (mm)",
  recruitment = "Recruits per individual",
  recruit_mean = "Recruit size (mm)"
)

png(
  file.path(artifact_directory, "vital_rate_curve_recovery.png"),
  width = 2100,
  height = 1700,
  res = 180
)
par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3, 1), oma = c(2.5, 0, 0, 0))
for (name in names(curve_draws)) {
  intervals <- apply(
    curve_draws[[name]],
    2L,
    quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  plot(
    sizes,
    truth_curves[[name]],
    type = "n",
    ylim = range(intervals, truth_curves[[name]]),
    xlab = "Size (mm)",
    ylab = ylabels[name],
    main = titles[name]
  )
  polygon(
    c(sizes, rev(sizes)),
    c(intervals[1L, ], rev(intervals[3L, ])),
    col = adjustcolor("#4C78A8", alpha.f = 0.25),
    border = NA
  )
  lines(sizes, intervals[2L, ], col = "#4C78A8", lwd = 3)
  lines(sizes, truth_curves[[name]], col = "#E45756", lwd = 3, lty = 2)
}
mtext(
  "Blue: posterior median and 95% interval; red dashed: simulation truth",
  side = 1,
  outer = TRUE,
  line = 1
)
dev.off()

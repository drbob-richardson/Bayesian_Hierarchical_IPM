# Robustness of the recruitment confounding to flexible vital-rate forms.
# Data are generated from the parametric truth; the recovery model uses
# natural-cubic-spline vital-rate functions (R/spline_ipm.R). If the confounding
# is structural (Proposition 1) rather than an artifact of the parametric kernel,
# the spline model should still recover survival, growth, and recruit size while
# leaving recruitment intensity poorly identified.

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/paired_benchmark.R")
source("R/monte_carlo.R")
source("R/spline_ipm.R")

mesh <- seq(25, 170, by = 5)
truth <- inverse_ipm_truth()
df <- 3L
B <- make_size_basis(mesh, df)
df1 <- ncol(B)
L <- spline_parameter_layout(df1)
results_directory <- "results/spline_confounding"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

# Prior: weak normals; intercepts at sensible baselines, spline terms shrunk to 0.
prior_mean <- numeric(4 * df1 + 3)
prior_sd <- numeric(4 * df1 + 3)
set_block <- function(idx, intercept_mean, intercept_sd, spline_sd) {
  prior_mean[idx[1]] <<- intercept_mean; prior_sd[idx[1]] <<- intercept_sd
  prior_mean[idx[-1]] <<- 0; prior_sd[idx[-1]] <<- spline_sd
}
set_block(L$survival, qlogis(0.6), 1.5, 1.5)
set_block(L$growth, 8, 6, 5)
set_block(L$recruitment, log(0.2), 1.5, 1.5)
set_block(L$recruit_mean, 45, 15, 12)
prior_mean[L$log_growth_sd] <- log(5); prior_sd[L$log_growth_sd] <- 0.7
prior_mean[L$log_recruit_sd] <- log(4); prior_sd[L$log_recruit_sd] <- 0.7
prior_mean[L$log_process_precision] <- log(50); prior_sd[L$log_process_precision] <- 2

# Moderately perturbed initial state (so survival/growth are identified while
# recruitment remains the confounded rate), large population, low noise.
initial_intensity <- 6000 * dnorm(mesh, 70, 22)
initial_profile <- round(initial_intensity / sum(initial_intensity) * 6000)
set.seed(20410001)
expected <- expected_profile_counts(truth, initial_profile, mesh, 8L)
counts <- expected
for (t in seq(2L, nrow(expected))) counts[t, ] <- rpois(ncol(expected), pmax(expected[t, ], 1e-9))
counts[1L, ] <- initial_profile

cat("Fitting spline inverse IPM (", length(prior_mean), "parameters )...\n")
fit <- fit_spline_laplace(list(counts), mesh, B, prior_mean, prior_sd, "poisson")
cat("convergence:", fit$convergence, "| PD:", fit$positive_definite, "\n")

# Posterior draws (Gaussian Laplace) -> vital-rate curves.
set.seed(20410002)
ev <- eigen((fit$covariance + t(fit$covariance)) / 2, symmetric = TRUE)
tr <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(ev$values))
draws <- sweep(matrix(rnorm(2000 * length(fit$map)), 2000) %*% t(tr), 2, fit$map, "+")

curve_band <- function(idx, link) {
  M <- apply(draws, 1, function(p) link(as.vector(B %*% p[idx])))
  t(apply(M, 1, quantile, c(.025, .5, .975)))
}
surv <- curve_band(L$survival, plogis)
grow <- curve_band(L$growth, function(x) pmax(0, x))
recr <- curve_band(L$recruitment, exp)
rmean <- curve_band(L$recruit_mean, identity)

truth_rates <- inverse_ipm_vital_rate_values(truth, mesh)
truth_curves <- list(
  survival = truth_rates$survival,
  growth = pmax(0, truth_rates$growth_mean - mesh),
  recruitment = truth_rates$recruitment,
  recruit_mean = truth_rates$recruit_mean
)

# Identifiability metric: mean posterior CV of each curve over sizes.
cv_curve <- function(idx, link) {
  M <- apply(draws, 1, function(p) link(as.vector(B %*% p[idx])))
  mean(apply(M, 1, sd) / pmax(abs(rowMeans(M)), 1e-8))
}
metrics <- data.frame(
  rate = c("survival", "growth", "recruitment", "recruit_mean"),
  mean_posterior_cv = c(cv_curve(L$survival, plogis), cv_curve(L$growth, function(x) pmax(0, x)),
                        cv_curve(L$recruitment, exp), cv_curve(L$recruit_mean, identity))
)
write.csv(metrics, file.path(results_directory, "curve_cv.csv"), row.names = FALSE)
saveRDS(list(fit = fit, draws = draws, mesh = mesh, B = B,
             curves = list(surv = surv, grow = grow, recr = recr, rmean = rmean),
             truth_curves = truth_curves, metrics = metrics),
        file.path(results_directory, "spline_fit.rds"), compress = "xz")

# Figure: recovered curves with 95% bands vs truth.
png("paper/figures/spline_confounding.png", width = 1900, height = 1450, res = 220)
par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.4, 1))
panel <- function(band, tru, ylab, main, col) {
  yl <- range(c(band, tru), na.rm = TRUE)
  plot(mesh, band[, 2], type = "n", ylim = yl, xlab = "size (mm)", ylab = ylab, main = main)
  polygon(c(mesh, rev(mesh)), c(band[, 1], rev(band[, 3])), col = adjustcolor(col, 0.2), border = NA)
  lines(mesh, band[, 2], col = col, lwd = 2.3)
  lines(mesh, tru, col = "black", lwd = 2, lty = 2)
  legend("topleft", c("spline posterior", "truth"), col = c(col, "black"),
         lwd = c(2.3, 2), lty = c(1, 2), bty = "n", cex = 0.85)
}
panel(surv, truth_curves$survival, "survival", "(a) Survival", "#1b9e77")
panel(grow, truth_curves$growth, "growth increment (mm)", "(b) Growth", "#7570b3")
panel(recr, truth_curves$recruitment, "recruitment intensity", "(c) Recruitment", "#d95f02")
panel(rmean, truth_curves$recruit_mean, "recruit mean size (mm)", "(d) Recruit size", "#e7298a")
dev.off()

cat("\n=== Spline (flexible) recovery; mean posterior CV by rate ===\n")
print(metrics, row.names = FALSE, digits = 3)
cat(sprintf("\nRecruitment CV / survival CV = %.1fx\n",
            metrics$mean_posterior_cv[3] / metrics$mean_posterior_cv[1]))
cat("Figure: paper/figures/spline_confounding.png\n")

# Kernel functionals for the West Brook brook trout fit: the model-implied
# stable size distribution and asymptotic growth rate, with posterior bands,
# compared to the observed average size profile. If the population is near its
# stable distribution (the observed profiles are near-stationary), the two
# should coincide -- validating the kernel AND making visible the Proposition-2
# condition under which recruitment is non-identified. We also check that the
# stable distribution, like lambda, is well-identified even though the vital
# rates are not (the iso-lambda / iso-w manifold).

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/kernel_functionals.R")

results_directory <- "results/west_brook_application"
profile <- readRDS(file.path(results_directory, "profile_fit.rds"))
mcmc <- readRDS(file.path(results_directory, "profile_mcmc.rds"))
mesh <- profile$mesh; pnames <- names(profile$prior$mean)

# Observed average size profile (proportion), over the 20 annual occasions.
counts <- profile$counts
obs_prop <- colMeans(counts / rowSums(counts))

# Posterior stable distribution + lambda from the (pooled) MCMC draws.
draws <- mcmc$draws
draws <- draws[round(seq(1, nrow(draws), length.out = 1200L)), , drop = FALSE]
W <- matrix(NA_real_, nrow(draws), length(mesh))
lam <- numeric(nrow(draws))
for (k in seq_len(nrow(draws))) {
  f <- inverse_ipm_functionals(setNames(draws[k, ], pnames), mesh)
  W[k, ] <- f$stable_distribution
  lam[k] <- f$lambda
}
w_band <- t(apply(W, 2, quantile, c(.025, .5, .975)))

# Identifiability: per-bin CV of the stable distribution vs the vital rates.
w_cv <- mean(apply(W, 2, sd) / pmax(colMeans(W), 1e-8))
surv_cv <- sd(plogis(draws[, "survival_at_80"])) / mean(plogis(draws[, "survival_at_80"]))

summary <- data.frame(
  quantity = c("lambda", "stable_distribution (mean per-bin CV)",
               "survival@80 (CV, for contrast)"),
  value = c(sprintf("%.3f [%.3f, %.3f]", median(lam), quantile(lam, .025), quantile(lam, .975)),
            sprintf("%.3f", w_cv), sprintf("%.3f", surv_cv))
)
write.csv(data.frame(size = mesh, observed = obs_prop,
                     stable_q025 = w_band[, 1], stable_med = w_band[, 2], stable_q975 = w_band[, 3]),
          file.path(results_directory, "stable_distribution.csv"), row.names = FALSE)

png("paper/figures/west_brook_stable.png", width = 1500, height = 1100, res = 230)
par(mar = c(4.5, 4.5, 2.5, 1))
yl <- c(0, max(c(obs_prop, w_band)))
plot(mesh, obs_prop, type = "n", ylim = yl, xlab = "size (mm)",
     ylab = "proportion", main = "West Brook: stable distribution vs observed profile")
polygon(c(mesh, rev(mesh)), c(w_band[, 1], rev(w_band[, 3])),
        col = adjustcolor("#1b9e77", 0.25), border = NA)
lines(mesh, w_band[, 2], col = "#1b9e77", lwd = 2.3)
points(mesh, obs_prop, pch = 19, col = "black", cex = 0.8)
lines(mesh, obs_prop, col = "black", lwd = 1, lty = 3)
legend("topright", c("model stable distribution (95%)", "observed average profile"),
       col = c("#1b9e77", "black"), lwd = c(2.3, 1), pch = c(NA, 19), lty = c(1, 3), bty = "n", cex = 0.85)
dev.off()

cat("=== West Brook kernel functionals ===\n")
print(summary, row.names = FALSE)
cat(sprintf("\nStable distribution mean per-bin CV = %.3f vs survival@80 CV = %.3f (%.0fx tighter)\n",
            w_cv, surv_cv, surv_cv / w_cv))
cat("Figure: paper/figures/west_brook_stable.png\n")

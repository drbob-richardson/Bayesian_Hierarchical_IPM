# Excitation experiment (test of Proposition 2). Generate matched profile data
# sets that differ ONLY in how far the initial state is perturbed from the
# stable size distribution w, from exactly stable (no excitation) to a strong
# small-size cohort pulse. Large population, low observation noise, so the result
# isolates structural identifiability. We measure whether the recruitment
# confounding (level/slope ridge) breaks and whether recruitment intensity and
# lambda are recovered as excitation increases. Residual confounding that
# survives excitation is the decomposition limit of Proposition 1.

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
library(parallel)

mesh <- seq(25, 170, by = 5)
truth <- inverse_ipm_truth()
process <- "poisson"
prior <- default_inverse_ipm_priors("weak", process)
transitions <- 6L
N0 <- 8000
alphas <- seq(0, 1, by = 0.1)        # perturbation fraction toward the cohort pulse
seeds_per <- 4L
results_directory <- "results/excitation"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

ipm <- build_centered_ipm(truth, mesh)
K <- ipm$delta * ipm$kernel
eg <- eigen(K)
lambda_true <- max(Re(eg$values))
w <- abs(Re(eg$vectors[, which.max(Re(eg$values))])); w <- w / sum(w)
truth_derived <- derived_inverse_ipm_quantities(truth, mesh)
pulse <- dnorm(mesh, mean = 45, sd = 8); pulse <- pulse / sum(pulse)   # small-size cohort pulse

run_one <- function(row) {
  alpha <- row$alpha
  mu0 <- (1 - alpha) * w + alpha * pulse; mu0 <- mu0 / sum(mu0)
  cosang <- sum(mu0 * w) / sqrt(sum(mu0^2) * sum(w^2))
  angle <- acos(pmin(1, cosang)) * 180 / pi
  init_counts <- round(N0 * mu0)
  set.seed(row$seed)
  expected <- expected_profile_counts(truth, init_counts, mesh, transitions)
  if (is.null(expected) || any(!is.finite(expected))) return(NULL)
  counts <- expected
  for (t in seq(2L, nrow(expected))) counts[t, ] <- rpois(ncol(expected), pmax(expected[t, ], 1e-9))
  counts[1L, ] <- init_counts
  fit <- tryCatch(fit_laplace_inverse_ipm(list(counts), mesh, prior, process),
                  error = function(e) NULL)
  if (is.null(fit) || !fit$positive_definite) return(NULL)
  corr <- cov2cor(fit$covariance)["log_recruitment_at_80", "recruitment_slope_20"]
  post_sd <- sqrt(diag(fit$covariance))
  contraction <- 1 - post_sd / prior$sd[names(post_sd)]
  d <- derived_inverse_ipm_quantities(fit$map, mesh)
  data.frame(
    alpha = alpha, angle = angle,
    ridge_corr = corr,
    recruit_contraction = contraction["log_recruitment_at_80"],
    survival_contraction = contraction["survival_at_80"],
    lambda_relerr = abs(d["lambda"] - truth_derived["lambda"]) / truth_derived["lambda"],
    recruit80_relerr = abs(d["recruitment_80"] - truth_derived["recruitment_80"]) /
      truth_derived["recruitment_80"]
  )
}

grid <- expand.grid(alpha = alphas, seed = 20400000 + seq_len(seeds_per))
raw <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))), run_one,
                               mc.cores = min(8L, seeds_per * 2L)))
agg <- aggregate(cbind(angle, ridge_corr, recruit_contraction, survival_contraction,
                       lambda_relerr, recruit80_relerr) ~ alpha, data = raw, FUN = mean)
write.csv(agg, file.path(results_directory, "excitation_summary.csv"), row.names = FALSE)

png("paper/figures/excitation.png", width = 1900, height = 800, res = 220)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.6, 1))
o <- order(agg$angle); a <- agg[o, ]
plot(a$angle, abs(a$ridge_corr), type = "b", pch = 19, col = "#d95f02", ylim = c(0, 1),
     xlab = "excitation: angle of initial state from stable w (deg)",
     ylab = "value", main = "(a) Confounding vs excitation")
lines(a$angle, a$recruit_contraction, type = "b", pch = 17, col = "#1b9e77")
lines(a$angle, a$survival_contraction, type = "b", pch = 15, col = "#7570b3")
legend("right", c("|cor(recruit level,slope)|", "recruitment contraction", "survival contraction"),
       col = c("#d95f02", "#1b9e77", "#7570b3"), pch = c(19, 17, 15), bty = "n", cex = 0.85)
grid(col = "grey90")
plot(a$angle, a$recruit80_relerr, type = "b", pch = 19, col = "#1b9e77",
     ylim = c(0, max(a$recruit80_relerr, a$lambda_relerr)),
     xlab = "excitation: angle of initial state from stable w (deg)",
     ylab = "relative error of posterior mode", main = "(b) Recovery vs excitation")
lines(a$angle, a$lambda_relerr, type = "b", pch = 17, col = "#e7298a")
legend("topright", c("recruitment @ 80", expression(lambda)),
       col = c("#1b9e77", "#e7298a"), pch = c(19, 17), bty = "n", cex = 0.9)
grid(col = "grey90")
dev.off()

cat("=== Excitation experiment (Proposition 2 test) ===\n")
cat(sprintf("True lambda = %.3f; stable w used as alpha=0 initial state.\n", lambda_true))
print(round(agg, 3), row.names = FALSE)
cat("\nFigure: paper/figures/excitation.png\n")

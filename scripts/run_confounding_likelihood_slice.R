# Numerical confirmation of the recruitment confounding proposition.
# Generates a large, low-noise dataset from the inverse-IPM model at truth,
# then (1) evaluates a 2-D profile log-likelihood slice over recruitment
# level and slope to expose the flat ridge, and (2) fits the full Laplace
# posterior to report the recruitment level/slope correlation and the global
# least-informed direction (smallest-eigenvalue eigenvector of the observed
# information). This makes the -0.99 ridge in the manuscript reproducible.

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

mesh <- seq(25, 160, by = 5)
truth <- inverse_ipm_truth()
process_model <- "poisson"
results_directory <- "results/confounding_slice"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

# Large, low-noise dataset so the slice reflects STRUCTURAL identifiability,
# not Monte Carlo noise: one long trajectory, full detection, big population.
transitions <- 8L
initial_intensity <- 6000 * dnorm(mesh, mean = 70, sd = 22)
initial_profile <- round(initial_intensity / sum(initial_intensity) * 6000)
set.seed(20310001)
expected <- expected_profile_counts(truth, initial_profile, mesh, transitions)
counts <- expected
for (t in seq(2L, nrow(expected))) {
  counts[t, ] <- rpois(ncol(expected), lambda = pmax(expected[t, ], 1e-9))
}
counts[1L, ] <- initial_profile
counts <- list(trajectory_1 = counts)

prior <- default_inverse_ipm_priors("weak", process_model)

log_likelihood_at <- function(parameter) {
  inverse_ipm_log_posterior(parameter, counts, mesh, prior, process_model) -
    inverse_ipm_log_prior(parameter, prior)
}

# --- (1) 2-D slice over recruitment level and slope, nuisance fixed at truth.
level_grid <- truth["log_recruitment_at_80"] + seq(-1.0, 1.0, length.out = 61)
slope_grid <- truth["recruitment_slope_20"] + seq(-1.0, 1.0, length.out = 61)
slice <- matrix(NA_real_, nrow = length(level_grid), ncol = length(slope_grid))
for (i in seq_along(level_grid)) {
  for (j in seq_along(slope_grid)) {
    parameter <- truth
    parameter["log_recruitment_at_80"] <- level_grid[i]
    parameter["recruitment_slope_20"] <- slope_grid[j]
    slice[i, j] <- log_likelihood_at(parameter)
  }
}
dimnames(slice) <- list(
  sprintf("%.4f", level_grid), sprintf("%.4f", slope_grid)
)
write.csv(slice, file.path(results_directory, "loglik_slice.csv"))

# Local curvature of the 2-D slice at truth -> flat (ridge) direction.
two_d_hessian <- function() {
  h <- 1e-3
  base <- c(truth["log_recruitment_at_80"], truth["recruitment_slope_20"])
  f <- function(a, b) {
    parameter <- truth
    parameter["log_recruitment_at_80"] <- a
    parameter["recruitment_slope_20"] <- b
    log_likelihood_at(parameter)
  }
  faa <- (f(base[1] + h, base[2]) - 2 * f(base[1], base[2]) + f(base[1] - h, base[2])) / h^2
  fbb <- (f(base[1], base[2] + h) - 2 * f(base[1], base[2]) + f(base[1], base[2] - h)) / h^2
  fab <- (f(base[1] + h, base[2] + h) - f(base[1] + h, base[2] - h) -
            f(base[1] - h, base[2] + h) + f(base[1] - h, base[2] - h)) / (4 * h^2)
  matrix(c(faa, fab, fab, fbb), 2, 2)
}
slice_hessian <- two_d_hessian()
slice_eigen <- eigen(-slice_hessian, symmetric = TRUE)
flat_direction <- slice_eigen$vectors[, which.min(slice_eigen$values)]
names(flat_direction) <- c("log_recruitment_at_80", "recruitment_slope_20")
implied_slice_correlation <- -slice_hessian[1, 2] /
  sqrt(slice_hessian[1, 1] * slice_hessian[2, 2])

# --- (2) Full Laplace fit: posterior correlation and global flat direction.
fit <- fit_laplace_inverse_ipm(counts, mesh, prior, process_model)
correlation <- cov2cor(fit$covariance)
recruit_corr <- correlation["log_recruitment_at_80", "recruitment_slope_20"]
information <- solve(fit$covariance)
info_eigen <- eigen((information + t(information)) / 2, symmetric = TRUE)
global_flat_direction <- info_eigen$vectors[, which.min(info_eigen$values)]
names(global_flat_direction) <- rownames(fit$covariance)

saveRDS(
  list(
    counts = counts, truth = truth, mesh = mesh,
    slice = slice, level_grid = level_grid, slope_grid = slope_grid,
    slice_hessian = slice_hessian, flat_direction = flat_direction,
    implied_slice_correlation = implied_slice_correlation,
    fit = fit, correlation = correlation,
    recruit_level_slope_correlation = recruit_corr,
    global_flat_direction = global_flat_direction,
    information_eigenvalues = info_eigen$values
  ),
  file.path(results_directory, "confounding_slice.rds"), compress = "xz"
)

cat("=== Recruitment confounding, numerical confirmation ===\n")
cat(sprintf(
  "Posterior corr(log_recruitment_at_80, recruitment_slope_20): %.4f\n",
  recruit_corr
))
cat(sprintf(
  "2-D slice implied correlation along ridge: %.4f\n", implied_slice_correlation
))
cat("Flat direction in (level, slope) plane (unit vector):\n")
print(round(flat_direction, 4))
cat("\nGlobal least-informed direction (largest-magnitude loadings):\n")
print(round(sort(global_flat_direction[order(-abs(global_flat_direction))][1:4], decreasing = TRUE), 4))
cat(sprintf(
  "\nInformation eigenvalue ratio (max/min): %.3e\n",
  max(info_eigen$values) / min(info_eigen$values)
))

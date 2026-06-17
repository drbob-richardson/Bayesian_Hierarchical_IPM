# Full MCMC for the West Brook profile-only inverse IPM, to replace the Laplace
# Gaussian approximation on the confounded directions (recruitment level/slope
# and lambda). SBC showed the Laplace posterior is overconfident along the
# recruitment ridge under diffuse priors; this quantifies the honest intervals.

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")
library(parallel)

chains_n <- as.integer(Sys.getenv("WB_MCMC_CHAINS", "4"))
iterations <- as.integer(Sys.getenv("WB_MCMC_ITERATIONS", "40000"))
warmup <- as.integer(Sys.getenv("WB_MCMC_WARMUP", "20000"))
results_directory <- "results/west_brook_application"

profile <- readRDS(file.path(results_directory, "profile_fit.rds"))
counts <- profile$counts; mesh <- profile$mesh; prior <- profile$prior
process_model <- "gamma_poisson"
pnames <- names(prior$mean)
map <- profile$fit$map[pnames]
cov <- profile$fit$covariance[pnames, pnames]
sds <- sqrt(diag(cov))

log_post <- function(parameter) {
  inverse_ipm_log_posterior(setNames(parameter, pnames), counts, mesh, prior, process_model)
}

run_one <- function(chain_index) {
  # overdispersed start from the Laplace mode, kept finite
  repeat {
    init <- map + rnorm(length(map), 0, 0.5 * sds)
    names(init) <- pnames
    if (is.finite(log_post(init))) break
  }
  run_adaptive_metropolis(log_post, init, cov, iterations, warmup,
                          seed = 20360000 + 101 * chain_index, target_acceptance = 0.234)
}

cat(sprintf("West Brook MCMC: %d chains x %d iters (warmup %d).\n", chains_n, iterations, warmup))
cores <- max(1L, min(chains_n, 8L))
fits <- mclapply(seq_len(chains_n), run_one, mc.cores = cores, mc.preschedule = FALSE)

chain_draws <- lapply(fits, `[[`, "draws")
rhat <- setNames(split_rhat(chain_draws), pnames)
ess <- setNames(effective_sample_size(chain_draws), pnames)
draws <- do.call(rbind, chain_draws)

# Derived quantities on a thinned set (lambda eigen-decomposition is the cost).
thin <- draws[round(seq(1, nrow(draws), length.out = 3000L)), , drop = FALSE]
derived <- t(apply(thin, 1, function(p) derived_inverse_ipm_quantities(setNames(p, pnames), mesh)))

q <- function(x) quantile(x, c(.025, .5, .975))
# Compare Laplace vs MCMC for the confounded directions and an identified one.
laplace_draws <- profile$draws
lap_der <- t(apply(laplace_draws[round(seq(1, nrow(laplace_draws), length.out = 3000L)), ], 1,
                   function(p) derived_inverse_ipm_quantities(setNames(p, pnames), mesh)))

compare <- function(name, mcmc_vec, lap_vec) {
  data.frame(quantity = name,
             laplace_lo = q(lap_vec)[1], laplace_med = q(lap_vec)[2], laplace_hi = q(lap_vec)[3],
             mcmc_lo = q(mcmc_vec)[1], mcmc_med = q(mcmc_vec)[2], mcmc_hi = q(mcmc_vec)[3],
             width_ratio = diff(q(mcmc_vec)[c(1,3)]) / diff(q(lap_vec)[c(1,3)]))
}
rows <- rbind(
  compare("lambda", derived[, "lambda"], lap_der[, "lambda"]),
  compare("recruitment_80", derived[, "recruitment_80"], lap_der[, "recruitment_80"]),
  compare("recruitment_110", derived[, "recruitment_110"], lap_der[, "recruitment_110"]),
  compare("log_recruitment_at_80", thin[, "log_recruitment_at_80"], laplace_draws[, "log_recruitment_at_80"]),
  compare("survival_80", derived[, "survival_80"], lap_der[, "survival_80"]),
  compare("growth_increment_80", derived[, "growth_increment_80"], lap_der[, "growth_increment_80"])
)

saveRDS(list(fits = fits, draws = draws, derived = derived, thin = thin,
             rhat = rhat, ess = ess, comparison = rows,
             acceptance = sapply(fits, `[[`, "acceptance")),
        file.path(results_directory, "profile_mcmc.rds"), compress = "xz")
write.csv(rows, file.path(results_directory, "laplace_vs_mcmc.csv"), row.names = FALSE)

cat(sprintf("Acceptance: %s\n", paste(round(sapply(fits, `[[`, "acceptance"), 3), collapse = ", ")))
cat(sprintf("Max Rhat: %.3f | Min ESS: %.0f\n", max(rhat), min(ess)))
cat("Laplace vs MCMC (median [95%]); width_ratio = MCMC/Laplace interval width:\n")
for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]
  cat(sprintf("  %-22s Laplace %.3f [%.3f, %.3f] | MCMC %.3f [%.3f, %.3f] | x%.2f\n",
              r$quantity, r$laplace_med, r$laplace_lo, r$laplace_hi,
              r$mcmc_med, r$mcmc_lo, r$mcmc_hi, r$width_ratio))
}

# Simulation-based calibration (SBC) of the Laplace-approximate inverse-IPM
# posterior. Data are generated from the recovery model's OWN data-generating
# process, so departures from rank uniformity isolate computational /
# approximation error (Laplace + Gaussian sampling) rather than model
# misspecification.

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

library(parallel)

regime <- Sys.getenv("SBC_REGIME", "informative")
process_model <- Sys.getenv("SBC_PROCESS", "gamma_poisson")
iterations <- as.integer(Sys.getenv("SBC_ITERATIONS", "200"))
posterior_draws <- as.integer(Sys.getenv("SBC_DRAWS", "99"))
n_trajectories <- as.integer(Sys.getenv("SBC_TRAJECTORIES", "3"))
transitions <- as.integer(Sys.getenv("SBC_TRANSITIONS", "4"))
core_cap <- as.integer(Sys.getenv("SBC_CORES", "12"))

mesh <- seq(25, 160, by = 5)
results_directory <- file.path("results", "sbc", paste0(regime, "_", process_model))
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

prior <- default_inverse_ipm_priors(regime, process_model)
parameter_names <- names(prior$mean)

# Fixed, conditioned initial profile (broad bump centered near 70 mm, ~1500 fish).
initial_intensity <- 1500 * dnorm(mesh, mean = 70, sd = 22)
initial_profile <- round(initial_intensity / sum(initial_intensity) * 1500)

draw_prior_parameter <- function(seed) {
  set.seed(seed)
  setNames(rnorm(length(prior$mean), prior$mean, prior$sd), parameter_names)
}

simulate_from_model <- function(parameter, seed) {
  set.seed(seed)
  expected <- expected_profile_counts(
    parameter, initial_profile, mesh, transitions
  )
  if (is.null(expected) || any(!is.finite(expected)) || any(expected < 0)) {
    return(NULL)
  }
  size <- if (process_model == "gamma_poisson") {
    exp(parameter["log_process_precision"])
  } else {
    NA_real_
  }
  lapply(seq_len(n_trajectories), function(j) {
    obs <- expected
    for (t in seq(2L, nrow(expected))) {
      mu <- pmax(expected[t, ], 1e-9)
      obs[t, ] <- if (process_model == "gamma_poisson") {
        rnbinom(length(mu), mu = mu, size = size)
      } else {
        rpois(length(mu), lambda = mu)
      }
    }
    obs[1L, ] <- initial_profile
    obs
  })
}

sbc_iteration <- function(index) {
  truth <- draw_prior_parameter(20290000 + index)
  counts <- simulate_from_model(truth, 20291000 + index)
  if (is.null(counts)) {
    return(list(index = index, ok = FALSE, reason = "degenerate_projection"))
  }
  out <- tryCatch({
    fit <- fit_laplace_inverse_ipm(counts, mesh, prior, process_model)
    draws <- sample_laplace_draws(fit, draws = posterior_draws, seed = 20292000 + index)
    if (nrow(draws) < posterior_draws) {
      return(list(index = index, ok = FALSE, reason = "insufficient_draws"))
    }
    ranks <- vapply(parameter_names, function(p) {
      sum(draws[, p] < truth[p])
    }, integer(1))
    list(
      index = index, ok = TRUE,
      positive_definite = fit$positive_definite,
      ranks = ranks
    )
  }, error = function(e) {
    list(index = index, ok = FALSE, reason = conditionMessage(e))
  })
  out
}

available_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN", intern = TRUE
)))
if (!is.finite(available_cores)) available_cores <- detectCores(logical = TRUE)
cores <- max(1L, min(core_cap, available_cores - 1L))

cat(sprintf(
  "SBC: regime=%s process=%s iterations=%d draws=%d on %d cores.\n",
  regime, process_model, iterations, posterior_draws, cores
))

results <- mclapply(
  seq_len(iterations), sbc_iteration,
  mc.cores = cores, mc.preschedule = FALSE
)

ok <- vapply(results, function(r) isTRUE(r$ok), logical(1))
n_ok <- sum(ok)
rank_matrix <- t(vapply(results[ok], function(r) r$ranks, integer(length(parameter_names))))
colnames(rank_matrix) <- parameter_names

# Per-parameter uniformity test: chi-square over (posterior_draws + 1) rank bins.
expected_per_bin <- n_ok / (posterior_draws + 1L)
uniformity <- data.frame(
  parameter = parameter_names,
  mean_rank = colMeans(rank_matrix),
  expected_mean_rank = posterior_draws / 2,
  chisq = vapply(parameter_names, function(p) {
    counts_per_bin <- tabulate(rank_matrix[, p] + 1L, nbins = posterior_draws + 1L)
    sum((counts_per_bin - expected_per_bin)^2 / expected_per_bin)
  }, numeric(1)),
  df = posterior_draws,
  row.names = NULL
)
uniformity$chisq_p <- pchisq(uniformity$chisq, df = uniformity$df, lower.tail = FALSE)

saveRDS(
  list(results = results, rank_matrix = rank_matrix, prior = prior,
       posterior_draws = posterior_draws, regime = regime,
       process_model = process_model),
  file.path(results_directory, "sbc_fit.rds"), compress = "xz"
)
write.csv(rank_matrix, file.path(results_directory, "ranks.csv"), row.names = FALSE)
write.csv(uniformity, file.path(results_directory, "uniformity.csv"), row.names = FALSE)

failures <- vapply(results[!ok], function(r) r$reason, character(1))
cat(sprintf(
  "SBC complete: %d/%d valid fits (%d failures).\n",
  n_ok, iterations, sum(!ok)
))
if (length(failures)) {
  cat("Failure reasons:\n")
  print(table(failures))
}
cat("Worst-calibrated parameters (smallest uniformity p-value):\n")
print(head(uniformity[order(uniformity$chisq_p), c("parameter", "mean_rank", "chisq_p")], 5))

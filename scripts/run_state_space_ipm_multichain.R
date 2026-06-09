source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/state_space_ipm.R")

library(parallel)

results_directory <- "results/state_space_ipm_pilot"
pilot <- readRDS(file.path(results_directory, "fit.rds"))
chains <- as.integer(Sys.getenv("STATE_SPACE_CHAINS", "4"))
particles <- as.integer(Sys.getenv("STATE_SPACE_PARTICLES", "1000"))
iterations <- as.integer(Sys.getenv("STATE_SPACE_ITERATIONS", "3000"))
warmup <- as.integer(Sys.getenv(
  "STATE_SPACE_WARMUP",
  as.character(floor(iterations / 2))
))

demographic_names <- names(inverse_ipm_truth())
initial <- c(
  pilot$initializer$map[demographic_names],
  log_state_fano = 0
)
proposal_covariance <- matrix(
  0,
  nrow = length(initial),
  ncol = length(initial),
  dimnames = list(names(initial), names(initial))
)
proposal_covariance[demographic_names, demographic_names] <-
  pilot$initializer$covariance[
    demographic_names,
    demographic_names,
    drop = FALSE
  ]
proposal_covariance["log_state_fano", "log_state_fano"] <- 0.5^2

available_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN",
  intern = TRUE
)))
if (!is.finite(available_cores)) {
  available_cores <- detectCores(logical = TRUE)
}
if (!is.finite(available_cores)) {
  available_cores <- 1L
}
cores <- max(1L, min(chains, available_cores))

cat(
  "Running", chains, "state-space PMMH chains with", particles,
  "particles on", cores, "cores.\n"
)
chain_results <- mclapply(seq_len(chains), function(chain_index) {
  run_state_space_pmmh(
    counts = pilot$counts,
    mesh = pilot$mesh,
    prior = pilot$prior,
    initial = initial,
    proposal_covariance = proposal_covariance,
    detection_probability = pilot$detection_probability,
    particles = particles,
    iterations = iterations,
    warmup = warmup,
    seed = 20282000 + 1000L * chain_index
  )
}, mc.cores = cores, mc.preschedule = FALSE)

retained_chains <- lapply(chain_results, `[[`, "draws")
rhat <- setNames(split_rhat(retained_chains), names(initial))
ess <- setNames(effective_sample_size(retained_chains), names(initial))
combined <- do.call(rbind, retained_chains)

truth <- c(inverse_ipm_truth(), log_state_fano = NA_real_)
summary <- summarize_draws_against_truth(combined, truth)
summary$prior_mean <- pilot$prior$mean[summary$quantity]
summary$prior_sd <- pilot$prior$sd[summary$quantity]
summary$contraction <- 1 - summary$posterior_sd / summary$prior_sd
summary$rhat <- rhat[summary$quantity]
summary$ess <- ess[summary$quantity]

derived_draws <- t(apply(combined[, demographic_names, drop = FALSE], 1L, function(
  parameter
) {
  derived_inverse_ipm_quantities(parameter, pilot$mesh)
}))
derived_truth <- derived_inverse_ipm_quantities(
  inverse_ipm_truth(),
  pilot$mesh
)
derived_summary <- summarize_draws_against_truth(derived_draws, derived_truth)

saveRDS(
  list(
    chains = chain_results,
    rhat = rhat,
    ess = ess,
    prior = pilot$prior,
    mesh = pilot$mesh,
    counts = pilot$counts,
    detection_probability = pilot$detection_probability
  ),
  file.path(results_directory, "multichain_fit.rds"),
  compress = "xz"
)
write.csv(
  summary,
  file.path(results_directory, "multichain_parameter_summary.csv"),
  row.names = FALSE
)
write.csv(
  derived_summary,
  file.path(results_directory, "multichain_derived_summary.csv"),
  row.names = FALSE
)

cat(
  "Acceptance:",
  paste(round(vapply(chain_results, `[[`, numeric(1), "acceptance"), 3),
    collapse = ", "
  ),
  "\nMax Rhat:", round(max(rhat), 3),
  "\nMin ESS:", round(min(ess), 1),
  "\nVital-rate coverage:",
  sum(summary$covers_truth[summary$quantity %in% demographic_names]),
  "of", length(demographic_names),
  "\n"
)

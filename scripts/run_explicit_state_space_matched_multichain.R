source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/state_space_ipm.R")
source("R/explicit_state_space_ipm.R")

library(parallel)

results_directory <- "results/explicit_state_space_matched_pilot"
pilot <- readRDS(file.path(results_directory, "fit.rds"))
chains <- as.integer(Sys.getenv("EXPLICIT_MATCHED_CHAINS", "4"))
particles <- as.integer(Sys.getenv("EXPLICIT_MATCHED_PARTICLES", "500"))
iterations <- as.integer(Sys.getenv("EXPLICIT_MATCHED_ITERATIONS", "2000"))
warmup <- as.integer(Sys.getenv(
  "EXPLICIT_MATCHED_WARMUP",
  as.character(floor(iterations / 2))
))

parameter_names <- names(inverse_ipm_truth())
initial <- pilot$initializer$map[parameter_names]
proposal_covariance <- pilot$initializer$covariance[
  parameter_names,
  parameter_names,
  drop = FALSE
]
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
  "Running", chains, "matched explicit-state PMMH chains with",
  particles, "particles on", cores, "cores.\n"
)
chain_results <- mclapply(seq_len(chains), function(chain_index) {
  run_explicit_state_space_pmmh(
    pilot$counts,
    pilot$mesh,
    pilot$prior,
    initial,
    proposal_covariance,
    pilot$detection_probability,
    particles,
    iterations,
    warmup,
    seed = 20285000 + 1000L * chain_index
  )
}, mc.cores = cores, mc.preschedule = FALSE)

retained_chains <- lapply(chain_results, `[[`, "draws")
rhat <- setNames(split_rhat(retained_chains), parameter_names)
ess <- setNames(effective_sample_size(retained_chains), parameter_names)
combined <- do.call(rbind, retained_chains)
summary <- summarize_draws_against_truth(combined, inverse_ipm_truth())
summary$prior_mean <- pilot$prior$mean[summary$quantity]
summary$prior_sd <- pilot$prior$sd[summary$quantity]
summary$contraction <- 1 - summary$posterior_sd / summary$prior_sd
summary$rhat <- rhat[summary$quantity]
summary$ess <- ess[summary$quantity]
derived_draws <- t(apply(combined, 1L, function(parameter) {
  derived_inverse_ipm_quantities(parameter, pilot$mesh)
}))
derived_summary <- summarize_draws_against_truth(
  derived_draws,
  derived_inverse_ipm_quantities(inverse_ipm_truth(), pilot$mesh)
)

saveRDS(
  list(
    chains = chain_results,
    rhat = rhat,
    ess = ess,
    prior = pilot$prior,
    counts = pilot$counts,
    mesh = pilot$mesh,
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
  "\nVital-rate coverage:", sum(summary$covers_truth), "of", nrow(summary),
  "\n"
)


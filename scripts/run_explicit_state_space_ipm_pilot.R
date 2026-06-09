source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")
source("R/state_space_ipm.R")
source("R/explicit_state_space_ipm.R")

particles <- as.integer(Sys.getenv("EXPLICIT_STATE_PARTICLES", "500"))
iterations <- as.integer(Sys.getenv("EXPLICIT_STATE_ITERATIONS", "2000"))
warmup <- as.integer(Sys.getenv(
  "EXPLICIT_STATE_WARMUP",
  as.character(floor(iterations / 2))
))
populations <- as.integer(Sys.getenv("EXPLICIT_STATE_POPULATIONS", "6"))
results_directory <- "results/explicit_state_space_ipm_pilot"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

gamma_pilot <- readRDS("results/state_space_ipm_pilot/fit.rds")
counts <- gamma_pilot$counts[seq_len(populations)]
mesh <- gamma_pilot$mesh
detection_probability <- gamma_pilot$detection_probability
prior <- default_inverse_ipm_priors("weak", "poisson")
parameter_names <- names(inverse_ipm_truth())
initial <- gamma_pilot$initializer$map[parameter_names]
proposal_covariance <- gamma_pilot$initializer$covariance[
  parameter_names,
  parameter_names,
  drop = FALSE
]

cat("Checking explicit particle-likelihood variance.\n")
likelihood_replicates <- vapply(seq_len(20L), function(index) {
  explicit_state_space_particle_filter(
    initial,
    counts,
    mesh,
    detection_probability,
    particles,
    seed = 20283000 + index
  )$log_likelihood
}, numeric(1))
write.csv(
  data.frame(log_likelihood = likelihood_replicates),
  file.path(results_directory, "likelihood_variance.csv"),
  row.names = FALSE
)
cat(
  "Explicit particle log-likelihood SD:",
  round(sd(likelihood_replicates), 3),
  "\n"
)

cat(
  "Running explicit PMMH with", populations, "populations,",
  particles, "particles, and", iterations, "iterations.\n"
)
fit <- run_explicit_state_space_pmmh(
  counts = counts,
  mesh = mesh,
  prior = prior,
  initial = initial,
  proposal_covariance = proposal_covariance,
  detection_probability = detection_probability,
  particles = particles,
  iterations = iterations,
  warmup = warmup,
  seed = 20283101
)

summary <- summarize_draws_against_truth(fit$draws, inverse_ipm_truth())
summary$prior_mean <- prior$mean[summary$quantity]
summary$prior_sd <- prior$sd[summary$quantity]
summary$contraction <- 1 - summary$posterior_sd / summary$prior_sd
derived_draws <- t(apply(fit$draws, 1L, function(parameter) {
  derived_inverse_ipm_quantities(parameter, mesh)
}))
derived_summary <- summarize_draws_against_truth(
  derived_draws,
  derived_inverse_ipm_quantities(inverse_ipm_truth(), mesh)
)

saveRDS(
  list(
    fit = fit,
    prior = prior,
    counts = counts,
    mesh = mesh,
    detection_probability = detection_probability,
    populations = populations
  ),
  file.path(results_directory, "fit.rds"),
  compress = "xz"
)
write.csv(
  summary,
  file.path(results_directory, "parameter_summary.csv"),
  row.names = FALSE
)
write.csv(
  derived_summary,
  file.path(results_directory, "derived_summary.csv"),
  row.names = FALSE
)

cat(
  "Explicit PMMH acceptance:", round(fit$acceptance, 3),
  "\nVital-rate coverage:", sum(summary$covers_truth), "of", nrow(summary),
  "\nMedian contraction:", round(median(summary$contraction), 3),
  "\n"
)


source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")
source("R/profile_process_alternatives.R")
source("R/state_space_ipm.R")
source("R/explicit_state_space_ipm.R")

particles <- as.integer(Sys.getenv("EXPLICIT_MATCHED_PARTICLES", "500"))
iterations <- as.integer(Sys.getenv("EXPLICIT_MATCHED_ITERATIONS", "2000"))
warmup <- as.integer(Sys.getenv(
  "EXPLICIT_MATCHED_WARMUP",
  as.character(floor(iterations / 2))
))
populations <- 6L
transitions <- 4L
detection_probability <- 0.25
mesh <- seq(25, 155, by = 10)
results_directory <- "results/explicit_state_space_matched_pilot"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

set.seed(20284001)
delta <- mean(diff(mesh))
breaks <- c(mesh - delta / 2, tail(mesh, 1L) + delta / 2)
centers <- seq(60, 117, length.out = populations)
initial_states <- t(vapply(centers, function(center) {
  counts <- hist(
    pmin(155, pmax(25, rnorm(500, center, 6))),
    breaks = breaks,
    plot = FALSE
  )$counts
  # Make the conditioned 25% baseline profile recover the true initial state.
  as.integer(4L * round(counts / 4))
}, integer(length(mesh))))

simulation <- simulate_explicit_state_space_profiles(
  initial_states,
  inverse_ipm_truth(),
  mesh,
  transitions,
  detection_probability,
  seed = 20284002,
  conditioned_baseline = TRUE
)
counts <- simulation$observed
saveRDS(
  list(
    simulation = simulation,
    initial_states = initial_states,
    mesh = mesh,
    detection_probability = detection_probability
  ),
  file.path(results_directory, "simulation.rds"),
  compress = "xz"
)

initializer_prior <- profile_process_prior("one_step_gamma_poisson", "weak")
initializer_log_posterior <- function(parameter) {
  profile_process_log_posterior(
    parameter,
    counts,
    mesh,
    initializer_prior,
    "one_step_gamma_poisson"
  )
}
initializer <- fit_laplace_custom(
  initializer_log_posterior,
  initializer_prior,
  "one_step_gamma_poisson",
  maxit = 1200L
)
parameter_names <- names(inverse_ipm_truth())
prior <- default_inverse_ipm_priors("weak", "poisson")
initial <- initializer$map[parameter_names]
proposal_covariance <- initializer$covariance[
  parameter_names,
  parameter_names,
  drop = FALSE
]

likelihood_replicates <- vapply(seq_len(20L), function(index) {
  explicit_state_space_particle_filter(
    initial,
    counts,
    mesh,
    detection_probability,
    particles,
    seed = 20284100 + index
  )$log_likelihood
}, numeric(1))
write.csv(
  data.frame(log_likelihood = likelihood_replicates),
  file.path(results_directory, "likelihood_variance.csv"),
  row.names = FALSE
)
cat(
  "Matched explicit particle log-likelihood SD:",
  round(sd(likelihood_replicates), 3),
  "\n"
)

fit <- run_explicit_state_space_pmmh(
  counts,
  mesh,
  prior,
  initial,
  proposal_covariance,
  detection_probability,
  particles,
  iterations,
  warmup,
  seed = 20284201
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
    initializer = initializer,
    counts = counts,
    mesh = mesh,
    detection_probability = detection_probability
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
  "Matched explicit PMMH acceptance:", round(fit$acceptance, 3),
  "\nVital-rate coverage:", sum(summary$covers_truth), "of", nrow(summary),
  "\nMedian contraction:", round(median(summary$contraction), 3),
  "\n"
)

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")
source("R/profile_process_alternatives.R")
source("R/state_space_ipm.R")

mesh_step <- as.numeric(Sys.getenv("STATE_SPACE_MESH_STEP", "10"))
particles <- as.integer(Sys.getenv("STATE_SPACE_PARTICLES", "300"))
iterations <- as.integer(Sys.getenv("STATE_SPACE_ITERATIONS", "2000"))
warmup <- as.integer(Sys.getenv(
  "STATE_SPACE_WARMUP",
  as.character(floor(iterations / 2))
))
mesh <- seq(25, 155, by = mesh_step)
detection_probability <- 0.25
populations <- 8L
transitions <- 4L
initial_n <- 500L
results_directory <- "results/state_space_ipm_pilot"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

make_parent_excitation_population <- function(population_index) {
  initial <- simulate_initial_population(n = initial_n, latent_quality_sd = 0)
  center <- seq(60, 140, length.out = populations)[population_index]
  initial$size <- pmin(145, pmax(35, rnorm(initial_n, center, 4)))
  initial$age <- pmax(0L, round((initial$size - 35) / 18))
  initial
}

set.seed(20280608)
simulations <- lapply(seq_len(populations), function(population_index) {
  simulate_population(
    years = transitions + 1L,
    initial_population = make_parent_excitation_population(population_index),
    vital_rates = default_fish_vital_rates(),
    environment = rep(0, transitions)
  )
})
counts <- lapply(simulations, function(simulation) {
  profiles <- sample_size_profiles(
    simulation$census,
    detection = function(size, year) {
      rep(detection_probability, length(size))
    },
    measurement_sd = 0
  )
  make_profile_count_matrix(profiles, mesh)
})
saveRDS(
  list(
    simulations = simulations,
    counts = counts,
    mesh = mesh,
    detection_probability = detection_probability
  ),
  file.path(results_directory, "simulation.rds"),
  compress = "xz"
)

cat("Fitting one-step initializer.\n")
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

demographic_names <- names(inverse_ipm_truth())
prior <- state_space_ipm_prior("weak")
initial <- c(
  initializer$map[demographic_names],
  log_state_fano = 0
)
proposal_covariance <- matrix(
  0,
  nrow = length(initial),
  ncol = length(initial),
  dimnames = list(names(initial), names(initial))
)
proposal_covariance[demographic_names, demographic_names] <-
  initializer$covariance[
    demographic_names,
    demographic_names,
    drop = FALSE
  ]
proposal_covariance["log_state_fano", "log_state_fano"] <- 0.5^2

cat("Checking particle-likelihood variance.\n")
likelihood_replicates <- vapply(seq_len(20L), function(index) {
  state_space_particle_filter(
    initial,
    counts,
    mesh,
    detection_probability,
    particles,
    seed = 20280900 + index
  )$log_likelihood
}, numeric(1))
write.csv(
  data.frame(log_likelihood = likelihood_replicates),
  file.path(results_directory, "likelihood_variance.csv"),
  row.names = FALSE
)
cat(
  "Particle log-likelihood SD:",
  round(sd(likelihood_replicates), 3),
  "\n"
)

cat(
  "Running PMMH with", particles, "particles and", iterations,
  "iterations.\n"
)
fit <- run_state_space_pmmh(
  counts = counts,
  mesh = mesh,
  prior = prior,
  initial = initial,
  proposal_covariance = proposal_covariance,
  detection_probability = detection_probability,
  particles = particles,
  iterations = iterations,
  warmup = warmup,
  seed = 20281001
)
saveRDS(
  list(
    fit = fit,
    prior = prior,
    initializer = initializer,
    mesh = mesh,
    counts = counts,
    detection_probability = detection_probability
  ),
  file.path(results_directory, "fit.rds"),
  compress = "xz"
)

truth <- c(inverse_ipm_truth(), log_state_fano = NA_real_)
summary <- summarize_draws_against_truth(fit$draws, truth)
summary$prior_mean <- prior$mean[summary$quantity]
summary$prior_sd <- prior$sd[summary$quantity]
summary$contraction <- 1 - summary$posterior_sd / summary$prior_sd
write.csv(
  summary,
  file.path(results_directory, "parameter_summary.csv"),
  row.names = FALSE
)

derived_draws <- t(apply(fit$draws[, demographic_names, drop = FALSE], 1L, function(
  parameter
) {
  derived_inverse_ipm_quantities(parameter, mesh)
}))
derived_truth <- derived_inverse_ipm_quantities(inverse_ipm_truth(), mesh)
derived_summary <- summarize_draws_against_truth(derived_draws, derived_truth)
write.csv(
  derived_summary,
  file.path(results_directory, "derived_summary.csv"),
  row.names = FALSE
)

representative <- state_space_particle_filter(
  setNames(colMeans(fit$draws), colnames(fit$draws)),
  counts,
  mesh,
  detection_probability,
  particles = max(1000L, particles),
  seed = 20281101,
  return_filtered = TRUE
)
saveRDS(
  representative$filtered,
  file.path(results_directory, "filtered_states.rds")
)

cat(
  "PMMH acceptance:", round(fit$acceptance, 3),
  "\nMedian vital-rate contraction:",
  round(median(summary$contraction[summary$quantity %in% demographic_names]), 3),
  "\n"
)


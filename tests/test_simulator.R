source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/identifiability.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/transportability.R")
source("R/correlated_process_ipm.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")
source("R/profile_process_alternatives.R")
source("R/maturity_structured_ipm.R")
source("R/state_space_ipm.R")
source("R/explicit_state_space_ipm.R")

set.seed(42)

initial <- simulate_initial_population(n = 200, latent_quality_sd = 0.2)
simulation <- simulate_population(
  years = 5,
  initial_population = initial,
  environment = rep(0, 4)
)

stopifnot(
  nrow(simulation$census) > 0,
  nrow(simulation$transitions) > 0,
  all(simulation$census$size >= 25),
  all(simulation$census$size <= 160),
  all(simulation$transitions$survival_probability >= 0),
  all(simulation$transitions$survival_probability <= 1)
)

profiles <- sample_size_profiles(simulation$census, measurement_sd = 0.5)
captures <- sample_mark_recapture(simulation$census, measurement_sd = 0.5)
stopifnot(nrow(profiles) > 0, nrow(captures) > 0)

mesh <- seq(25, 160, by = 1)
ipm <- make_ipm_kernel(
  mesh = mesh,
  survival = function(size) inv_logit(-2.2 + 0.032 * size),
  growth_mean = function(size) pmax(size, size + 18 - 0.12 * size),
  growth_sd = function(size) 4.5,
  recruitment = function(size) 0.5 * exp(-7 + 0.075 * size),
  recruit_mean = function(size) 43 + 0.02 * size,
  recruit_sd = function(size) 3.5
)

stopifnot(
  all(ipm$kernel >= 0),
  max(abs(colSums(ipm$projection) * ipm$delta - ipm$survival)) < 1e-10,
  max(abs(colSums(ipm$fertility) * ipm$delta - ipm$recruitment)) < 1e-10
)

projection <- project_ipm(rep(1, length(mesh)), ipm, transitions = 3)
stopifnot(all(projection >= 0), nrow(projection) == 4)

jacobian <- finite_difference_jacobian(
  c(a = 1, b = 2),
  function(parameter) c(parameter["a"] + parameter["b"], parameter["a"]^2)
)
diagnostic <- sensitivity_svd(jacobian)
stopifnot(
  all(dim(jacobian) == c(2, 2)),
  diagnostic$numerical_rank == 2,
  is.finite(diagnostic$condition_number)
)

truth <- inverse_ipm_truth()
prior <- default_inverse_ipm_priors("weak")
process_prior <- default_inverse_ipm_priors("weak", "gamma_poisson")
centered_ipm <- build_centered_ipm(truth, mesh)
expected <- expected_profile_counts(truth, rep(1, length(mesh)), mesh, 2)
stopifnot(
  all(centered_ipm$kernel >= 0),
  all(names(truth) == names(prior$mean)),
  "log_process_precision" %in% names(process_prior$mean),
  all(dim(expected) == c(3, length(mesh)))
)

external_prior <- make_external_information_prior(
  "survival_growth",
  "gamma_poisson"
)
biased_prior <- make_external_information_prior(
  "all_vital_rates_biased",
  "poisson"
)
stopifnot(
  external_prior$sd["survival_at_80"] < prior$sd["survival_at_80"],
  external_prior$sd["log_recruitment_at_80"] ==
    process_prior$sd["log_recruitment_at_80"],
  biased_prior$mean["survival_at_80"] != truth["survival_at_80"]
)

fecundity_prior <- make_fecundity_study_prior(
  sample_size = 50,
  accuracy = "biased",
  borrowing = "robust",
  seed = 8
)
stopifnot(
  length(fecundity_prior$blocks) == 1L,
  length(fecundity_prior$alternative_starts) == 2L,
  is.finite(inverse_ipm_log_prior(truth, fecundity_prior)),
  external_block_responsibility(
    truth,
    fecundity_prior$blocks$fecundity_study
  ) >= 0,
  external_block_responsibility(
    truth,
    fecundity_prior$blocks$fecundity_study
  ) <= 1
)

process_basis <- make_smooth_process_basis(mesh, rank = 4)
stopifnot(
  all(dim(process_basis) == c(length(mesh), 4)),
  all(is.finite(process_basis))
)

laplace_counts <- make_profile_count_matrix(simulation$census, mesh, 2)
laplace_fit <- fit_laplace_inverse_ipm(
  laplace_counts,
  mesh,
  default_inverse_ipm_priors("weak"),
  "poisson",
  maxit = 100
)
laplace_draws <- sample_laplace_draws(laplace_fit, draws = 20, seed = 4)
stopifnot(
  nrow(laplace_draws) == 20,
  all(colnames(laplace_draws) == names(inverse_ipm_truth()))
)

robust_fit <- fit_robust_laplace_inverse_ipm(
  laplace_counts,
  mesh,
  fecundity_prior,
  "poisson",
  maxit = 100
)
robust_draws <- sample_robust_laplace_draws(
  robust_fit,
  draws = 20,
  seed = 5
)
stopifnot(
  robust_fit$posterior_external_weight >= 0,
  robust_fit$posterior_external_weight <= 1,
  nrow(robust_draws) == 20
)

truths <- transportability_truths()
alternative_rates <- fish_vital_rates_from_inverse_ipm(truths$low_survival)
alternative_survival <- alternative_rates$survival(
  80, 1, "F", 0, 0, 500
)
stopifnot(
  length(truths) == 3L,
  abs(alternative_survival - 0.40) < 1e-10,
  all(names(fecundity_transport_bias(1)) %in% names(truth))
)

benchmark_events <- prepare_capture_recapture_events(
  captures,
  final_year = max(simulation$census$year),
  recapture_detection = 0.35
)
complete_history <- prepare_complete_history(simulation)
stopifnot(
  nrow(benchmark_events) > 0,
  all(benchmark_events$recaptured == !is.na(benchmark_events$to_size)),
  nrow(complete_history$transitions) == nrow(simulation$transitions),
  sum(complete_history$transitions$offspring_count) == nrow(simulation$recruits),
  is.finite(capture_recapture_log_likelihood(truth, benchmark_events)),
  is.finite(complete_history_log_likelihood(truth, complete_history))
)

alternative_counts <- make_profile_count_matrix(simulation$census, mesh, 2)
shared_prior <- profile_process_prior("one_step_shared_gamma_poisson")
finite_prior <- profile_process_prior("one_step_finite_population")
finite_overdispersed_prior <- profile_process_prior(
  "one_step_finite_population_overdispersed"
)
finite_moments <- finite_population_transition_moments(
  alternative_counts[1L, ],
  centered_ipm,
  detection_probability = 1
)
stopifnot(
  all(one_step_expected_counts(
    truth,
    alternative_counts[1L, ],
    mesh
  ) >= 0),
  is.finite(profile_process_log_posterior(
    truth,
    alternative_counts,
    mesh,
    prior,
    "one_step_poisson"
  )),
  is.finite(profile_process_log_posterior(
    shared_prior$mean,
    alternative_counts,
    mesh,
    shared_prior,
    "one_step_shared_gamma_poisson"
  )),
  max(abs(
    finite_moments$mean -
      one_step_expected_counts(truth, alternative_counts[1L, ], mesh)
  )) < 1e-8,
  max(abs(finite_moments$covariance - t(finite_moments$covariance))) < 1e-8,
  is.finite(profile_process_log_posterior(
    finite_prior$mean,
    alternative_counts,
    mesh,
    finite_prior,
    "one_step_finite_population"
  )),
  is.finite(profile_process_log_posterior(
    finite_overdispersed_prior$mean,
    alternative_counts,
    mesh,
    finite_overdispersed_prior,
    "one_step_finite_population_overdispersed",
    detection_probability = 0.6
  ))
)

maturity_ipm <- build_maturity_centered_ipm(truth, mesh)
maturity_rates <- fish_vital_rates_from_maturity_ipm(truth)
stopifnot(
  all(maturity_ipm$kernel >= 0),
  maturity_rates$realized_recruitment(50, 1, "F", 0, 0, 500) <
    maturity_rates$realized_recruitment(80, 1, "F", 0, 0, 500),
  abs(
    maturity_inverse_ipm_quantities(truth, mesh)["recruitment_80"] -
      derived_inverse_ipm_quantities(truth, mesh)["recruitment_80"]
  ) < 1e-10
)

state_space_prior <- state_space_ipm_prior("weak")
state_space_parameter <- c(truth, log_state_fano = 0)
state_space_filter <- state_space_particle_filter(
  state_space_parameter,
  alternative_counts,
  mesh,
  detection_probability = 1,
  particles = 50,
  seed = 17,
  return_filtered = TRUE,
  proposal = "adapted"
)
stopifnot(
  all(state_space_transition_matrix(state_space_parameter, mesh) >= 0),
  is.finite(state_space_filter$log_likelihood),
  length(state_space_filter$filtered) == 1L,
  all(dim(state_space_filter$filtered[[1L]]) == dim(alternative_counts)),
  is.finite(state_space_ipm_log_posterior_estimate(
    state_space_parameter,
    alternative_counts,
    mesh,
    state_space_prior,
    particles = 50,
    seed = 18,
    proposal = "adapted"
  ))
)

explicit_components <- explicit_transition_components(truth, mesh)
explicit_states <- matrix(
  rep(alternative_counts[1L, ], each = 50L),
  nrow = 50L
)
set.seed(19)
explicit_next <- simulate_explicit_ipm_transition(
  explicit_states,
  explicit_components
)
explicit_filter <- explicit_state_space_particle_filter(
  truth,
  alternative_counts,
  mesh,
  detection_probability = 1,
  particles = 50,
  seed = 20,
  return_filtered = TRUE
)
explicit_simulation <- simulate_explicit_state_space_profiles(
  explicit_states[1L, , drop = FALSE],
  truth,
  mesh,
  transitions = 2L,
  detection_probability = 1,
  seed = 22
)
stopifnot(
  all(explicit_next >= 0),
  all(explicit_next == round(explicit_next)),
  is.finite(explicit_filter$log_likelihood),
  length(explicit_filter$filtered) == 1L,
  length(explicit_filter$effective_sample_sizes) ==
    nrow(alternative_counts) - 1L,
  length(explicit_simulation$latent) == 1L,
  all(dim(explicit_simulation$latent[[1L]]) == c(3L, length(mesh))),
  all(explicit_simulation$latent[[1L]] >= 0),
  is.finite(explicit_state_space_log_posterior_estimate(
    truth,
    alternative_counts,
    mesh,
    prior,
    particles = 50,
    seed = 21
  ))
)

cat("All simulator scaffold checks passed.\n")

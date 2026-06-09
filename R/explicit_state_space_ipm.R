# Explicit finite-population state-space IPM.

explicit_transition_components <- function(
  parameter,
  mesh,
  kernel_builder = build_centered_ipm
) {
  ipm <- kernel_builder(parameter, mesh)
  projection <- ipm$delta * ipm$projection
  fertility <- ipm$delta * ipm$fertility
  if (
    any(!is.finite(projection)) ||
      any(!is.finite(fertility)) ||
      any(projection < 0) ||
      any(fertility < 0) ||
      any(colSums(projection) > 1 + 1e-8)
  ) {
    stop("Invalid finite-population transition components.")
  }
  list(
    projection = projection,
    fertility = fertility,
    transition = projection + fertility
  )
}

simulate_explicit_ipm_transition <- function(states, components) {
  stopifnot(
    is.matrix(states),
    all(is.finite(states)),
    all(states >= 0),
    all(states == round(states))
  )
  particles <- nrow(states)
  bins <- ncol(states)
  stopifnot(all(dim(components$projection) == c(bins, bins)))

  survivors <- matrix(0, nrow = particles, ncol = bins)
  recruits <- matrix(0, nrow = particles, ncol = bins)

  for (source in seq_len(bins)) {
    remaining_count <- states[, source]
    remaining_probability <- 1
    for (destination in seq_len(bins)) {
      probability <- components$projection[destination, source]
      conditional_probability <- if (remaining_probability > 0) {
        probability / remaining_probability
      } else {
        0
      }
      conditional_probability <- pmin(1, pmax(0, conditional_probability))
      moved <- rbinom(
        particles,
        size = remaining_count,
        prob = conditional_probability
      )
      survivors[, destination] <- survivors[, destination] + moved
      remaining_count <- remaining_count - moved
      remaining_probability <- remaining_probability - probability
    }

    # Marginalize unobserved sex composition within each source size bin.
    # Females produce at twice the population-average fertility rate.
    females <- rbinom(particles, states[, source], 0.5)
    recruit_rate <- 2 * outer(
      females,
      components$fertility[, source]
    )
    recruits <- recruits + matrix(
      rpois(length(recruit_rate), lambda = as.vector(recruit_rate)),
      nrow = particles
    )
  }

  survivors + recruits
}

simulate_guided_explicit_ipm_transition <- function(
  states,
  components,
  observed,
  detection_probability,
  guide_power = 0.75,
  guide_bounds = c(0.1, 10)
) {
  stopifnot(
    is.matrix(states),
    length(observed) == ncol(states),
    length(guide_bounds) == 2L,
    guide_bounds[1L] > 0,
    guide_bounds[2L] >= guide_bounds[1L]
  )
  particles <- nrow(states)
  bins <- ncol(states)
  predicted_mean <- states %*% t(components$transition)
  guide <- (
    (rep(observed, each = particles) + 0.5) /
      (detection_probability * predicted_mean + 0.5)
  )^guide_power
  guide <- matrix(
    pmin(guide_bounds[2L], pmax(guide_bounds[1L], guide)),
    nrow = particles,
    ncol = bins
  )

  survivors <- matrix(0, nrow = particles, ncol = bins)
  recruits <- matrix(0, nrow = particles, ncol = bins)
  log_target_over_proposal <- numeric(particles)

  for (source in seq_len(bins)) {
    target_probability <- components$projection[, source]
    target_death <- pmax(0, 1 - sum(target_probability))
    proposal_mass <- sweep(guide, 2L, target_probability, "*")
    normalizer <- rowSums(proposal_mass) + target_death
    proposal_probability <- proposal_mass / normalizer
    proposal_death <- target_death / normalizer

    remaining_count <- states[, source]
    remaining_probability <- rep(1, particles)
    for (destination in seq_len(bins)) {
      probability <- proposal_probability[, destination]
      conditional_probability <- ifelse(
        remaining_probability > 0,
        probability / remaining_probability,
        0
      )
      conditional_probability <- pmin(1, pmax(0, conditional_probability))
      moved <- rbinom(
        particles,
        size = remaining_count,
        prob = conditional_probability
      )
      survivors[, destination] <- survivors[, destination] + moved
      if (target_probability[destination] > 0) {
        log_target_over_proposal <- log_target_over_proposal +
          moved * (
            log(target_probability[destination]) -
              log(pmax(probability, 1e-300))
          )
      }
      remaining_count <- remaining_count - moved
      remaining_probability <- remaining_probability - probability
    }
    if (target_death > 0) {
      log_target_over_proposal <- log_target_over_proposal +
        remaining_count * (
          log(target_death) -
            log(pmax(proposal_death, 1e-300))
        )
    }

    females <- rbinom(particles, states[, source], 0.5)
    target_rate <- 2 * outer(
      females,
      components$fertility[, source]
    )
    proposal_rate <- target_rate * guide
    proposed_recruits <- matrix(
      rpois(length(proposal_rate), lambda = as.vector(proposal_rate)),
      nrow = particles
    )
    recruits <- recruits + proposed_recruits
    positive <- target_rate > 0
    log_target_over_proposal <- log_target_over_proposal +
      rowSums(
        ifelse(
          positive,
          proposed_recruits * (-log(guide)) +
            target_rate * (guide - 1),
          0
        )
      )
  }

  list(
    states = survivors + recruits,
    log_target_over_proposal = log_target_over_proposal
  )
}

explicit_transition_marginal_moments <- function(states, components) {
  expected <- states %*% t(components$transition)
  projection_variance <- components$projection *
    (1 - components$projection)
  # If female status is unobserved, recruitment variance from source j to
  # destination i is N_j(F_ij + F_ij^2).
  fertility_variance <- components$fertility + components$fertility^2
  variance <- states %*% t(projection_variance + fertility_variance)
  list(mean = expected, variance = variance)
}

simulate_explicit_state_space_profiles <- function(
  initial_states,
  parameter,
  mesh,
  transitions,
  detection_probability = 1,
  seed = 1L,
  conditioned_baseline = TRUE,
  kernel_builder = build_centered_ipm
) {
  if (is.vector(initial_states)) {
    initial_states <- matrix(initial_states, nrow = 1L)
  }
  stopifnot(
    is.matrix(initial_states),
    ncol(initial_states) == length(mesh),
    all(initial_states >= 0),
    all(initial_states == round(initial_states)),
    transitions >= 1L
  )
  set.seed(seed)
  components <- explicit_transition_components(
    parameter,
    mesh,
    kernel_builder
  )

  latent <- vector("list", nrow(initial_states))
  observed <- vector("list", nrow(initial_states))
  for (population in seq_len(nrow(initial_states))) {
    states <- matrix(
      NA_real_,
      nrow = transitions + 1L,
      ncol = length(mesh)
    )
    profiles <- states
    states[1L, ] <- initial_states[population, ]
    profiles[1L, ] <- if (conditioned_baseline) {
      round(detection_probability * states[1L, ])
    } else {
      rpois(length(mesh), detection_probability * states[1L, ])
    }

    for (transition in seq_len(transitions)) {
      states[transition + 1L, ] <- simulate_explicit_ipm_transition(
        states[transition, , drop = FALSE],
        components
      )
      profiles[transition + 1L, ] <- rpois(
        length(mesh),
        detection_probability * states[transition + 1L, ]
      )
    }
    latent[[population]] <- states
    observed[[population]] <- profiles
  }
  list(latent = latent, observed = observed)
}

explicit_observation_log_likelihood <- function(
  states,
  observed,
  detection_probability
) {
  particles <- nrow(states)
  rowSums(matrix(
    dpois(
      rep(observed, each = particles),
      lambda = pmax(detection_probability * states, 1e-12),
      log = TRUE
    ),
    nrow = particles
  ))
}

explicit_predictive_log_likelihood <- function(
  states,
  observed,
  components,
  detection_probability
) {
  moments <- explicit_transition_marginal_moments(states, components)
  predictive_mean <- detection_probability * moments$mean
  predictive_variance <- predictive_mean +
    detection_probability^2 * moments$variance
  particles <- nrow(states)
  log_likelihood <- matrix(0, nrow = particles, ncol = ncol(states))

  for (bin in seq_len(ncol(states))) {
    overdispersed <- predictive_variance[, bin] >
      predictive_mean[, bin] * (1 + 1e-8)
    if (any(overdispersed)) {
      mean <- pmax(predictive_mean[overdispersed, bin], 1e-12)
      variance <- predictive_variance[overdispersed, bin]
      size <- mean^2 / pmax(variance - mean, 1e-12)
      log_likelihood[overdispersed, bin] <- dnbinom(
        observed[bin],
        size = size,
        mu = mean,
        log = TRUE
      )
    }
    if (any(!overdispersed)) {
      log_likelihood[!overdispersed, bin] <- dpois(
        observed[bin],
        lambda = pmax(predictive_mean[!overdispersed, bin], 1e-12),
        log = TRUE
      )
    }
  }
  rowSums(log_likelihood)
}

normalized_log_weights <- function(log_weight) {
  maximum <- max(log_weight)
  if (!is.finite(maximum)) {
    return(NULL)
  }
  weight <- exp(log_weight - maximum)
  weight / sum(weight)
}

explicit_state_space_particle_filter <- function(
  parameter,
  counts,
  mesh,
  detection_probability = 1,
  particles = 1000L,
  seed = NULL,
  return_filtered = FALSE,
  guide_power = 0,
  guide_bounds = c(0.1, 10),
  kernel_builder = build_centered_ipm
) {
  counts <- as_profile_trajectory_list(counts)
  if (
    particles < 2L ||
      any(!is.finite(parameter)) ||
      length(detection_probability) != 1L ||
      !is.finite(detection_probability) ||
      detection_probability <= 0 ||
      detection_probability > 1
  ) {
    return(list(log_likelihood = -Inf, filtered = NULL))
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  components <- tryCatch(
    explicit_transition_components(parameter, mesh, kernel_builder),
    error = function(error) NULL
  )
  if (is.null(components)) {
    return(list(log_likelihood = -Inf, filtered = NULL))
  }

  log_likelihood <- 0
  filtered <- if (return_filtered) vector("list", length(counts)) else NULL
  effective_sample_sizes <- numeric()

  for (trajectory_index in seq_along(counts)) {
    trajectory <- counts[[trajectory_index]]
    baseline <- round(trajectory[1L, ] / detection_probability)
    state <- matrix(
      rep(baseline, each = particles),
      nrow = particles
    )
    if (return_filtered) {
      filtered_trajectory <- matrix(
        NA_real_,
        nrow = nrow(trajectory),
        ncol = ncol(trajectory)
      )
      filtered_trajectory[1L, ] <- baseline
    }

    for (transition in seq_len(nrow(trajectory) - 1L)) {
      observed <- trajectory[transition + 1L, ]

      # Auxiliary look-ahead resampling directs particles toward source states
      # likely to produce the next observed profile.
      lookahead <- explicit_predictive_log_likelihood(
        state,
        observed,
        components,
        detection_probability
      )
      lookahead_weights <- normalized_log_weights(lookahead)
      if (is.null(lookahead_weights)) {
        return(list(log_likelihood = -Inf, filtered = NULL))
      }
      first_increment <- particle_log_mean_weight(lookahead)
      ancestors <- systematic_resample(lookahead_weights)
      proposal <- simulate_guided_explicit_ipm_transition(
        state[ancestors, , drop = FALSE],
        components,
        observed,
        detection_probability,
        guide_power,
        guide_bounds
      )
      state <- proposal$states

      observation_log_likelihood <- explicit_observation_log_likelihood(
        state,
        observed,
        detection_probability
      )
      correction <- observation_log_likelihood +
        proposal$log_target_over_proposal -
        lookahead[ancestors]
      correction_weights <- normalized_log_weights(correction)
      if (is.null(correction_weights)) {
        return(list(log_likelihood = -Inf, filtered = NULL))
      }
      log_likelihood <- log_likelihood +
        first_increment +
        particle_log_mean_weight(correction)
      effective_sample_sizes <- c(
        effective_sample_sizes,
        1 / sum(correction_weights^2)
      )

      if (return_filtered) {
        filtered_trajectory[transition + 1L, ] <-
          colSums(state * correction_weights)
      }
      state <- state[
        systematic_resample(correction_weights),
        ,
        drop = FALSE
      ]
    }
    if (return_filtered) {
      filtered[[trajectory_index]] <- filtered_trajectory
    }
  }

  list(
    log_likelihood = log_likelihood,
    filtered = filtered,
    effective_sample_sizes = effective_sample_sizes
  )
}

explicit_state_space_log_posterior_estimate <- function(
  parameter,
  counts,
  mesh,
  prior,
  detection_probability = 1,
  particles = 1000L,
  seed = NULL,
  kernel_builder = build_centered_ipm
) {
  log_prior <- inverse_ipm_log_prior(parameter, prior)
  if (!is.finite(log_prior)) {
    return(-Inf)
  }
  filter <- explicit_state_space_particle_filter(
    parameter,
    counts,
    mesh,
    detection_probability,
    particles,
    seed,
    FALSE,
    kernel_builder = kernel_builder
  )
  filter$log_likelihood + log_prior
}

run_explicit_state_space_pmmh <- function(
  counts,
  mesh,
  prior,
  initial,
  proposal_covariance,
  detection_probability = 1,
  particles = 1000L,
  iterations = 4000L,
  warmup = 2000L,
  seed = 1L,
  target_acceptance = 0.12,
  kernel_builder = build_centered_ipm
) {
  stopifnot(iterations > warmup)
  set.seed(seed)
  dimension <- length(initial)
  parameter_names <- names(initial)
  base_chol <- chol(
    proposal_covariance[parameter_names, parameter_names, drop = FALSE] +
      diag(1e-10, dimension)
  )
  log_scale <- log(0.35 * 2.38 / sqrt(dimension))

  current <- initial
  current_seed <- sample.int(.Machine$integer.max, 1L)
  current_log_posterior <- explicit_state_space_log_posterior_estimate(
    current,
    counts,
    mesh,
    prior,
    detection_probability,
    particles,
    current_seed,
    kernel_builder
  )
  if (!is.finite(current_log_posterior)) {
    stop("Initial explicit-state PMMH value is not finite.")
  }

  draws <- matrix(
    NA_real_,
    nrow = iterations,
    ncol = dimension,
    dimnames = list(NULL, parameter_names)
  )
  log_posterior <- numeric(iterations)
  accepted <- logical(iterations)

  for (iteration in seq_len(iterations)) {
    candidate <- current + as.vector(
      rnorm(dimension) %*% base_chol
    ) * exp(log_scale)
    names(candidate) <- parameter_names
    candidate_seed <- sample.int(.Machine$integer.max, 1L)
    candidate_log_posterior <- explicit_state_space_log_posterior_estimate(
      candidate,
      counts,
      mesh,
      prior,
      detection_probability,
      particles,
      candidate_seed,
      kernel_builder
    )
    if (
      is.finite(candidate_log_posterior) &&
        log(runif(1)) < candidate_log_posterior - current_log_posterior
    ) {
      current <- candidate
      current_seed <- candidate_seed
      current_log_posterior <- candidate_log_posterior
      accepted[iteration] <- TRUE
    }
    draws[iteration, ] <- current
    log_posterior[iteration] <- current_log_posterior

    if (iteration <= warmup) {
      learning_rate <- min(0.03, 1 / sqrt(iteration))
      log_scale <- log_scale +
        learning_rate * (accepted[iteration] - target_acceptance)
    }
  }

  retained <- (warmup + 1L):iterations
  list(
    draws = draws[retained, , drop = FALSE],
    log_posterior = log_posterior[retained],
    acceptance = mean(accepted[retained]),
    warmup_acceptance = mean(accepted[seq_len(warmup)]),
    final_scale = exp(log_scale),
    particles = particles,
    detection_probability = detection_probability
  )
}

# Particle-filtered state-space IPM with Poisson profile observations.

state_space_ipm_prior <- function(
  regime = c("weak", "informative", "informative_wrong")
) {
  regime <- match.arg(regime)
  demographic <- default_inverse_ipm_priors(regime, "poisson")
  list(
    regime = regime,
    process_model = "gamma_state_poisson_observation",
    mean = c(demographic$mean, log_state_fano = 0),
    sd = c(demographic$sd, log_state_fano = 1)
  )
}

state_space_transition_matrix <- function(
  parameter,
  mesh,
  kernel_builder = build_centered_ipm
) {
  ipm <- kernel_builder(parameter, mesh)
  ipm$delta * ipm$kernel
}

systematic_resample <- function(weights) {
  count <- length(weights)
  positions <- (runif(1) + seq.int(0, count - 1L)) / count
  findInterval(positions, cumsum(weights)) + 1L
}

particle_log_mean_weight <- function(log_weight) {
  maximum <- max(log_weight)
  if (!is.finite(maximum)) {
    return(-Inf)
  }
  maximum + log(mean(exp(log_weight - maximum)))
}

state_space_particle_filter <- function(
  parameter,
  counts,
  mesh,
  detection_probability = 1,
  particles = 500L,
  seed = NULL,
  return_filtered = FALSE,
  proposal = c("adapted", "bootstrap"),
  kernel_builder = build_centered_ipm
) {
  proposal <- match.arg(proposal)
  counts <- as_profile_trajectory_list(counts)
  if (
    particles < 2L ||
      any(!is.finite(parameter)) ||
      !"log_state_fano" %in% names(parameter) ||
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

  transition_matrix <- tryCatch(
    state_space_transition_matrix(parameter, mesh, kernel_builder),
    error = function(error) NULL
  )
  state_fano <- exp(parameter["log_state_fano"])
  if (
    is.null(transition_matrix) ||
      any(!is.finite(transition_matrix)) ||
      any(transition_matrix < 0) ||
      !is.finite(state_fano) ||
      state_fano <= 0
  ) {
    return(list(log_likelihood = -Inf, filtered = NULL))
  }

  log_likelihood <- 0
  filtered <- if (return_filtered) vector("list", length(counts)) else NULL

  for (trajectory_index in seq_along(counts)) {
    trajectory <- counts[[trajectory_index]]
    state <- matrix(
      rep(
        pmax(trajectory[1L, ] / detection_probability, 1e-8),
        each = particles
      ),
      nrow = particles
    )
    if (return_filtered) {
      filtered_trajectory <- matrix(
        NA_real_,
        nrow = nrow(trajectory),
        ncol = ncol(trajectory)
      )
      filtered_trajectory[1L, ] <- colMeans(state)
    }

    for (transition in seq_len(nrow(trajectory) - 1L)) {
      expected_state <- state %*% t(transition_matrix)
      expected_state <- pmax(expected_state, 1e-10)
      observed <- trajectory[transition + 1L, ]
      transition_shape <- expected_state / state_fano
      if (proposal == "adapted") {
        # The Gamma transition and Poisson observation are conjugate. This
        # fully adapted proposal uses the exact predictive likelihood and then
        # samples the latent state conditional on the new profile.
        log_weight <- rowSums(matrix(
          dnbinom(
            rep(observed, each = particles),
            size = as.vector(transition_shape),
            mu = as.vector(detection_probability * expected_state),
            log = TRUE
          ),
          nrow = particles
        ))
        state <- matrix(
          rgamma(
            length(expected_state),
            shape = as.vector(transition_shape) +
              rep(observed, each = particles),
            rate = 1 / state_fano + detection_probability
          ),
          nrow = particles
        )
      } else {
        state <- matrix(
          rgamma(
            length(expected_state),
            shape = as.vector(transition_shape),
            scale = state_fano
          ),
          nrow = particles
        )
        log_weight <- rowSums(matrix(
          dpois(
            rep(observed, each = particles),
            lambda = pmax(detection_probability * state, 1e-12),
            log = TRUE
          ),
          nrow = particles
        ))
      }
      increment <- particle_log_mean_weight(log_weight)
      if (!is.finite(increment)) {
        return(list(log_likelihood = -Inf, filtered = NULL))
      }
      log_likelihood <- log_likelihood + increment

      maximum <- max(log_weight)
      weight <- exp(log_weight - maximum)
      weight <- weight / sum(weight)
      if (return_filtered) {
        filtered_trajectory[transition + 1L, ] <-
          colSums(state * weight)
      }
      state <- state[systematic_resample(weight), , drop = FALSE]
    }

    if (return_filtered) {
      filtered[[trajectory_index]] <- filtered_trajectory
    }
  }

  list(log_likelihood = log_likelihood, filtered = filtered)
}

state_space_ipm_log_prior <- function(parameter, prior) {
  if (
    any(!names(parameter) %in% names(prior$mean)) ||
      any(!is.finite(parameter))
  ) {
    return(-Inf)
  }
  sum(dnorm(
    parameter,
    mean = prior$mean[names(parameter)],
    sd = prior$sd[names(parameter)],
    log = TRUE
  ))
}

state_space_ipm_log_posterior_estimate <- function(
  parameter,
  counts,
  mesh,
  prior,
  detection_probability = 1,
  particles = 500L,
  seed = NULL,
  proposal = "adapted",
  kernel_builder = build_centered_ipm
) {
  log_prior <- state_space_ipm_log_prior(parameter, prior)
  if (!is.finite(log_prior)) {
    return(-Inf)
  }
  filter <- state_space_particle_filter(
    parameter,
    counts,
    mesh,
    detection_probability,
    particles,
    seed,
    FALSE,
    proposal,
    kernel_builder
  )
  filter$log_likelihood + log_prior
}

run_state_space_pmmh <- function(
  counts,
  mesh,
  prior,
  initial,
  proposal_covariance,
  detection_probability = 1,
  particles = 500L,
  iterations = 4000L,
  warmup = 2000L,
  seed = 1L,
  target_acceptance = 0.15,
  proposal = "adapted",
  kernel_builder = build_centered_ipm
) {
  stopifnot(
    iterations > warmup,
    all(names(initial) %in% names(prior$mean)),
    all(dim(proposal_covariance) == length(initial))
  )
  set.seed(seed)
  dimension <- length(initial)
  parameter_names <- names(initial)
  base_chol <- chol(
    proposal_covariance[parameter_names, parameter_names, drop = FALSE] +
      diag(1e-10, dimension)
  )
  log_scale <- log(0.5 * 2.38 / sqrt(dimension))

  current <- initial
  current_seed <- sample.int(.Machine$integer.max, 1L)
  current_log_posterior <- state_space_ipm_log_posterior_estimate(
    current,
    counts,
    mesh,
    prior,
    detection_probability,
    particles,
    current_seed,
    proposal,
    kernel_builder
  )
  if (!is.finite(current_log_posterior)) {
    stop("Initial state-space PMMH value is not finite.")
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
    candidate_log_posterior <- state_space_ipm_log_posterior_estimate(
      candidate,
      counts,
      mesh,
      prior,
      detection_probability,
      particles,
      candidate_seed,
      proposal,
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

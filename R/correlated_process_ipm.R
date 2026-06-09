# Positive low-rank log-Gaussian process discrepancy for profile trajectories.

make_smooth_process_basis <- function(mesh, rank = 4L) {
  stopifnot(rank >= 1L, length(mesh) > rank)
  scaled <- (mesh - min(mesh)) / (max(mesh) - min(mesh))
  basis <- matrix(1, nrow = length(mesh), ncol = rank)

  if (rank >= 2L) {
    for (index in 2:rank) {
      frequency <- ceiling((index - 1L) / 2)
      basis[, index] <- if ((index - 1L) %% 2L == 1L) {
        cos(pi * frequency * scaled)
      } else {
        sin(pi * frequency * scaled)
      }
    }
  }

  basis <- qr.Q(qr(basis))
  basis <- sweep(basis, 2, sqrt(colMeans(basis^2)), "/")
  colnames(basis) <- paste0("process_basis_", seq_len(rank))
  basis
}

correlated_expected_counts <- function(
  parameter,
  initial_counts,
  mesh,
  latent_z,
  basis
) {
  ipm <- tryCatch(
    build_centered_ipm(parameter, mesh),
    error = function(error) NULL
  )
  if (is.null(ipm)) {
    return(NULL)
  }

  process_sd <- exp(parameter["log_process_sd"])
  marginal_variance <- process_sd^2 * rowSums(basis^2)
  current_intensity <- initial_counts / ipm$delta
  projected_counts <- matrix(
    NA_real_,
    nrow = nrow(latent_z) + 1L,
    ncol = length(mesh)
  )
  projected_counts[1L, ] <- initial_counts

  for (transition in seq_len(nrow(latent_z))) {
    deterministic_intensity <- as.vector(
      ipm$delta * ipm$kernel %*% current_intensity
    )
    multiplier <- exp(
      process_sd * as.vector(basis %*% latent_z[transition, ]) -
        0.5 * marginal_variance
    )
    current_intensity <- deterministic_intensity * multiplier
    projected_counts[transition + 1L, ] <- current_intensity * ipm$delta
  }
  projected_counts
}

correlated_process_log_likelihood <- function(
  parameter,
  counts,
  latent_z,
  mesh,
  basis
) {
  log_likelihood <- 0
  for (index in seq_along(counts)) {
    expected <- correlated_expected_counts(
      parameter,
      counts[[index]][1L, ],
      mesh,
      latent_z[[index]],
      basis
    )
    if (is.null(expected) || any(!is.finite(expected)) || any(expected < 0)) {
      return(-Inf)
    }
    log_likelihood <- log_likelihood + sum(dpois(
      counts[[index]][-1L, , drop = FALSE],
      lambda = pmax(expected[-1L, , drop = FALSE], 1e-12),
      log = TRUE
    ))
  }
  log_likelihood
}

correlated_process_log_prior <- function(parameter, prior) {
  sum(dnorm(parameter, prior$mean, prior$sd, log = TRUE))
}

run_correlated_process_chain <- function(
  counts,
  mesh,
  basis,
  prior,
  initial_parameter,
  proposal_covariance,
  iterations,
  warmup,
  seed
) {
  set.seed(seed)
  parameter_names <- names(initial_parameter)
  demographic_names <- setdiff(parameter_names, "log_process_sd")
  dimension <- length(demographic_names)
  demographic_chol <- chol(
    proposal_covariance[demographic_names, demographic_names, drop = FALSE] +
      diag(1e-10, dimension)
  )

  parameter <- initial_parameter
  latent_z <- lapply(counts, function(trajectory) {
    matrix(0, nrow = nrow(trajectory) - 1L, ncol = ncol(basis))
  })
  log_likelihood <- correlated_process_log_likelihood(
    parameter, counts, latent_z, mesh, basis
  )
  log_prior <- correlated_process_log_prior(parameter, prior)

  draws <- matrix(
    NA_real_,
    nrow = iterations,
    ncol = length(parameter),
    dimnames = list(NULL, parameter_names)
  )
  demographic_accepted <- logical(iterations)
  process_sd_accepted <- logical(iterations)
  latent_accepted <- matrix(
    FALSE,
    nrow = iterations,
    ncol = length(counts)
  )

  log_demographic_scale <- log(2.38 / sqrt(dimension))
  log_process_sd_scale <- log(0.15)
  latent_logit_beta <- qlogis(0.20)

  for (iteration in seq_len(iterations)) {
    demographic_candidate <- parameter
    demographic_candidate[demographic_names] <-
      parameter[demographic_names] +
      as.vector(rnorm(dimension) %*% demographic_chol) *
      exp(log_demographic_scale)
    candidate_log_likelihood <- correlated_process_log_likelihood(
      demographic_candidate, counts, latent_z, mesh, basis
    )
    candidate_log_prior <- correlated_process_log_prior(
      demographic_candidate, prior
    )
    if (
      is.finite(candidate_log_likelihood) &&
        log(runif(1)) <
          candidate_log_likelihood + candidate_log_prior -
          log_likelihood - log_prior
    ) {
      parameter <- demographic_candidate
      log_likelihood <- candidate_log_likelihood
      log_prior <- candidate_log_prior
      demographic_accepted[iteration] <- TRUE
    }

    process_sd_candidate <- parameter
    process_sd_candidate["log_process_sd"] <-
      parameter["log_process_sd"] + rnorm(1, sd = exp(log_process_sd_scale))
    candidate_log_likelihood <- correlated_process_log_likelihood(
      process_sd_candidate, counts, latent_z, mesh, basis
    )
    candidate_log_prior <- correlated_process_log_prior(
      process_sd_candidate, prior
    )
    if (
      is.finite(candidate_log_likelihood) &&
        log(runif(1)) <
          candidate_log_likelihood + candidate_log_prior -
          log_likelihood - log_prior
    ) {
      parameter <- process_sd_candidate
      log_likelihood <- candidate_log_likelihood
      log_prior <- candidate_log_prior
      process_sd_accepted[iteration] <- TRUE
    }

    beta <- plogis(latent_logit_beta)
    for (trajectory_index in seq_along(latent_z)) {
      latent_candidate <- latent_z
      latent_candidate[[trajectory_index]] <-
        sqrt(1 - beta^2) * latent_z[[trajectory_index]] +
        beta * matrix(
          rnorm(length(latent_z[[trajectory_index]])),
          nrow = nrow(latent_z[[trajectory_index]])
        )
      candidate_log_likelihood <- correlated_process_log_likelihood(
        parameter, counts, latent_candidate, mesh, basis
      )
      if (
        is.finite(candidate_log_likelihood) &&
          log(runif(1)) < candidate_log_likelihood - log_likelihood
      ) {
        latent_z <- latent_candidate
        log_likelihood <- candidate_log_likelihood
        latent_accepted[iteration, trajectory_index] <- TRUE
      }
    }

    draws[iteration, ] <- parameter

    if (iteration <= warmup) {
      learning_rate <- min(0.03, 1 / sqrt(iteration))
      log_demographic_scale <- log_demographic_scale +
        learning_rate * (demographic_accepted[iteration] - 0.234)
      log_process_sd_scale <- log_process_sd_scale +
        learning_rate * (process_sd_accepted[iteration] - 0.44)
      latent_logit_beta <- latent_logit_beta +
        learning_rate * (mean(latent_accepted[iteration, ]) - 0.30)
    }
  }

  retained <- (warmup + 1L):iterations
  list(
    draws = draws[retained, , drop = FALSE],
    demographic_acceptance = mean(demographic_accepted[retained]),
    process_sd_acceptance = mean(process_sd_accepted[retained]),
    latent_acceptance = colMeans(latent_accepted[retained, , drop = FALSE])
  )
}

fit_correlated_process_ipm <- function(
  counts,
  mesh,
  prior_regime = "weak",
  basis_rank = 4L,
  chains = 4L,
  iterations = 6000L,
  warmup = 3000L,
  seed = 1L
) {
  counts <- as_profile_trajectory_list(counts)
  basis <- make_smooth_process_basis(mesh, basis_rank)
  demographic_prior <- default_inverse_ipm_priors(prior_regime, "poisson")
  prior <- list(
    mean = c(demographic_prior$mean, log_process_sd = log(0.10)),
    sd = c(demographic_prior$sd, log_process_sd = 1.0),
    regime = prior_regime,
    process_model = "correlated_log_gaussian"
  )

  deterministic_map <- optim(
    par = demographic_prior$mean,
    fn = function(parameter) {
      -inverse_ipm_log_posterior(
        setNames(parameter, names(demographic_prior$mean)),
        counts,
        mesh,
        demographic_prior,
        "poisson"
      )
    },
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 1000, reltol = 1e-9)
  )
  demographic_covariance <- regularized_inverse_hessian(
    deterministic_map$hessian
  )
  rownames(demographic_covariance) <- colnames(demographic_covariance) <-
    names(demographic_prior$mean)
  proposal_covariance <- matrix(
    0,
    nrow = length(prior$mean),
    ncol = length(prior$mean),
    dimnames = list(names(prior$mean), names(prior$mean))
  )
  proposal_covariance[
    names(demographic_prior$mean),
    names(demographic_prior$mean)
  ] <- demographic_covariance
  proposal_covariance["log_process_sd", "log_process_sd"] <- 0.1^2

  initial <- c(
    setNames(deterministic_map$par, names(demographic_prior$mean)),
    log_process_sd = log(0.10)
  )

  chain_results <- lapply(seq_len(chains), function(chain_index) {
    run_correlated_process_chain(
      counts = counts,
      mesh = mesh,
      basis = basis,
      prior = prior,
      initial_parameter = initial,
      proposal_covariance = proposal_covariance,
      iterations = iterations,
      warmup = warmup,
      seed = seed + 1000L * chain_index
    )
  })
  retained_chains <- lapply(chain_results, `[[`, "draws")

  list(
    counts = counts,
    mesh = mesh,
    basis = basis,
    prior = prior,
    process_model = "correlated_log_gaussian",
    chains = retained_chains,
    demographic_acceptance = vapply(
      chain_results, `[[`, numeric(1), "demographic_acceptance"
    ),
    process_sd_acceptance = vapply(
      chain_results, `[[`, numeric(1), "process_sd_acceptance"
    ),
    latent_acceptance = lapply(chain_results, `[[`, "latent_acceptance"),
    rhat = setNames(split_rhat(retained_chains), names(initial)),
    ess = setNames(effective_sample_size(retained_chains), names(initial))
  )
}

posterior_correlated_process_summary <- function(
  fit,
  truth = inverse_ipm_truth()
) {
  posterior_summary(fit, truth)
}

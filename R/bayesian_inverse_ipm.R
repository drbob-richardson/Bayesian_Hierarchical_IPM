# Bayesian profile-only inverse IPM using base-R adaptive Metropolis sampling.

inverse_ipm_truth <- function() {
  c(
    survival_at_80 = -2.2 + 0.032 * 80,
    survival_slope_20 = 0.032 * 20,
    growth_increment_80 = 18 - 0.12 * 80,
    growth_slope_20 = -0.12 * 20,
    log_growth_sd = log(4.5),
    log_recruitment_at_80 = log(0.5) - 7 + 0.075 * 80,
    recruitment_slope_20 = 0.075 * 20,
    recruit_mean_80 = 43 + 0.02 * 80,
    recruit_mean_slope_20 = 0.02 * 20,
    log_recruit_sd = log(3.5)
  )
}

default_inverse_ipm_priors <- function(regime = c(
  "weak", "informative", "informative_wrong"
), process_model = c("poisson", "gamma_poisson")) {
  regime <- match.arg(regime)
  process_model <- match.arg(process_model)
  truth <- inverse_ipm_truth()

  if (regime == "weak") {
    mean <- c(
      survival_at_80 = qlogis(0.60),
      survival_slope_20 = 0.5,
      growth_increment_80 = 8,
      growth_slope_20 = -2,
      log_growth_sd = log(5),
      log_recruitment_at_80 = log(0.20),
      recruitment_slope_20 = 1,
      recruit_mean_80 = 45,
      recruit_mean_slope_20 = 0,
      log_recruit_sd = log(4)
    )
    sd <- c(1.5, 1.0, 6.0, 3.0, 0.7, 1.5, 1.0, 10.0, 2.0, 0.7)
  } else {
    mean <- truth
    sd <- c(0.30, 0.20, 2.0, 1.0, 0.20, 0.30, 0.30, 2.0, 0.50, 0.20)

    if (regime == "informative_wrong") {
      mean[c(
        "survival_at_80",
        "growth_increment_80",
        "log_recruitment_at_80"
      )] <- mean[c(
        "survival_at_80",
        "growth_increment_80",
        "log_recruitment_at_80"
      )] + c(0.75, -3.0, 0.75)
    }
  }

  sd <- setNames(sd, names(mean))
  if (process_model == "gamma_poisson") {
    mean <- c(mean, log_process_precision = log(50))
    sd <- c(sd, log_process_precision = 2)
  }

  list(
    regime = regime,
    process_model = process_model,
    mean = mean,
    sd = sd
  )
}

log_sum_exp <- function(value) {
  maximum <- max(value)
  maximum + log(sum(exp(value - maximum)))
}

log_multivariate_normal <- function(value, mean, covariance) {
  centered <- value - mean
  decomposition <- chol(covariance)
  standardized <- backsolve(decomposition, centered, transpose = TRUE)
  -0.5 * (
    length(value) * log(2 * pi) +
      2 * sum(log(diag(decomposition))) +
      sum(standardized^2)
  )
}

inverse_ipm_log_prior <- function(parameter, prior) {
  block_parameters <- character()
  if (!is.null(prior$blocks)) {
    block_parameters <- unique(unlist(lapply(
      prior$blocks,
      `[[`,
      "parameters"
    )))
  }

  independent <- setdiff(names(parameter), block_parameters)
  log_prior <- sum(dnorm(
    parameter[independent],
    mean = prior$mean[independent],
    sd = prior$sd[independent],
    log = TRUE
  ))

  if (is.null(prior$blocks)) {
    return(log_prior)
  }

  for (block in prior$blocks) {
    value <- parameter[block$parameters]
    external_density <- log_multivariate_normal(
      value,
      block$mean,
      block$covariance
    )
    if (is.null(block$robust_weight)) {
      log_prior <- log_prior + external_density
    } else {
      weak_density <- log_multivariate_normal(
        value,
        block$weak_mean,
        block$weak_covariance
      )
      log_prior <- log_prior + log_sum_exp(c(
        log(block$robust_weight) + external_density,
        log1p(-block$robust_weight) + weak_density
      ))
    }
  }
  log_prior
}

build_centered_ipm <- function(parameter, mesh) {
  required <- names(inverse_ipm_truth())
  stopifnot(all(required %in% names(parameter)))

  z <- function(size) (size - 80) / 20

  make_ipm_kernel(
    mesh = mesh,
    survival = function(size) {
      inv_logit(
        parameter["survival_at_80"] +
          parameter["survival_slope_20"] * z(size)
      )
    },
    growth_mean = function(size) {
      pmax(
        size,
        size + parameter["growth_increment_80"] +
          parameter["growth_slope_20"] * z(size)
      )
    },
    growth_sd = function(size) exp(parameter["log_growth_sd"]),
    recruitment = function(size) {
      exp(
        parameter["log_recruitment_at_80"] +
          parameter["recruitment_slope_20"] * z(size)
      )
    },
    recruit_mean = function(size) {
      parameter["recruit_mean_80"] +
        parameter["recruit_mean_slope_20"] * z(size)
    },
    recruit_sd = function(size) exp(parameter["log_recruit_sd"])
  )
}

make_profile_count_matrix <- function(census, mesh, transitions = NULL) {
  delta <- mean(diff(mesh))
  breaks <- c(mesh - delta / 2, tail(mesh, 1) + delta / 2)
  years <- sort(unique(census$year))

  if (!is.null(transitions)) {
    years <- years[years <= min(years) + transitions]
  }

  counts <- vapply(years, function(year) {
    hist(
      census$size[census$year == year],
      breaks = breaks,
      plot = FALSE
    )$counts
  }, integer(length(mesh)))

  counts <- t(counts)
  rownames(counts) <- years
  colnames(counts) <- mesh
  counts
}

expected_profile_counts <- function(parameter, initial_counts, mesh, transitions) {
  ipm <- tryCatch(
    build_centered_ipm(parameter, mesh),
    error = function(error) NULL
  )
  if (is.null(ipm)) {
    return(NULL)
  }

  initial_intensity <- initial_counts / ipm$delta
  projected <- project_ipm(initial_intensity, ipm, transitions = transitions)
  projected * ipm$delta
}

as_profile_trajectory_list <- function(counts) {
  if (is.matrix(counts)) {
    counts <- list(trajectory_1 = counts)
  }
  stopifnot(
    is.list(counts),
    length(counts) >= 1L,
    all(vapply(counts, is.matrix, logical(1))),
    all(vapply(counts, nrow, integer(1)) >= 2L)
  )
  counts
}

inverse_ipm_log_posterior <- function(
  parameter,
  counts,
  mesh,
  prior,
  process_model = c("poisson", "gamma_poisson")
) {
  process_model <- match.arg(process_model)
  counts <- as_profile_trajectory_list(counts)
  if (any(!is.finite(parameter))) {
    return(-Inf)
  }

  log_likelihood <- 0
  for (trajectory in counts) {
    expected <- expected_profile_counts(
      parameter,
      initial_counts = trajectory[1L, ],
      mesh = mesh,
      transitions = nrow(trajectory) - 1L
    )
    if (is.null(expected) || any(!is.finite(expected)) || any(expected < 0)) {
      return(-Inf)
    }
    expected <- pmax(expected[-1L, , drop = FALSE], 1e-12)
    observed <- trajectory[-1L, , drop = FALSE]

    if (process_model == "poisson") {
      log_likelihood <- log_likelihood +
        sum(dpois(observed, lambda = expected, log = TRUE))
    } else {
      process_precision <- exp(parameter["log_process_precision"])
      if (!is.finite(process_precision) || process_precision <= 0) {
        return(-Inf)
      }
      log_likelihood <- log_likelihood + sum(dnbinom(
        observed,
        mu = expected,
        size = process_precision,
        log = TRUE
      ))
    }
  }
  log_prior <- inverse_ipm_log_prior(parameter, prior)

  log_likelihood + log_prior
}

regularized_inverse_hessian <- function(hessian, minimum_eigenvalue = 1e-5) {
  information <- (hessian + t(hessian)) / 2
  decomposition <- eigen(information, symmetric = TRUE)
  values <- pmax(decomposition$values, minimum_eigenvalue)
  covariance <- decomposition$vectors %*%
    diag(1 / values, nrow = length(values)) %*%
    t(decomposition$vectors)
  (covariance + t(covariance)) / 2
}

run_adaptive_metropolis <- function(
  log_posterior,
  initial,
  proposal_covariance,
  iterations,
  warmup,
  seed,
  target_acceptance = 0.234
) {
  set.seed(seed)
  dimension <- length(initial)
  draws <- matrix(NA_real_, nrow = iterations, ncol = dimension)
  colnames(draws) <- names(initial)

  current <- initial
  current_log_posterior <- log_posterior(current)
  if (!is.finite(current_log_posterior)) {
    stop("Initial chain state has non-finite log posterior.")
  }

  base_chol <- chol(proposal_covariance + diag(1e-10, dimension))
  log_scale <- log(2.38 / sqrt(dimension))
  accepted <- logical(iterations)

  for (iteration in seq_len(iterations)) {
    candidate <- current + as.vector(
      rnorm(dimension) %*% base_chol
    ) * exp(log_scale)
    names(candidate) <- names(initial)
    candidate_log_posterior <- log_posterior(candidate)

    if (
      is.finite(candidate_log_posterior) &&
        log(runif(1)) < candidate_log_posterior - current_log_posterior
    ) {
      current <- candidate
      current_log_posterior <- candidate_log_posterior
      accepted[iteration] <- TRUE
    }
    draws[iteration, ] <- current

    if (iteration <= warmup) {
      learning_rate <- min(0.05, 1 / sqrt(iteration))
      log_scale <- log_scale +
        learning_rate * (accepted[iteration] - target_acceptance)
    }
  }

  list(
    draws = draws[(warmup + 1L):iterations, , drop = FALSE],
    acceptance = mean(accepted[(warmup + 1L):iterations]),
    warmup_acceptance = mean(accepted[seq_len(warmup)]),
    final_scale = exp(log_scale)
  )
}

split_rhat <- function(chains) {
  chain_length <- min(vapply(chains, nrow, integer(1)))
  half <- floor(chain_length / 2)
  split_chains <- unlist(lapply(chains, function(chain) {
    list(
      chain[seq_len(half), , drop = FALSE],
      chain[(chain_length - half + 1L):chain_length, , drop = FALSE]
    )
  }), recursive = FALSE)

  means <- vapply(split_chains, colMeans, numeric(ncol(split_chains[[1L]])))
  variances <- vapply(
    split_chains,
    function(chain) apply(chain, 2, var),
    numeric(ncol(split_chains[[1L]]))
  )
  between <- half * apply(means, 1, var)
  within <- rowMeans(variances)
  variance_hat <- ((half - 1) / half) * within + between / half
  sqrt(variance_hat / within)
}

effective_sample_size <- function(chains, max_lag = 500L) {
  combined <- do.call(rbind, chains)
  n_total <- nrow(combined)

  vapply(seq_len(ncol(combined)), function(index) {
    acf_values <- acf(
      combined[, index],
      lag.max = min(max_lag, n_total - 1L),
      plot = FALSE,
      demean = TRUE
    )$acf[-1L]
    first_negative <- which(acf_values < 0)[1L]
    if (is.na(first_negative)) {
      positive <- acf_values
    } else {
      positive <- acf_values[seq_len(max(1L, first_negative - 1L))]
    }
    n_total / (1 + 2 * sum(positive))
  }, numeric(1))
}

fit_bayesian_inverse_ipm <- function(
  counts,
  mesh,
  prior = default_inverse_ipm_priors("weak"),
  process_model = prior$process_model,
  chains = 4L,
  iterations = 6000L,
  warmup = 3000L,
  seed = 1L
) {
  counts <- as_profile_trajectory_list(counts)
  stopifnot(iterations > warmup, chains >= 2L)
  process_model <- match.arg(process_model, c("poisson", "gamma_poisson"))

  log_posterior <- function(parameter) {
    inverse_ipm_log_posterior(parameter, counts, mesh, prior, process_model)
  }

  map_fit <- optim(
    par = prior$mean,
    fn = function(parameter) -log_posterior(setNames(parameter, names(prior$mean))),
    method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-9),
    hessian = TRUE
  )
  map <- setNames(map_fit$par, names(prior$mean))
  proposal_covariance <- regularized_inverse_hessian(map_fit$hessian)

  chain_results <- lapply(seq_len(chains), function(chain_index) {
    set.seed(seed + chain_index)
    initial <- map
    for (attempt in seq_len(20L)) {
      jitter_scale <- 0.25 / sqrt(attempt)
      candidate <- as.vector(map + t(chol(proposal_covariance)) %*%
        rnorm(length(map), sd = jitter_scale))
      names(candidate) <- names(map)
      if (is.finite(log_posterior(candidate))) {
        initial <- candidate
        break
      }
    }

    run_adaptive_metropolis(
      log_posterior = log_posterior,
      initial = initial,
      proposal_covariance = proposal_covariance,
      iterations = iterations,
      warmup = warmup,
      seed = seed + 1000L * chain_index
    )
  })

  retained_chains <- lapply(chain_results, `[[`, "draws")
  names(retained_chains) <- paste0("chain_", seq_len(chains))

  list(
    counts = counts,
    mesh = mesh,
    prior = prior,
    process_model = process_model,
    map = map,
    map_convergence = map_fit$convergence,
    proposal_covariance = proposal_covariance,
    chains = retained_chains,
    acceptance = vapply(chain_results, `[[`, numeric(1), "acceptance"),
    warmup_acceptance = vapply(
      chain_results,
      `[[`,
      numeric(1),
      "warmup_acceptance"
    ),
    rhat = setNames(split_rhat(retained_chains), names(map)),
    ess = setNames(effective_sample_size(retained_chains), names(map))
  )
}

posterior_summary <- function(fit, truth = inverse_ipm_truth()) {
  draws <- do.call(rbind, fit$chains)
  quantiles <- t(apply(draws, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  posterior_sd <- apply(draws, 2, sd)

  aligned_truth <- truth[colnames(draws)]
  finite_truth <- is.finite(aligned_truth)
  coverage <- rep(NA, length(aligned_truth))
  coverage[finite_truth] <- quantiles[finite_truth, 1L] <=
    aligned_truth[finite_truth] &
    quantiles[finite_truth, 3L] >= aligned_truth[finite_truth]
  standardized_bias <- rep(NA_real_, length(aligned_truth))
  standardized_bias[finite_truth] <- (
    colMeans(draws)[finite_truth] - aligned_truth[finite_truth]
  ) / posterior_sd[finite_truth]

  data.frame(
    parameter = colnames(draws),
    truth = aligned_truth,
    prior_mean = fit$prior$mean[colnames(draws)],
    prior_sd = fit$prior$sd[colnames(draws)],
    posterior_mean = colMeans(draws),
    posterior_sd = posterior_sd,
    q025 = quantiles[, 1L],
    median = quantiles[, 2L],
    q975 = quantiles[, 3L],
    covers_truth = coverage,
    contraction = 1 - posterior_sd / fit$prior$sd[colnames(draws)],
    standardized_bias = standardized_bias,
    rhat = fit$rhat[colnames(draws)],
    ess = fit$ess[colnames(draws)],
    row.names = NULL
  )
}

derived_inverse_ipm_quantities <- function(
  parameter,
  mesh,
  evaluation_sizes = c(50, 80, 110)
) {
  z <- function(size) (size - 80) / 20
  survival <- inv_logit(
    parameter["survival_at_80"] +
      parameter["survival_slope_20"] * z(evaluation_sizes)
  )
  growth_increment <- pmax(
    0,
    parameter["growth_increment_80"] +
      parameter["growth_slope_20"] * z(evaluation_sizes)
  )
  recruitment <- exp(
    parameter["log_recruitment_at_80"] +
      parameter["recruitment_slope_20"] * z(evaluation_sizes)
  )
  recruit_mean <- parameter["recruit_mean_80"] +
    parameter["recruit_mean_slope_20"] * z(evaluation_sizes)

  ipm <- build_centered_ipm(parameter, mesh)
  dominant_eigenvalue <- max(Re(eigen(
    ipm$delta * ipm$kernel,
    only.values = TRUE
  )$values))

  c(
    setNames(survival, paste0("survival_", evaluation_sizes)),
    setNames(growth_increment, paste0("growth_increment_", evaluation_sizes)),
    growth_sd = exp(parameter["log_growth_sd"]),
    setNames(recruitment, paste0("recruitment_", evaluation_sizes)),
    setNames(recruit_mean, paste0("recruit_mean_", evaluation_sizes)),
    recruit_sd = exp(parameter["log_recruit_sd"]),
    lambda = dominant_eigenvalue
  )
}

posterior_derived_summary <- function(
  fit,
  truth = inverse_ipm_truth(),
  evaluation_sizes = c(50, 80, 110),
  max_draws = 2000L
) {
  draws <- do.call(rbind, fit$chains)
  if (nrow(draws) > max_draws) {
    index <- unique(round(seq(1, nrow(draws), length.out = max_draws)))
    draws <- draws[index, , drop = FALSE]
  }

  derived_draws <- t(apply(draws, 1, function(parameter) {
    derived_inverse_ipm_quantities(
      setNames(parameter, colnames(draws)),
      fit$mesh,
      evaluation_sizes
    )
  }))
  truth_derived <- derived_inverse_ipm_quantities(
    truth,
    fit$mesh,
    evaluation_sizes
  )
  quantiles <- t(apply(
    derived_draws,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975)
  ))

  data.frame(
    quantity = colnames(derived_draws),
    truth = truth_derived[colnames(derived_draws)],
    posterior_mean = colMeans(derived_draws),
    posterior_sd = apply(derived_draws, 2, sd),
    q025 = quantiles[, 1L],
    median = quantiles[, 2L],
    q975 = quantiles[, 3L],
    covers_truth = quantiles[, 1L] <= truth_derived[colnames(derived_draws)] &
      quantiles[, 3L] >= truth_derived[colnames(derived_draws)],
    row.names = NULL
  )
}

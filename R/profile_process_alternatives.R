# Alternative profile-process likelihoods that preserve a positive IPM kernel.

one_step_expected_counts <- function(
  parameter,
  source_counts,
  mesh,
  kernel_builder = build_centered_ipm
) {
  ipm <- tryCatch(
    kernel_builder(parameter, mesh),
    error = function(error) NULL
  )
  if (is.null(ipm)) {
    return(NULL)
  }
  as.vector(ipm$delta * ipm$kernel %*% source_counts)
}

shared_gamma_poisson_log_likelihood <- function(
  observed,
  expected,
  precision
) {
  total_observed <- sum(observed)
  total_expected <- sum(expected)
  if (
    !is.finite(precision) ||
      precision <= 0 ||
      !is.finite(total_expected) ||
      total_expected <= 0 ||
      any(!is.finite(expected)) ||
      any(expected <= 0)
  ) {
    return(-Inf)
  }

  lgamma(precision + total_observed) -
    lgamma(precision) +
    precision * log(precision) -
    (precision + total_observed) * log(precision + total_expected) +
    sum(observed * log(expected) - lgamma(observed + 1))
}

multivariate_normal_log_likelihood <- function(
  observed,
  expected,
  covariance,
  minimum_variance = 0.25
) {
  if (
    any(!is.finite(observed)) ||
      any(!is.finite(expected)) ||
      any(!is.finite(covariance)) ||
      any(dim(covariance) != c(length(observed), length(observed)))
  ) {
    return(-Inf)
  }

  covariance <- (covariance + t(covariance)) / 2
  diagonal <- diag(covariance)
  covariance[cbind(seq_along(diagonal), seq_along(diagonal))] <-
    pmax(diagonal, minimum_variance)
  decomposition <- tryCatch(
    chol(covariance),
    error = function(error) NULL
  )
  if (is.null(decomposition)) {
    minimum_eigenvalue <- min(eigen(covariance, symmetric = TRUE)$values)
    covariance <- covariance +
      diag(
        max(minimum_variance, -minimum_eigenvalue + minimum_variance),
        nrow = length(observed)
      )
    decomposition <- tryCatch(
      chol(covariance),
      error = function(error) NULL
    )
    if (is.null(decomposition)) {
      return(-Inf)
    }
  }

  residual <- observed - expected
  standardized <- backsolve(
    decomposition,
    residual,
    transpose = TRUE
  )
  -0.5 * (
    length(observed) * log(2 * pi) +
      2 * sum(log(diag(decomposition))) +
      sum(standardized^2)
  )
}

finite_population_transition_moments <- function(
  source_counts,
  ipm,
  detection_probability = 1,
  residual_precision = Inf,
  fecundity_cv2 = 1
) {
  if (
    length(detection_probability) != 1L ||
      !is.finite(detection_probability) ||
      detection_probability <= 0 ||
      detection_probability > 1 ||
      !is.finite(fecundity_cv2) ||
      fecundity_cv2 < 0 ||
      any(!is.finite(source_counts)) ||
      any(source_counts < 0)
  ) {
    return(NULL)
  }

  projection_probability <- ipm$delta * ipm$projection
  fertility_mean <- ipm$delta * ipm$fertility
  transition_mean <- projection_probability + fertility_mean

  estimated_source <- source_counts / detection_probability
  source_variance <- source_counts * (1 - detection_probability) /
    detection_probability^2
  latent_mean <- as.vector(transition_mean %*% estimated_source)

  # A survivor enters one size bin or dies, producing negative cross-bin
  # covariance. Poisson offspring contribute independent variance by bin.
  process_covariance <- diag(latent_mean, nrow = length(latent_mean)) -
    sweep(
      projection_probability,
      2L,
      estimated_source,
      "*"
    ) %*% t(projection_probability)

  # Marginalizing unobserved reproductive heterogeneity adds positive
  # covariance among offspring bins. A 1:1 mixture of non-reproductive
  # males and females with twice the population-average fecundity has CV^2=1.
  process_covariance <- process_covariance +
    fecundity_cv2 * sweep(
      fertility_mean,
      2L,
      estimated_source,
      "*"
    ) %*% t(fertility_mean)

  # Propagate uncertainty about the unobserved members of an incompletely
  # detected source profile through the transition.
  source_covariance <- sweep(
    transition_mean,
    2L,
    source_variance,
    "*"
  ) %*% t(transition_mean)
  latent_covariance <- process_covariance + source_covariance

  observed_mean <- detection_probability * latent_mean
  observed_covariance <- detection_probability^2 * latent_covariance +
    diag(
      detection_probability * (1 - detection_probability) * latent_mean,
      nrow = length(latent_mean)
    )

  if (is.finite(residual_precision)) {
    if (residual_precision <= 0) {
      return(NULL)
    }
    observed_covariance <- observed_covariance +
      diag(observed_mean^2 / residual_precision, nrow = length(observed_mean))
  }

  list(
    mean = observed_mean,
    covariance = observed_covariance
  )
}

one_step_profile_log_posterior <- function(
  parameter,
  counts,
  mesh,
  prior,
  process_model = c(
    "poisson",
    "gamma_poisson",
    "shared_gamma_poisson",
    "finite_population",
    "finite_population_overdispersed"
  ),
  kernel_builder = build_centered_ipm,
  detection_probability = 1
) {
  process_model <- match.arg(process_model)
  counts <- as_profile_trajectory_list(counts)
  if (any(!is.finite(parameter))) {
    return(-Inf)
  }
  ipm <- tryCatch(
    kernel_builder(parameter, mesh),
    error = function(error) NULL
  )
  if (is.null(ipm)) {
    return(-Inf)
  }

  log_likelihood <- 0
  for (trajectory in counts) {
    for (transition in seq_len(nrow(trajectory) - 1L)) {
      source_counts <- trajectory[transition, ]
      expected <- as.vector(ipm$delta * ipm$kernel %*% source_counts)
      if (
        is.null(expected) ||
          any(!is.finite(expected)) ||
          any(expected < 0)
      ) {
        return(-Inf)
      }
      expected <- pmax(expected, 1e-12)
      observed <- trajectory[transition + 1L, ]

      if (process_model == "poisson") {
        log_likelihood <- log_likelihood + sum(dpois(
          observed,
          lambda = expected,
          log = TRUE
        ))
      } else if (process_model == "gamma_poisson") {
        precision <- exp(parameter["log_process_precision"])
        if (!is.finite(precision) || precision <= 0) {
          return(-Inf)
        }
        log_likelihood <- log_likelihood + sum(dnbinom(
          observed,
          mu = expected,
          size = precision,
          log = TRUE
        ))
      } else if (process_model == "shared_gamma_poisson") {
        precision <- exp(parameter["log_process_precision"])
        log_likelihood <- log_likelihood +
          shared_gamma_poisson_log_likelihood(
            observed,
            expected,
            precision
          )
      } else {
        residual_precision <- if (
          process_model == "finite_population_overdispersed"
        ) {
          exp(parameter["log_process_precision"])
        } else {
          Inf
        }
        moments <- finite_population_transition_moments(
          source_counts,
          ipm,
          detection_probability,
          residual_precision
        )
        if (is.null(moments)) {
          return(-Inf)
        }
        log_likelihood <- log_likelihood +
          multivariate_normal_log_likelihood(
            observed,
            moments$mean,
            moments$covariance
          )
      }
    }
  }

  log_likelihood + inverse_ipm_log_prior(parameter, prior)
}

profile_process_method <- function(method) {
  choices <- c(
    "recursive_poisson",
    "recursive_gamma_poisson",
    "one_step_poisson",
    "one_step_gamma_poisson",
    "one_step_shared_gamma_poisson",
    "one_step_finite_population",
    "one_step_finite_population_overdispersed",
    "one_step_gamma_poisson_shape_informed",
    "one_step_gamma_poisson_fixed_recruitment_slope"
  )
  match.arg(method, choices)
}

profile_process_base_method <- function(method) {
  method <- profile_process_method(method)
  if (method %in% c(
    "one_step_gamma_poisson_shape_informed",
    "one_step_gamma_poisson_fixed_recruitment_slope"
  )) {
    return("one_step_gamma_poisson")
  }
  method
}

profile_process_fixed_parameters <- function(method) {
  method <- profile_process_method(method)
  if (method == "one_step_gamma_poisson_fixed_recruitment_slope") {
    return(inverse_ipm_truth()["recruitment_slope_20"])
  }
  numeric()
}

expand_profile_process_parameter <- function(parameter, fixed) {
  if (length(fixed) == 0L) {
    return(parameter)
  }
  full <- c(parameter, fixed)
  ordered <- c(
    names(inverse_ipm_truth()),
    setdiff(names(full), names(inverse_ipm_truth()))
  )
  full[ordered]
}

remove_fixed_from_prior <- function(prior, fixed) {
  if (length(fixed) == 0L) {
    return(prior)
  }
  keep <- setdiff(names(prior$mean), names(fixed))
  prior$mean <- prior$mean[keep]
  prior$sd <- prior$sd[keep]
  if (!is.null(prior$alternative_starts)) {
    prior$alternative_starts <- lapply(
      prior$alternative_starts,
      function(start) start[keep]
    )
  }
  prior
}

profile_process_prior <- function(
  method,
  regime = c("weak", "informative", "informative_wrong")
) {
  method <- profile_process_method(method)
  regime <- match.arg(regime)
  if (profile_process_base_method(method) %in% c(
    "recursive_gamma_poisson",
    "one_step_gamma_poisson",
    "one_step_shared_gamma_poisson",
    "one_step_finite_population_overdispersed"
  )) {
    prior <- default_inverse_ipm_priors(regime, "gamma_poisson")
  } else {
    prior <- default_inverse_ipm_priors(regime, "poisson")
  }
  if (method == "one_step_gamma_poisson_shape_informed") {
    prior$mean["recruitment_slope_20"] <-
      inverse_ipm_truth()["recruitment_slope_20"]
    prior$sd["recruitment_slope_20"] <- 0.20
  }
  prior
}

profile_process_log_posterior <- function(
  parameter,
  counts,
  mesh,
  prior,
  method,
  kernel_builder = build_centered_ipm,
  detection_probability = 1
) {
  method <- profile_process_base_method(method)
  if (method == "recursive_poisson") {
    return(inverse_ipm_log_posterior(
      parameter,
      counts,
      mesh,
      prior,
      "poisson"
    ))
  }
  if (method == "recursive_gamma_poisson") {
    return(inverse_ipm_log_posterior(
      parameter,
      counts,
      mesh,
      prior,
      "gamma_poisson"
    ))
  }

  one_step_profile_log_posterior(
    parameter,
    counts,
    mesh,
    prior,
    sub("^one_step_", "", method),
    kernel_builder,
    detection_probability
  )
}

profile_process_parameter_group <- function(quantity) {
  if (quantity == "log_process_precision") {
    return("process")
  }
  parameter_group(quantity)
}

# Biologically interpretable prior-information regimes.

inverse_ipm_information_groups <- function() {
  list(
    survival_growth = c(
      "survival_at_80",
      "survival_slope_20",
      "growth_increment_80",
      "growth_slope_20",
      "log_growth_sd"
    ),
    fecundity = c(
      "log_recruitment_at_80",
      "recruitment_slope_20"
    ),
    recruit_size = c(
      "recruit_mean_80",
      "recruit_mean_slope_20",
      "log_recruit_sd"
    )
  )
}

make_external_information_prior <- function(
  regime = c(
    "weak",
    "survival_growth",
    "recruit_size",
    "fecundity",
    "all_vital_rates",
    "all_vital_rates_biased"
  ),
  process_model = c("poisson", "gamma_poisson"),
  uncertainty_multiplier = 2
) {
  regime <- match.arg(regime)
  process_model <- match.arg(process_model)
  stopifnot(
    length(uncertainty_multiplier) == 1L,
    is.finite(uncertainty_multiplier),
    uncertainty_multiplier > 0
  )

  prior <- default_inverse_ipm_priors("weak", process_model)
  prior$regime <- regime
  prior$information_source <- "diffuse"
  prior$uncertainty_multiplier <- uncertainty_multiplier

  if (regime == "weak") {
    return(prior)
  }

  truth <- inverse_ipm_truth()
  oracle <- default_inverse_ipm_priors("informative", process_model)
  groups <- inverse_ipm_information_groups()
  selected <- switch(
    regime,
    survival_growth = groups$survival_growth,
    recruit_size = groups$recruit_size,
    fecundity = groups$fecundity,
    all_vital_rates = names(truth),
    all_vital_rates_biased = names(truth)
  )

  prior$mean[selected] <- truth[selected]
  prior$sd[selected] <- uncertainty_multiplier * oracle$sd[selected]
  prior$information_source <- paste(selected, collapse = ",")

  if (regime == "all_vital_rates_biased") {
    # A related-population prior: broadly right, but systematically displaced.
    bias <- c(
      survival_at_80 = 0.50,
      survival_slope_20 = -0.20,
      growth_increment_80 = -2.00,
      growth_slope_20 = 0,
      log_growth_sd = 0.15,
      log_recruitment_at_80 = 0.50,
      recruitment_slope_20 = -0.30,
      recruit_mean_80 = 3.00,
      recruit_mean_slope_20 = 0.20,
      log_recruit_sd = 0.15
    )
    prior$mean[names(bias)] <- prior$mean[names(bias)] + bias
    prior$bias <- bias
  }

  prior
}

prior_effective_sample_size_normal <- function(prior_sd, observation_sd) {
  stopifnot(
    length(prior_sd) == length(observation_sd),
    all(prior_sd > 0),
    all(observation_sd > 0)
  )
  (observation_sd / prior_sd)^2
}

make_fecundity_study_prior <- function(
  sample_size,
  accuracy = c("correct", "biased"),
  borrowing = c("full", "robust"),
  process_model = c("poisson", "gamma_poisson"),
  robust_weight = 0.8,
  truth = inverse_ipm_truth(),
  bias = NULL,
  seed = 1L
) {
  accuracy <- match.arg(accuracy)
  borrowing <- match.arg(borrowing)
  process_model <- match.arg(process_model)
  stopifnot(
    sample_size >= 2L,
    robust_weight > 0,
    robust_weight < 1
  )

  prior <- default_inverse_ipm_priors("weak", process_model)
  parameters <- c("log_recruitment_at_80", "recruitment_slope_20")
  weak_mean <- prior$mean[parameters]
  weak_covariance <- diag(prior$sd[parameters]^2)
  dimnames(weak_covariance) <- list(parameters, parameters)

  truth <- truth[parameters]
  study_truth <- truth
  if (accuracy == "biased") {
    if (is.null(bias)) {
      bias <- c(
        log_recruitment_at_80 = 0.50,
        recruitment_slope_20 = -0.30
      )
    }
    stopifnot(all(parameters %in% names(bias)))
    study_truth <- study_truth + bias[parameters]
  }

  set.seed(seed)
  z <- seq(-1.5, 1.5, length.out = sample_size)
  expected <- exp(study_truth[1L] + study_truth[2L] * z)
  offspring <- rpois(sample_size, expected)
  design <- cbind(1, z)

  log_external_posterior <- function(coefficient) {
    linear_predictor <- as.vector(design %*% coefficient)
    sum(dpois(offspring, exp(linear_predictor), log = TRUE)) +
      log_multivariate_normal(coefficient, weak_mean, weak_covariance)
  }
  fit <- optim(
    par = weak_mean,
    fn = function(coefficient) -log_external_posterior(coefficient),
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 1000L, reltol = 1e-10)
  )
  external_mean <- setNames(fit$par, parameters)
  external_covariance <- regularized_inverse_hessian(fit$hessian, 1e-7)
  dimnames(external_covariance) <- list(parameters, parameters)

  block <- list(
    parameters = parameters,
    mean = external_mean,
    covariance = external_covariance,
    robust_weight = if (borrowing == "robust") robust_weight else NULL,
    weak_mean = weak_mean,
    weak_covariance = weak_covariance
  )
  prior$blocks <- list(fecundity_study = block)
  prior$regime <- paste("fecundity_study", accuracy, borrowing, sep = "_")
  prior$external_study <- list(
    sample_size = sample_size,
    accuracy = accuracy,
    borrowing = borrowing,
    z = z,
    offspring = offspring,
    truth = study_truth,
    posterior_mean = external_mean,
    posterior_covariance = external_covariance
  )

  if (borrowing == "full") {
    prior$mean[parameters] <- external_mean
    prior$sd[parameters] <- sqrt(diag(external_covariance))
  } else {
    mixture_mean <- robust_weight * external_mean +
      (1 - robust_weight) * weak_mean
    mean_difference <- external_mean - weak_mean
    mixture_covariance <- robust_weight * external_covariance +
      (1 - robust_weight) * weak_covariance +
      robust_weight * (1 - robust_weight) *
        tcrossprod(mean_difference)
    prior$mean[parameters] <- mixture_mean
    prior$sd[parameters] <- sqrt(diag(mixture_covariance))
    external_start <- prior$mean
    external_start[parameters] <- external_mean
    weak_start <- prior$mean
    weak_start[parameters] <- weak_mean
    prior$alternative_starts <- list(external_start, weak_start)
  }

  prior
}

external_block_responsibility <- function(parameter, block) {
  stopifnot(!is.null(block$robust_weight))
  value <- parameter[block$parameters]
  external <- log(block$robust_weight) + log_multivariate_normal(
    value,
    block$mean,
    block$covariance
  )
  weak <- log1p(-block$robust_weight) + log_multivariate_normal(
    value,
    block$weak_mean,
    block$weak_covariance
  )
  exp(external - log_sum_exp(c(external, weak)))
}

robust_prior_component <- function(prior, component = c("external", "weak")) {
  component <- match.arg(component)
  stopifnot(
    length(prior$blocks) == 1L,
    !is.null(prior$blocks[[1L]]$robust_weight)
  )

  result <- prior
  block <- result$blocks[[1L]]
  if (component == "weak") {
    block$mean <- block$weak_mean
    block$covariance <- block$weak_covariance
  }
  block$robust_weight <- NULL
  result$blocks[[1L]] <- block
  result$mean[block$parameters] <- block$mean
  result$sd[block$parameters] <- sqrt(diag(block$covariance))
  result$alternative_starts <- NULL
  result$component <- component
  result
}

laplace_log_evidence <- function(fit) {
  dimension <- length(fit$map)
  log_determinant <- as.numeric(determinant(
    fit$covariance,
    logarithm = TRUE
  )$modulus)
  fit$log_posterior_at_map +
    0.5 * dimension * log(2 * pi) +
    0.5 * log_determinant
}

fit_robust_laplace_inverse_ipm <- function(
  counts,
  mesh,
  prior,
  process_model = prior$process_model,
  maxit = 1000L
) {
  stopifnot(
    length(prior$blocks) == 1L,
    !is.null(prior$blocks[[1L]]$robust_weight)
  )
  mixture_prior_weight <- prior$blocks[[1L]]$robust_weight
  external_prior <- robust_prior_component(prior, "external")
  weak_prior <- robust_prior_component(prior, "weak")

  external_fit <- fit_laplace_inverse_ipm(
    counts,
    mesh,
    external_prior,
    process_model,
    maxit
  )
  weak_fit <- fit_laplace_inverse_ipm(
    counts,
    mesh,
    weak_prior,
    process_model,
    maxit
  )
  log_weight <- c(
    external = log(mixture_prior_weight) + laplace_log_evidence(external_fit),
    weak = log1p(-mixture_prior_weight) + laplace_log_evidence(weak_fit)
  )
  posterior_external_weight <- exp(
    log_weight["external"] - log_sum_exp(log_weight)
  )

  list(
    external_fit = external_fit,
    weak_fit = weak_fit,
    posterior_external_weight = posterior_external_weight,
    mixture_prior_weight = mixture_prior_weight,
    positive_definite = external_fit$positive_definite &&
      weak_fit$positive_definite,
    condition_number = max(
      external_fit$condition_number,
      weak_fit$condition_number
    ),
    convergence = max(external_fit$convergence, weak_fit$convergence),
    mesh = mesh,
    process_model = process_model
  )
}

sample_robust_laplace_draws <- function(
  fit,
  draws = 800L,
  seed = 1L
) {
  set.seed(seed)
  external_draws <- rbinom(1L, draws, fit$posterior_external_weight)
  weak_draws <- draws - external_draws
  samples <- list()
  if (external_draws > 0L) {
    samples[[length(samples) + 1L]] <- sample_laplace_draws(
      fit$external_fit,
      external_draws,
      seed + 1L
    )
  }
  if (weak_draws > 0L) {
    samples[[length(samples) + 1L]] <- sample_laplace_draws(
      fit$weak_fit,
      weak_draws,
      seed + 2L
    )
  }
  result <- do.call(rbind, samples)
  result[sample.int(nrow(result)), , drop = FALSE]
}

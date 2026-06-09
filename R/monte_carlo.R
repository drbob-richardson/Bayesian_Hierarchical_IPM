# Fast, resumable Monte Carlo tools using Laplace posterior approximations.

fit_laplace_inverse_ipm <- function(
  counts,
  mesh,
  prior = default_inverse_ipm_priors("weak"),
  process_model = prior$process_model,
  maxit = 1000L
) {
  counts <- as_profile_trajectory_list(counts)
  process_model <- match.arg(process_model, c("poisson", "gamma_poisson"))

  log_posterior <- function(parameter) {
    inverse_ipm_log_posterior(
      setNames(parameter, names(prior$mean)),
      counts,
      mesh,
      prior,
      process_model
    )
  }

  starts <- c(list(prior$mean), prior$alternative_starts)
  if (process_model == "gamma_poisson") {
    starts <- c(
      starts,
      lapply(c(log(5), log(25), log(100), log(500)), function(value) {
        start <- prior$mean
        start["log_process_precision"] <- value
        start
      })
    )
  }

  fits <- lapply(starts, function(start) {
    tryCatch(
      optim(
        par = start,
        fn = function(parameter) -log_posterior(parameter),
        method = "BFGS",
        control = list(maxit = maxit, reltol = 1e-9),
        hessian = TRUE
      ),
      error = function(error) NULL
    )
  })
  valid <- vapply(fits, function(fit) {
    !is.null(fit) && is.finite(fit$value) && all(is.finite(fit$par))
  }, logical(1))
  if (!any(valid)) {
    stop("No finite Laplace posterior mode was found.")
  }
  fits <- fits[valid]
  fit <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]

  information <- (fit$hessian + t(fit$hessian)) / 2
  eig <- eigen(information, symmetric = TRUE)
  positive_definite <- all(eig$values > 1e-7)
  condition_number <- max(abs(eig$values)) / max(min(abs(eig$values)), 1e-12)
  covariance <- regularized_inverse_hessian(information, 1e-5)
  rownames(covariance) <- colnames(covariance) <- names(prior$mean)

  list(
    counts = counts,
    mesh = mesh,
    prior = prior,
    process_model = process_model,
    map = setNames(fit$par, names(prior$mean)),
    log_posterior_at_map = -fit$value,
    convergence = fit$convergence,
    covariance = covariance,
    information_eigenvalues = eig$values,
    positive_definite = positive_definite,
    condition_number = condition_number
  )
}

sample_laplace_draws <- function(
  fit,
  draws = 800L,
  seed = 1L,
  maximum_attempt_multiplier = 20L
) {
  set.seed(seed)
  eig <- eigen(fit$covariance, symmetric = TRUE)
  transformation <- eig$vectors %*%
    diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
  accepted <- matrix(
    NA_real_,
    nrow = 0,
    ncol = length(fit$map),
    dimnames = list(NULL, names(fit$map))
  )

  attempts <- 0L
  while (
    nrow(accepted) < draws &&
      attempts < maximum_attempt_multiplier * draws
  ) {
    batch_size <- min(2L * (draws - nrow(accepted)), draws)
    candidate <- matrix(
      rnorm(batch_size * length(fit$map)),
      nrow = batch_size
    ) %*% t(transformation)
    candidate <- sweep(candidate, 2, fit$map, "+")
    colnames(candidate) <- names(fit$map)

    valid <- apply(candidate, 1, function(parameter) {
      is.finite(inverse_ipm_log_posterior(
        parameter,
        fit$counts,
        fit$mesh,
        fit$prior,
        fit$process_model
      ))
    })
    accepted <- rbind(accepted, candidate[valid, , drop = FALSE])
    attempts <- attempts + batch_size
  }

  if (nrow(accepted) < draws) {
    warning("Fewer finite Laplace draws were available than requested.")
  }
  accepted[seq_len(min(draws, nrow(accepted))), , drop = FALSE]
}

summarize_draws_against_truth <- function(draws, truth) {
  aligned_truth <- truth[colnames(draws)]
  finite_truth <- is.finite(aligned_truth)
  quantiles <- t(apply(draws, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  posterior_sd <- apply(draws, 2, sd)
  coverage <- rep(NA, ncol(draws))
  coverage[finite_truth] <- quantiles[finite_truth, 1L] <=
    aligned_truth[finite_truth] &
    quantiles[finite_truth, 3L] >= aligned_truth[finite_truth]

  data.frame(
    quantity = colnames(draws),
    truth = aligned_truth,
    posterior_mean = colMeans(draws),
    posterior_sd = posterior_sd,
    q025 = quantiles[, 1L],
    median = quantiles[, 2L],
    q975 = quantiles[, 3L],
    covers_truth = coverage,
    standardized_bias = ifelse(
      finite_truth,
      (colMeans(draws) - aligned_truth) / posterior_sd,
      NA_real_
    ),
    row.names = NULL
  )
}

monte_carlo_fit_summary <- function(
  fit,
  draws = 800L,
  seed = 1L,
  evaluation_sizes = c(50, 80, 110),
  truth = inverse_ipm_truth()
) {
  posterior_draws <- sample_laplace_draws(fit, draws = draws, seed = seed)
  summary <- summarize_inverse_ipm_draws(
    posterior_draws,
    fit$mesh,
    evaluation_sizes,
    truth
  )

  list(
    parameter = summary$parameter,
    derived = summary$derived,
    diagnostics = data.frame(
      process_model = fit$process_model,
      convergence = fit$convergence,
      positive_definite = fit$positive_definite,
      condition_number = fit$condition_number,
      laplace_draws = nrow(posterior_draws),
      parameter_coverage = mean(
        summary$parameter$covers_truth,
        na.rm = TRUE
      ),
      derived_coverage = mean(summary$derived$covers_truth, na.rm = TRUE),
      parameter_mean_abs_z = mean(
        abs(summary$parameter$standardized_bias),
        na.rm = TRUE
      ),
      derived_mean_abs_z = mean(
        abs(summary$derived$standardized_bias),
        na.rm = TRUE
      )
    )
  )
}

summarize_inverse_ipm_draws <- function(
  posterior_draws,
  mesh,
  evaluation_sizes = c(50, 80, 110),
  truth = inverse_ipm_truth()
) {
  parameter_summary <- summarize_draws_against_truth(
    posterior_draws,
    truth
  )

  derived_draws <- t(apply(posterior_draws, 1, function(parameter) {
    derived_inverse_ipm_quantities(
      setNames(parameter, colnames(posterior_draws)),
      mesh,
      evaluation_sizes
    )
  }))
  derived_truth <- derived_inverse_ipm_quantities(
    truth,
    mesh,
    evaluation_sizes
  )
  derived_summary <- summarize_draws_against_truth(
    derived_draws,
    derived_truth
  )

  list(
    parameter = parameter_summary,
    derived = derived_summary
  )
}

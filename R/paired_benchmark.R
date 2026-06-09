# Paired profile, capture-recapture, integrated, and complete-history benchmark.

inverse_ipm_vital_rate_values <- function(parameter, size) {
  z <- (size - 80) / 20
  list(
    survival = inv_logit(
      parameter["survival_at_80"] +
        parameter["survival_slope_20"] * z
    ),
    growth_mean = pmax(
      size,
      size + parameter["growth_increment_80"] +
        parameter["growth_slope_20"] * z
    ),
    growth_sd = exp(parameter["log_growth_sd"]),
    recruitment = exp(
      parameter["log_recruitment_at_80"] +
        parameter["recruitment_slope_20"] * z
    ),
    recruit_mean = parameter["recruit_mean_80"] +
      parameter["recruit_mean_slope_20"] * z,
    recruit_sd = exp(parameter["log_recruit_sd"])
  )
}

prepare_capture_recapture_events <- function(
  captures,
  final_year = max(captures$year),
  recapture_detection = NULL
) {
  stopifnot(
    all(c("year", "id", "observed_size", "capture_probability") %in%
      names(captures))
  )
  if (is.null(recapture_detection)) {
    probabilities <- unique(captures$capture_probability)
    if (length(probabilities) != 1L) {
      stop("Supply recapture_detection when capture probability is not constant.")
    }
    recapture_detection <- probabilities
  }
  stopifnot(
    length(recapture_detection) == 1L,
    recapture_detection > 0,
    recapture_detection <= 1
  )

  origins <- captures[captures$year < final_year, c(
    "year", "id", "observed_size"
  )]
  names(origins)[names(origins) == "observed_size"] <- "from_size"
  next_key <- paste(captures$id, captures$year, sep = ":")
  origin_next_key <- paste(origins$id, origins$year + 1L, sep = ":")
  next_index <- match(origin_next_key, next_key)

  origins$recaptured <- !is.na(next_index)
  origins$to_size <- NA_real_
  origins$to_size[origins$recaptured] <-
    captures$observed_size[next_index[origins$recaptured]]
  origins$recapture_detection <- recapture_detection
  origins
}

capture_recapture_log_likelihood <- function(parameter, events) {
  rates <- inverse_ipm_vital_rate_values(parameter, events$from_size)
  recapture_probability <- rates$survival * events$recapture_detection
  if (
    any(!is.finite(recapture_probability)) ||
      any(recapture_probability <= 0) ||
      any(recapture_probability >= 1) ||
      !is.finite(rates$growth_sd) ||
      rates$growth_sd <= 0
  ) {
    return(-Inf)
  }

  log_likelihood <- sum(dbinom(
    events$recaptured,
    size = 1L,
    prob = recapture_probability,
    log = TRUE
  ))
  if (any(events$recaptured)) {
    log_likelihood <- log_likelihood + sum(dnorm(
      events$to_size[events$recaptured],
      mean = rates$growth_mean[events$recaptured],
      sd = rates$growth_sd,
      log = TRUE
    ))
  }
  log_likelihood
}

prepare_complete_history <- function(simulation) {
  transitions <- merge(
    simulation$transitions,
    simulation$census[, c("year", "id", "sex")],
    by = c("year", "id"),
    all.x = TRUE,
    sort = FALSE
  )
  transition_key <- paste(transitions$year, transitions$id, sep = ":")

  if (nrow(simulation$recruits) > 0L) {
    recruit_parent_key <- paste(
      simulation$recruits$year - 1L,
      simulation$recruits$parent_id,
      sep = ":"
    )
    offspring_table <- table(recruit_parent_key)
    transitions$offspring_count <- as.integer(
      offspring_table[transition_key]
    )
    transitions$offspring_count[is.na(transitions$offspring_count)] <- 0L
  } else {
    transitions$offspring_count <- 0L
  }

  list(
    transitions = transitions,
    recruits = simulation$recruits
  )
}

complete_history_log_likelihood <- function(parameter, histories) {
  if (!is.list(histories) || !is.null(histories$transitions)) {
    histories <- list(histories)
  }

  log_likelihood <- 0
  for (history in histories) {
    transitions <- history$transitions
    rates <- inverse_ipm_vital_rate_values(parameter, transitions$from_size)
    if (
      any(!is.finite(rates$survival)) ||
        any(rates$survival <= 0) ||
        any(rates$survival >= 1) ||
        !is.finite(rates$growth_sd) ||
        rates$growth_sd <= 0 ||
        any(!is.finite(rates$recruitment))
    ) {
      return(-Inf)
    }

    log_likelihood <- log_likelihood + sum(dbinom(
      transitions$survived,
      size = 1L,
      prob = rates$survival,
      log = TRUE
    ))
    if (any(transitions$survived)) {
      log_likelihood <- log_likelihood + sum(dnorm(
        transitions$to_size[transitions$survived],
        mean = rates$growth_mean[transitions$survived],
        sd = rates$growth_sd,
        log = TRUE
      ))
    }

    offspring_mean <- ifelse(
      transitions$sex == "F",
      2 * rates$recruitment,
      0
    )
    log_likelihood <- log_likelihood + sum(dpois(
      transitions$offspring_count,
      lambda = offspring_mean,
      log = TRUE
    ))

    recruits <- history$recruits
    if (nrow(recruits) > 0L) {
      recruit_rates <- inverse_ipm_vital_rate_values(
        parameter,
        recruits$parent_size
      )
      if (
        !is.finite(recruit_rates$recruit_sd) ||
          recruit_rates$recruit_sd <= 0
      ) {
        return(-Inf)
      }
      log_likelihood <- log_likelihood + sum(dnorm(
        recruits$size,
        mean = recruit_rates$recruit_mean,
        sd = recruit_rates$recruit_sd,
        log = TRUE
      ))
    }
  }
  log_likelihood
}

fit_laplace_custom <- function(
  log_posterior,
  prior,
  analysis,
  maxit = 1000L
) {
  starts <- c(list(prior$mean), prior$alternative_starts)
  fits <- lapply(starts, function(start) {
    tryCatch(
      optim(
        par = start,
        fn = function(parameter) {
          -log_posterior(setNames(parameter, names(prior$mean)))
        },
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
    stop("No finite custom Laplace posterior mode was found.")
  }
  fits <- fits[valid]
  fit <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]

  information <- (fit$hessian + t(fit$hessian)) / 2
  eig <- eigen(information, symmetric = TRUE)
  covariance <- regularized_inverse_hessian(information, 1e-5)
  rownames(covariance) <- colnames(covariance) <- names(prior$mean)

  list(
    analysis = analysis,
    prior = prior,
    map = setNames(fit$par, names(prior$mean)),
    log_posterior_at_map = -fit$value,
    convergence = fit$convergence,
    covariance = covariance,
    information_eigenvalues = eig$values,
    positive_definite = all(eig$values > 1e-7),
    condition_number = max(abs(eig$values)) /
      max(min(abs(eig$values)), 1e-12)
  )
}

sample_custom_laplace_draws <- function(
  fit,
  log_posterior,
  valid_parameter = function(parameter) TRUE,
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
    nrow = 0L,
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
    candidate <- sweep(candidate, 2L, fit$map, "+")
    colnames(candidate) <- names(fit$map)
    valid <- apply(candidate, 1L, function(parameter) {
      is.finite(log_posterior(parameter)) && isTRUE(valid_parameter(parameter))
    })
    accepted <- rbind(accepted, candidate[valid, , drop = FALSE])
    attempts <- attempts + batch_size
  }

  if (nrow(accepted) < draws) {
    warning("Fewer finite custom Laplace draws were available than requested.")
  }
  accepted[seq_len(min(draws, nrow(accepted))), , drop = FALSE]
}

parameter_group <- function(quantity) {
  if (grepl("process_precision", quantity)) {
    return("process")
  }
  if (grepl("^survival", quantity)) {
    return("survival")
  }
  if (grepl("^growth|^log_growth", quantity)) {
    return("growth")
  }
  if (grepl("^recruitment|^log_recruitment", quantity)) {
    return("recruitment")
  }
  if (grepl("^recruit_mean|^recruit_sd|^log_recruit_sd", quantity)) {
    return("recruit_size")
  }
  if (quantity == "lambda") {
    return("lambda")
  }
  "other"
}

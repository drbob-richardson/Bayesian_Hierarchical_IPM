# Flexible (spline-basis) inverse IPM, to test whether the recruitment
# confounding is an artifact of the parametric vital-rate forms. Survival,
# growth increment, log-recruitment, and recruit mean are each modeled as
# natural-cubic-spline functions of size; growth/recruit SD and the process
# precision are scalars. Reuses normal_bin_density() and normalize_columns()
# from R/numerical_ipm.R.

library(splines)

make_size_basis <- function(mesh, df = 3L) {
  B <- cbind(1, splines::ns(mesh, df = df))   # intercept + df natural-spline cols
  attr(B, "df1") <- ncol(B)
  B
}

spline_parameter_layout <- function(df1) {
  list(
    survival = 1:df1,
    growth = df1 + 1:df1,
    recruitment = 2 * df1 + 1:df1,
    recruit_mean = 3 * df1 + 1:df1,
    log_growth_sd = 4 * df1 + 1L,
    log_recruit_sd = 4 * df1 + 2L,
    log_process_precision = 4 * df1 + 3L
  )
}

build_spline_ipm <- function(parameter, mesh, B) {
  df1 <- ncol(B)
  L <- spline_parameter_layout(df1)
  survival <- plogis(as.vector(B %*% parameter[L$survival]))
  increment <- pmax(0, as.vector(B %*% parameter[L$growth]))
  growth_mean <- mesh + increment
  growth_sd <- exp(parameter[L$log_growth_sd])
  recruitment <- exp(as.vector(B %*% parameter[L$recruitment]))
  recruit_mean <- as.vector(B %*% parameter[L$recruit_mean])
  recruit_sd <- exp(parameter[L$log_recruit_sd])
  if (!all(is.finite(c(survival, growth_mean, growth_sd, recruitment, recruit_mean, recruit_sd))) ||
      growth_sd <= 0 || recruit_sd <= 0) {
    stop("Invalid spline vital rates.")
  }
  delta <- mean(diff(mesh))
  growth_density <- vapply(seq_along(mesh),
    function(j) normal_bin_density(mesh, growth_mean[j], growth_sd), numeric(length(mesh)))
  recruit_density <- vapply(seq_along(mesh),
    function(j) normal_bin_density(mesh, recruit_mean[j], recruit_sd), numeric(length(mesh)))
  growth_density <- normalize_columns(growth_density, delta)
  recruit_density <- normalize_columns(recruit_density, delta)
  projection <- growth_density * rep(survival, each = length(mesh))
  fertility <- recruit_density * rep(recruitment, each = length(mesh))
  list(mesh = mesh, delta = delta, projection = projection, fertility = fertility,
       kernel = projection + fertility, survival = survival, recruitment = recruitment,
       growth_increment = increment, recruit_mean = recruit_mean)
}

spline_expected_counts <- function(parameter, initial_counts, mesh, transitions, B) {
  ipm <- tryCatch(build_spline_ipm(parameter, mesh, B), error = function(e) NULL)
  if (is.null(ipm)) return(NULL)
  intensity <- initial_counts / ipm$delta
  out <- matrix(NA_real_, transitions + 1L, length(mesh)); out[1L, ] <- intensity
  for (s in seq_len(transitions)) out[s + 1L, ] <- as.vector(ipm$delta * ipm$kernel %*% out[s, ])
  out * ipm$delta
}

spline_log_posterior <- function(parameter, counts, mesh, B, prior_mean, prior_sd,
                                 process = c("poisson", "gamma_poisson")) {
  process <- match.arg(process)
  if (any(!is.finite(parameter))) return(-Inf)
  counts <- as_profile_trajectory_list(counts)
  log_prior <- sum(dnorm(parameter, prior_mean, prior_sd, log = TRUE))
  if (!is.finite(log_prior)) return(-Inf)
  log_lik <- 0
  L <- spline_parameter_layout(ncol(B))
  for (trajectory in counts) {
    expected <- spline_expected_counts(parameter, trajectory[1L, ], mesh,
                                       nrow(trajectory) - 1L, B)
    if (is.null(expected) || any(!is.finite(expected)) || any(expected < 0)) return(-Inf)
    expected <- pmax(expected[-1L, , drop = FALSE], 1e-12)
    observed <- trajectory[-1L, , drop = FALSE]
    if (process == "poisson") {
      log_lik <- log_lik + sum(dpois(observed, expected, log = TRUE))
    } else {
      size <- exp(parameter[L$log_process_precision])
      log_lik <- log_lik + sum(dnbinom(observed, mu = expected, size = size, log = TRUE))
    }
  }
  log_lik + log_prior
}

fit_spline_laplace <- function(counts, mesh, B, prior_mean, prior_sd,
                               process = "poisson", init = NULL, maxit = 2000L) {
  if (is.null(init)) init <- prior_mean
  neg <- function(p) -spline_log_posterior(p, counts, mesh, B, prior_mean, prior_sd, process)
  opt <- optim(init, neg, method = "BFGS",
               control = list(maxit = maxit, reltol = 1e-9), hessian = TRUE)
  information <- (opt$hessian + t(opt$hessian)) / 2
  eig <- eigen(information, symmetric = TRUE)
  covariance <- regularized_inverse_hessian(information, 1e-5)
  list(map = opt$par, covariance = covariance, convergence = opt$convergence,
       log_posterior_at_map = -opt$value, B = B,
       positive_definite = all(eig$values > 1e-7))
}

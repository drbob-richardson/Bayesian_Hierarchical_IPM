# Positive numerical-quadrature IPM scaffold.

normalize_columns <- function(matrix, delta) {
  column_mass <- colSums(matrix) * delta
  if (any(!is.finite(column_mass)) || any(column_mass <= 0)) {
    stop("Each conditional density must have positive finite mass on the mesh.")
  }
  matrix / rep(column_mass, each = nrow(matrix))
}

normal_bin_density <- function(mesh, mean, sd) {
  delta <- mean(diff(mesh))
  lower <- mesh - delta / 2
  upper <- mesh + delta / 2
  (
    pnorm(upper, mean = mean, sd = sd) -
      pnorm(lower, mean = mean, sd = sd)
  ) / delta
}

make_ipm_kernel <- function(
  mesh,
  survival,
  growth_mean,
  growth_sd,
  recruitment,
  recruit_mean,
  recruit_sd
) {
  stopifnot(length(mesh) >= 2L, all(diff(mesh) > 0))
  delta <- mean(diff(mesh))
  parent_size <- mesh

  survival_probability <- pmin(1, pmax(0, survival(parent_size)))
  realized_recruitment <- pmax(0, recruitment(parent_size))

  growth_density <- vapply(
    parent_size,
    function(size) {
      normal_bin_density(mesh, growth_mean(size), growth_sd(size))
    },
    numeric(length(mesh))
  )
  recruit_density <- vapply(
    parent_size,
    function(size) {
      normal_bin_density(mesh, recruit_mean(size), recruit_sd(size))
    },
    numeric(length(mesh))
  )

  growth_density <- normalize_columns(growth_density, delta)
  recruit_density <- normalize_columns(recruit_density, delta)

  projection <- growth_density * rep(survival_probability, each = length(mesh))
  fertility <- recruit_density * rep(realized_recruitment, each = length(mesh))

  list(
    mesh = mesh,
    delta = delta,
    projection = projection,
    fertility = fertility,
    kernel = projection + fertility,
    survival = survival_probability,
    recruitment = realized_recruitment
  )
}

project_ipm <- function(initial_intensity, ipm, transitions = 1L) {
  stopifnot(
    length(initial_intensity) == length(ipm$mesh),
    all(initial_intensity >= 0),
    transitions >= 0L
  )

  result <- matrix(
    NA_real_,
    nrow = transitions + 1L,
    ncol = length(initial_intensity)
  )
  result[1L, ] <- initial_intensity

  if (transitions > 0L) {
    for (step in seq_len(transitions)) {
      result[step + 1L, ] <- as.vector(
        ipm$delta * ipm$kernel %*% result[step, ]
      )
    }
  }

  result
}

poisson_profile_loglik <- function(observed_counts, expected_counts) {
  stopifnot(
    length(observed_counts) == length(expected_counts),
    all(observed_counts >= 0),
    all(expected_counts >= 0)
  )
  sum(dpois(observed_counts, lambda = expected_counts, log = TRUE))
}

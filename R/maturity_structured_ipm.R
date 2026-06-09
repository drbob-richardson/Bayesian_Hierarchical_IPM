# Recruitment structure using a known maturation curve.

maturity_probability <- function(size, midpoint = 75, scale = 5) {
  plogis((size - midpoint) / scale)
}

maturity_multiplier <- function(size, midpoint = 75, scale = 5) {
  maturity_probability(size, midpoint, scale) /
    maturity_probability(80, midpoint, scale)
}

build_maturity_centered_ipm <- function(
  parameter,
  mesh,
  maturity_midpoint = 75,
  maturity_scale = 5
) {
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
      maturity_multiplier(size, maturity_midpoint, maturity_scale) *
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

fish_vital_rates_from_maturity_ipm <- function(
  parameter,
  maturity_midpoint = 75,
  maturity_scale = 5
) {
  rates <- fish_vital_rates_from_inverse_ipm(parameter)
  z <- function(size) (size - 80) / 20
  rates$realized_recruitment <- function(
    size, age, sex, quality, environment, density
  ) {
    (sex == "F") * 2 *
      maturity_multiplier(size, maturity_midpoint, maturity_scale) *
      exp(
        parameter["log_recruitment_at_80"] +
          parameter["recruitment_slope_20"] * z(size)
      )
  }
  rates
}

maturity_inverse_ipm_quantities <- function(
  parameter,
  mesh,
  evaluation_sizes = c(50, 80, 110),
  maturity_midpoint = 75,
  maturity_scale = 5
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
  recruitment <- maturity_multiplier(
    evaluation_sizes,
    maturity_midpoint,
    maturity_scale
  ) * exp(
    parameter["log_recruitment_at_80"] +
      parameter["recruitment_slope_20"] * z(evaluation_sizes)
  )
  recruit_mean <- parameter["recruit_mean_80"] +
    parameter["recruit_mean_slope_20"] * z(evaluation_sizes)

  ipm <- build_maturity_centered_ipm(
    parameter,
    mesh,
    maturity_midpoint,
    maturity_scale
  )
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

summarize_draws_with_derived_target <- function(
  posterior_draws,
  mesh,
  draw_derived,
  truth_derived,
  truth = inverse_ipm_truth()
) {
  parameter_summary <- summarize_draws_against_truth(posterior_draws, truth)
  derived_draws <- t(apply(posterior_draws, 1L, function(parameter) {
    draw_derived(setNames(parameter, colnames(posterior_draws)), mesh)
  }))
  list(
    parameter = parameter_summary,
    derived = summarize_draws_against_truth(derived_draws, truth_derived)
  )
}

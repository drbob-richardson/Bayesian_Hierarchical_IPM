source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/identifiability.R")

mesh <- seq(25, 160, by = 2)
initial_intensity <- dnorm(mesh, 55, 10) * 500

parameter <- c(
  survival_intercept = -2.2,
  survival_slope = 0.032,
  growth_increment = 18,
  growth_slope = -0.12,
  log_growth_sd = log(4.5),
  recruitment_intercept = -7,
  recruitment_slope = 0.075,
  recruit_mean = 43,
  log_recruit_sd = log(3.5)
)

forward_profiles <- function(parameter, transitions = 3L) {
  ipm <- make_ipm_kernel(
    mesh = mesh,
    survival = function(size) {
      inv_logit(
        parameter["survival_intercept"] +
          parameter["survival_slope"] * size
      )
    },
    growth_mean = function(size) {
      pmax(
        size,
        size + parameter["growth_increment"] +
          parameter["growth_slope"] * size
      )
    },
    growth_sd = function(size) exp(parameter["log_growth_sd"]),
    recruitment = function(size) {
      0.5 * exp(
        parameter["recruitment_intercept"] +
          parameter["recruitment_slope"] * size
      )
    },
    recruit_mean = function(size) parameter["recruit_mean"] + 0.02 * size,
    recruit_sd = function(size) exp(parameter["log_recruit_sd"])
  )

  as.vector(project_ipm(initial_intensity, ipm, transitions = transitions)[-1, ])
}

for (transitions in c(1L, 2L, 4L, 8L)) {
  jacobian <- finite_difference_jacobian(
    parameter,
    function(candidate) forward_profiles(candidate, transitions)
  )
  diagnostic <- sensitivity_svd(jacobian)

  cat("\nTransitions:", transitions, "\n")
  cat("Numerical rank:", diagnostic$numerical_rank, "of", length(parameter), "\n")
  cat("Condition number:", format(diagnostic$condition_number, digits = 3), "\n")
  cat(
    "Smallest relative singular values:",
    paste(
      format(tail(diagnostic$relative_singular_values, 3), digits = 3),
      collapse = ", "
    ),
    "\n"
  )
  cat("Weakest direction:\n")
  weakest <- diagnostic$weak_directions[, 1L]
  print(sort(weakest, decreasing = TRUE))
}

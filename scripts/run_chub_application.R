# Real-data application: profile-only inverse IPM applied to two southern
# leatherside chub populations (Salina, Lost Creek). The populations CANNOT be
# pooled, so each is fitted separately with its own vital-rate parameters.
# This is a diagnostic application: we report what the short profile series
# identify (prior-to-posterior contraction) and what remains confounded
# (recruitment level/slope ridge, information conditioning), rather than
# claiming recovery.

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

process_model <- Sys.getenv("CHUB_PROCESS", "gamma_poisson")
mesh <- seq(25, 160, by = 5)
results_directory <- "results/chub_application"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

chub <- read.csv("data/chub.csv", stringsAsFactors = FALSE)
chub$size <- chub$Size
prior <- default_inverse_ipm_priors("weak", process_model)
evaluation_sizes <- c(50, 80, 110)

fit_population <- function(location) {
  census <- chub[chub$Location == location, c("size", "year")]
  counts <- make_profile_count_matrix(census, mesh)
  transitions <- nrow(counts) - 1L

  fit <- fit_laplace_inverse_ipm(counts, mesh, prior, process_model)
  draws <- sample_laplace_draws(fit, draws = 2000L, seed = 20330001)

  # Parameter-level posterior summary and prior-to-posterior contraction.
  param_summary <- data.frame(
    parameter = colnames(draws),
    map = fit$map[colnames(draws)],
    posterior_mean = colMeans(draws),
    posterior_sd = apply(draws, 2, sd),
    q025 = apply(draws, 2, quantile, 0.025),
    q500 = apply(draws, 2, quantile, 0.5),
    q975 = apply(draws, 2, quantile, 0.975),
    prior_sd = prior$sd[colnames(draws)],
    row.names = NULL
  )
  param_summary$contraction <- 1 - param_summary$posterior_sd / param_summary$prior_sd

  # Derived vital-rate curves and lambda with posterior intervals.
  derived_draws <- t(apply(draws, 1, function(p) {
    derived_inverse_ipm_quantities(setNames(p, colnames(draws)), mesh, evaluation_sizes)
  }))
  derived_summary <- data.frame(
    quantity = colnames(derived_draws),
    posterior_mean = colMeans(derived_draws),
    posterior_sd = apply(derived_draws, 2, sd),
    q025 = apply(derived_draws, 2, quantile, 0.025),
    q500 = apply(derived_draws, 2, quantile, 0.5),
    q975 = apply(derived_draws, 2, quantile, 0.975),
    row.names = NULL
  )

  # Confounding diagnostics.
  correlation <- cov2cor(fit$covariance)
  recruit_ridge <- correlation["log_recruitment_at_80", "recruitment_slope_20"]
  information <- solve(fit$covariance)
  info_eigen <- eigen((information + t(information)) / 2, symmetric = TRUE)
  flat_direction <- info_eigen$vectors[, which.min(info_eigen$values)]
  names(flat_direction) <- rownames(fit$covariance)
  condition_number <- max(info_eigen$values) / max(min(info_eigen$values), 1e-12)

  list(
    location = location, transitions = transitions, n_by_year = table(census$year),
    counts = counts, fit = fit, draws = draws,
    param_summary = param_summary, derived_summary = derived_summary,
    recruit_level_slope_correlation = recruit_ridge,
    flat_direction = flat_direction, condition_number = condition_number,
    positive_definite = fit$positive_definite
  )
}

results <- list()
for (location in unique(chub$Location)) {
  cat("==== Fitting", location, "(", process_model, ") ====\n")
  res <- fit_population(location)
  results[[location]] <- res
  safe <- gsub("[^A-Za-z0-9]+", "_", location)
  write.csv(res$param_summary,
    file.path(results_directory, sprintf("%s_parameters.csv", safe)), row.names = FALSE)
  write.csv(res$derived_summary,
    file.path(results_directory, sprintf("%s_derived.csv", safe)), row.names = FALSE)

  cat("Transitions:", res$transitions, "| PD Hessian:", res$positive_definite,
      "| info condition number:", format(res$condition_number, digits = 3), "\n")
  cat("Posterior corr(recruitment level, slope):",
      round(res$recruit_level_slope_correlation, 3), "\n")
  cat("Most- and least-identified parameters by contraction:\n")
  ps <- res$param_summary[order(-res$param_summary$contraction), ]
  print(head(ps[, c("parameter", "posterior_mean", "contraction")], 3), row.names = FALSE, digits = 3)
  print(tail(ps[, c("parameter", "posterior_mean", "contraction")], 3), row.names = FALSE, digits = 3)
  cat("Derived survival / recruitment / lambda (median [95%]):\n")
  show <- res$derived_summary[res$derived_summary$quantity %in%
    c("survival_80", "growth_increment_80", "recruitment_80", "recruit_mean_80", "lambda"), ]
  for (i in seq_len(nrow(show))) {
    cat(sprintf("  %-20s %.3f [%.3f, %.3f]\n", show$quantity[i],
                show$q500[i], show$q025[i], show$q975[i]))
  }
  cat("\n")
}

saveRDS(
  list(results = results, prior = prior, mesh = mesh, process_model = process_model),
  file.path(results_directory, "chub_fits.rds"), compress = "xz"
)
cat("Saved chub application to", results_directory, "\n")

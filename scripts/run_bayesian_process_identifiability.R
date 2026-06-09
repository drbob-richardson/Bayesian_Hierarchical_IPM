source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")

simulation <- readRDS("results/bayesian_identifiability_simulation.rds")
mesh <- seq(25, 160, by = 5)
transitions_to_fit <- c(2L, 4L, 8L, 12L)

fits <- list()
summaries <- list()

for (index in seq_along(transitions_to_fit)) {
  transitions <- transitions_to_fit[index]
  cat("\nFitting latent Gamma-Poisson model with", transitions, "transitions...\n")
  counts <- make_profile_count_matrix(
    simulation$census,
    mesh,
    transitions = transitions
  )
  prior <- default_inverse_ipm_priors("weak", process_model = "gamma_poisson")
  fit <- fit_bayesian_inverse_ipm(
    counts = counts,
    mesh = mesh,
    prior = prior,
    process_model = "gamma_poisson",
    chains = 4,
    iterations = 7000,
    warmup = 3500,
    seed = 20260707 + index * 100
  )
  summary <- posterior_summary(fit)
  summary$transitions <- transitions
  summary$prior_regime <- "weak"
  summary$process_model <- "gamma_poisson"
  fits[[paste0("gamma_poisson_", transitions)]] <- fit
  summaries[[index]] <- summary

  identifiable <- !is.na(summary$covers_truth)
  cat(
    "Acceptance:", paste(round(fit$acceptance, 2), collapse = ", "),
    "\nMax Rhat:", round(max(fit$rhat), 3),
    "\nMin ESS:", round(min(fit$ess)),
    "\nCoverage:", sum(summary$covers_truth[identifiable]), "of",
    sum(identifiable),
    "\nMedian contraction:",
    round(median(summary$contraction[identifiable]), 2),
    "\nProcess precision posterior:",
    paste(round(
      summary[summary$parameter == "log_process_precision", c(
        "q025", "median", "q975"
      )],
      2
    ), collapse = ", "),
    "\n"
  )
}

summary_table <- do.call(rbind, summaries)
write.csv(
  summary_table,
  "results/bayesian_process_identifiability_summary.csv",
  row.names = FALSE
)
saveRDS(
  fits,
  "results/bayesian_process_identifiability_fits.rds",
  compress = "xz"
)

derived <- do.call(rbind, lapply(seq_along(fits), function(index) {
  summary <- posterior_derived_summary(fits[[index]])
  summary$transitions <- transitions_to_fit[index]
  summary$prior_regime <- "weak"
  summary$process_model <- "gamma_poisson"
  summary
}))
write.csv(
  derived,
  "results/bayesian_process_identifiability_derived.csv",
  row.names = FALSE
)

cat("\nGamma-Poisson parameter coverage:\n")
print(aggregate(
  covers_truth ~ transitions,
  summary_table,
  mean,
  na.rm = TRUE
))
cat("\nGamma-Poisson derived-quantity coverage:\n")
print(aggregate(covers_truth ~ transitions, derived, mean))

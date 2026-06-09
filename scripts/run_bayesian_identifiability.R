source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")

set.seed(20260607)

mesh <- seq(25, 160, by = 5)
transitions_to_fit <- c(2L, 4L, 8L, 12L)
prior_regimes <- c("weak", "informative", "informative_wrong")

initial <- simulate_initial_population(
  n = 2500,
  min_size = 35,
  max_size = 145,
  latent_quality_sd = 0
)
simulation <- simulate_population(
  years = max(transitions_to_fit) + 1L,
  initial_population = initial,
  vital_rates = default_fish_vital_rates(),
  environment = rep(0, max(transitions_to_fit))
)

dir.create("results", showWarnings = FALSE)
fits <- list()
summaries <- list()
run_index <- 1L

for (prior_regime in prior_regimes) {
  for (transitions in transitions_to_fit) {
    cat(
      "\nFitting", transitions, "transitions with", prior_regime, "priors...\n"
    )
    counts <- make_profile_count_matrix(
      simulation$census,
      mesh,
      transitions = transitions
    )
    prior <- default_inverse_ipm_priors(prior_regime)
    fit <- fit_bayesian_inverse_ipm(
      counts = counts,
      mesh = mesh,
      prior = prior,
      chains = 4,
      iterations = 5000,
      warmup = 2500,
      seed = 20260607 + run_index * 100
    )
    summary <- posterior_summary(fit)
    summary$transitions <- transitions
    summary$prior_regime <- prior_regime
    summaries[[run_index]] <- summary
    fits[[paste(prior_regime, transitions, sep = "_")]] <- fit

    cat(
      "Acceptance:",
      paste(round(fit$acceptance, 2), collapse = ", "),
      "\nMax Rhat:", round(max(fit$rhat), 3),
      "\nMin ESS:", round(min(fit$ess)),
      "\nCoverage:", sum(summary$covers_truth), "of", nrow(summary),
      "\nMedian contraction:", round(median(summary$contraction), 2),
      "\n"
    )
    run_index <- run_index + 1L
  }
}

summary_table <- do.call(rbind, summaries)
write.csv(
  summary_table,
  "results/bayesian_identifiability_summary.csv",
  row.names = FALSE
)
saveRDS(fits, "results/bayesian_identifiability_fits.rds", compress = "xz")
saveRDS(simulation, "results/bayesian_identifiability_simulation.rds")

aggregate_summary <- do.call(rbind, lapply(
  split(summary_table, list(summary_table$transitions, summary_table$prior_regime)),
  function(group) {
    data.frame(
      transitions = group$transitions[1L],
      prior_regime = group$prior_regime[1L],
      coverage = mean(group$covers_truth),
      mean_contraction = mean(group$contraction),
      median_contraction = median(group$contraction),
      mean_abs_standardized_bias = mean(abs(group$standardized_bias)),
      max_abs_standardized_bias = max(abs(group$standardized_bias)),
      max_rhat = max(group$rhat),
      min_ess = min(group$ess)
    )
  }
))
write.csv(
  aggregate_summary,
  "results/bayesian_identifiability_aggregate.csv",
  row.names = FALSE
)

cat("\nAggregate posterior-identifiability summary:\n")
print(aggregate_summary)

cat("\nStrongest posterior correlations by fit:\n")
for (fit_name in names(fits)) {
  draws <- do.call(rbind, fits[[fit_name]]$chains)
  correlation <- cor(draws)
  correlation[lower.tri(correlation, diag = TRUE)] <- NA
  index <- which.max(abs(correlation))
  location <- arrayInd(index, dim(correlation))
  cat(
    fit_name, ":",
    rownames(correlation)[location[1L]], "vs",
    colnames(correlation)[location[2L]], "=",
    round(correlation[index], 3), "\n"
  )
}

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/correlated_process_ipm.R")

set.seed(20260807)

mesh <- seq(25, 160, by = 5)
designs <- data.frame(
  design = c("one_by_twelve", "three_by_four", "six_by_two"),
  populations = c(1L, 3L, 6L),
  transitions = c(12L, 4L, 2L),
  stringsAsFactors = FALSE
)

simulate_design <- function(populations, transitions, initial_n = 2500) {
  lapply(seq_len(populations), function(population_index) {
    initial <- simulate_initial_population(
      n = initial_n,
      min_size = 35,
      max_size = 145,
      latent_quality_sd = 0
    )
    simulation <- simulate_population(
      years = transitions + 1L,
      initial_population = initial,
      vital_rates = default_fish_vital_rates(),
      environment = rep(0, transitions)
    )
    make_profile_count_matrix(simulation$census, mesh)
  })
}

counts_by_design <- lapply(seq_len(nrow(designs)), function(index) {
  simulate_design(
    populations = designs$populations[index],
    transitions = designs$transitions[index]
  )
})
names(counts_by_design) <- designs$design

dir.create("results", showWarnings = FALSE)
saveRDS(counts_by_design, "results/replicate_design_counts.rds")

fits <- list()
summaries <- list()
derived_summaries <- list()
run_index <- 1L

for (design_index in seq_len(nrow(designs))) {
  design <- designs$design[design_index]
  counts <- counts_by_design[[design]]
  cat(
    "\nDesign:", design, "with", designs$populations[design_index],
    "population(s) and", designs$transitions[design_index],
    "transition(s) each.\n"
  )

  for (process_model in c("poisson", "gamma_poisson")) {
    cat("Fitting", process_model, "model...\n")
    prior <- default_inverse_ipm_priors(
      "weak",
      process_model = process_model
    )
    fit <- fit_bayesian_inverse_ipm(
      counts = counts,
      mesh = mesh,
      prior = prior,
      process_model = process_model,
      chains = 4,
      iterations = 6000,
      warmup = 3000,
      seed = 20260807 + run_index * 100
    )
    summary <- posterior_summary(fit)
    derived <- posterior_derived_summary(fit)

    summary$design <- design
    summary$process_model <- process_model
    derived$design <- design
    derived$process_model <- process_model
    fits[[paste(design, process_model, sep = "__")]] <- fit
    summaries[[run_index]] <- summary
    derived_summaries[[run_index]] <- derived

    identifiable <- !is.na(summary$covers_truth)
    cat(
      "Coverage:", mean(summary$covers_truth[identifiable]),
      "Derived coverage:", mean(derived$covers_truth),
      "Max Rhat:", round(max(fit$rhat), 3),
      "Min ESS:", round(min(fit$ess)), "\n"
    )
    run_index <- run_index + 1L
  }

  cat("Fitting correlated log-Gaussian process model...\n")
  fit <- fit_correlated_process_ipm(
    counts = counts,
    mesh = mesh,
    prior_regime = "weak",
    basis_rank = 4,
    chains = 4,
    iterations = 7000,
    warmup = 3500,
    seed = 20260807 + run_index * 100
  )
  summary <- posterior_correlated_process_summary(fit)
  derived <- posterior_derived_summary(fit)
  summary$design <- design
  summary$process_model <- "correlated_log_gaussian"
  derived$design <- design
  derived$process_model <- "correlated_log_gaussian"
  fits[[paste(design, "correlated_log_gaussian", sep = "__")]] <- fit
  summaries[[run_index]] <- summary
  derived_summaries[[run_index]] <- derived

  identifiable <- !is.na(summary$covers_truth)
  cat(
    "Coverage:", mean(summary$covers_truth[identifiable]),
    "Derived coverage:", mean(derived$covers_truth),
    "Max Rhat:", round(max(fit$rhat), 3),
    "Min ESS:", round(min(fit$ess)),
    "Process SD median:",
    round(summary$median[summary$parameter == "log_process_sd"], 3),
    "\n"
  )
  run_index <- run_index + 1L
}

summary_table <- do.call(rbind, summaries)
derived_table <- do.call(rbind, derived_summaries)
write.csv(
  summary_table,
  "results/replicate_identifiability_summary.csv",
  row.names = FALSE
)
write.csv(
  derived_table,
  "results/replicate_identifiability_derived.csv",
  row.names = FALSE
)
saveRDS(fits, "results/replicate_identifiability_fits.rds", compress = "xz")

aggregate_results <- do.call(rbind, lapply(
  split(summary_table, list(summary_table$design, summary_table$process_model)),
  function(group) {
    identified <- !is.na(group$covers_truth)
    data.frame(
      design = group$design[1L],
      process_model = group$process_model[1L],
      coverage = mean(group$covers_truth[identified]),
      mean_contraction = mean(group$contraction[identified]),
      mean_abs_standardized_bias = mean(
        abs(group$standardized_bias[identified])
      ),
      max_rhat = max(group$rhat),
      min_ess = min(group$ess)
    )
  }
))
derived_aggregate <- aggregate(
  covers_truth ~ design + process_model,
  derived_table,
  mean
)
names(derived_aggregate)[3L] <- "derived_coverage"
aggregate_results <- merge(
  aggregate_results,
  derived_aggregate,
  by = c("design", "process_model")
)
write.csv(
  aggregate_results,
  "results/replicate_identifiability_aggregate.csv",
  row.names = FALSE
)

cat("\nReplicate-design comparison:\n")
print(aggregate_results[order(
  aggregate_results$process_model,
  aggregate_results$design
), ])

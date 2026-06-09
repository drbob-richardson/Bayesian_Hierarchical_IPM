source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
scenario_grid <- read.csv("results/monte_carlo/scenario_grid.csv")

selected <- scenario_grid[
  scenario_grid$initial_profile == "baseline" &
    scenario_grid$replicate == 1 &
    (
      (
        scenario_grid$initial_n == 500 &
          scenario_grid$profile_detection == 0.25
      ) |
        (
          scenario_grid$initial_n == 2500 &
            scenario_grid$profile_detection == 1
        )
    ),
]

simulate_initial_for_scenario <- function(n, initial_profile) {
  if (initial_profile == "baseline") {
    return(simulate_initial_population(n = n, latent_quality_sd = 0))
  }
  population <- simulate_initial_population(n = n, latent_quality_sd = 0)
  pulse <- sample.int(3, n, replace = TRUE, prob = c(0.88, 0.10, 0.02))
  population$size <- rnorm(
    n,
    mean = c(43, 78, 108)[pulse],
    sd = c(3.5, 9, 11)[pulse]
  )
  population$size <- pmin(145, pmax(35, population$size))
  population$age <- pmax(0L, round((population$size - 35) / 18))
  population
}

regenerate_counts <- function(scenario) {
  set.seed(20261000 + scenario$scenario_id)
  lapply(seq_len(scenario$populations), function(population_index) {
    initial <- simulate_initial_for_scenario(
      scenario$initial_n,
      scenario$initial_profile
    )
    simulation <- simulate_population(
      years = scenario$transitions + 1L,
      initial_population = initial,
      vital_rates = default_fish_vital_rates(),
      environment = rep(0, scenario$transitions)
    )
    profiles <- sample_size_profiles(
      simulation$census,
      detection = function(size, year) {
        rep(scenario$profile_detection, length(size))
      },
      measurement_sd = 0
    )
    make_profile_count_matrix(profiles, mesh)
  })
}

audit_scenario <- function(scenario) {
  output <- sprintf(
    "results/monte_carlo/mcmc_audit_scenario_%04d.rds",
    scenario$scenario_id
  )
  if (file.exists(output)) {
    return(output)
  }

  counts <- regenerate_counts(scenario)
  fit <- fit_bayesian_inverse_ipm(
    counts = counts,
    mesh = mesh,
    prior = default_inverse_ipm_priors("weak", "poisson"),
    process_model = "poisson",
    chains = 4,
    iterations = 5000,
    warmup = 2500,
    seed = 20263000 + scenario$scenario_id
  )
  mcmc_summary <- posterior_summary(fit)
  mcmc_derived <- posterior_derived_summary(fit)

  task <- readRDS(sprintf(
    "results/monte_carlo/tasks/scenario_%04d.rds",
    scenario$scenario_id
  ))
  laplace_summary <- task$models$poisson$summary$parameter
  laplace_derived <- task$models$poisson$summary$derived

  result <- list(
    scenario = scenario,
    fit = fit,
    mcmc_summary = mcmc_summary,
    mcmc_derived = mcmc_derived,
    laplace_summary = laplace_summary,
    laplace_derived = laplace_derived
  )
  saveRDS(result, output, compress = "xz")
  output
}

paths <- mclapply(
  split(selected, seq_len(nrow(selected))),
  audit_scenario,
  mc.cores = min(6L, nrow(selected)),
  mc.preschedule = FALSE
)

audits <- lapply(paths, readRDS)
comparison <- do.call(rbind, lapply(audits, function(audit) {
  scenario <- audit$scenario
  merge_table <- merge(
    audit$mcmc_summary[
      , c("parameter", "posterior_mean", "posterior_sd", "covers_truth")
    ],
    audit$laplace_summary[
      , c("quantity", "posterior_mean", "posterior_sd", "covers_truth")
    ],
    by.x = "parameter",
    by.y = "quantity",
    suffixes = c("_mcmc", "_laplace")
  )
  data.frame(
    scenario_id = scenario$scenario_id,
    design = scenario$design,
    initial_n = scenario$initial_n,
    profile_detection = scenario$profile_detection,
    parameter = merge_table$parameter,
    posterior_mean_mcmc = merge_table$posterior_mean_mcmc,
    posterior_mean_laplace = merge_table$posterior_mean_laplace,
    posterior_sd_mcmc = merge_table$posterior_sd_mcmc,
    posterior_sd_laplace = merge_table$posterior_sd_laplace,
    mean_difference_in_mcmc_sd = (
      merge_table$posterior_mean_laplace -
        merge_table$posterior_mean_mcmc
    ) / merge_table$posterior_sd_mcmc,
    sd_ratio_laplace_to_mcmc = merge_table$posterior_sd_laplace /
      merge_table$posterior_sd_mcmc,
    coverage_mcmc = merge_table$covers_truth_mcmc,
    coverage_laplace = merge_table$covers_truth_laplace,
    max_rhat = max(audit$fit$rhat),
    min_ess = min(audit$fit$ess)
  )
}))
write.csv(
  comparison,
  "results/monte_carlo/mcmc_audit_comparison.csv",
  row.names = FALSE
)

cat("MCMC audit scenarios:", length(audits), "\n")
cat("Maximum R-hat:", max(comparison$max_rhat), "\n")
cat("Minimum ESS:", min(comparison$min_ess), "\n")
cat(
  "Median absolute mean difference in MCMC SD units:",
  median(abs(comparison$mean_difference_in_mcmc_sd)),
  "\n"
)
cat(
  "Median Laplace/MCMC posterior SD ratio:",
  median(comparison$sd_ratio_laplace_to_mcmc),
  "\n"
)
cat(
  "Parameter coverage agreement:",
  mean(comparison$coverage_mcmc == comparison$coverage_laplace),
  "\n"
)

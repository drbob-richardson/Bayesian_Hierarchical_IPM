source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")

task_files <- list.files(
  "results/monte_carlo/tasks",
  pattern = "\\.rds$",
  full.names = TRUE
)
results <- lapply(task_files, readRDS)

diagnostics <- list()
parameters <- list()
derived <- list()
index <- 1L

for (result in results) {
  scenario <- result$scenario
  for (model in result$models) {
    base <- data.frame(
      scenario_id = scenario$scenario_id,
      design = scenario$design,
      populations = scenario$populations,
      transitions = scenario$transitions,
      initial_n = scenario$initial_n,
      profile_detection = scenario$profile_detection,
      initial_profile = scenario$initial_profile,
      replicate = scenario$replicate,
      process_model = model$model,
      failed = !is.null(model$error),
      error = ifelse(is.null(model$error), NA, model$error),
      stringsAsFactors = FALSE
    )
    if (is.null(model$error)) {
      diagnostics[[index]] <- cbind(base, model$summary$diagnostics)
      scientific_base <- cbind(
        base,
        positive_definite = model$fit$positive_definite,
        condition_number = model$fit$condition_number
      )
      parameters[[index]] <- cbind(scientific_base, model$summary$parameter)
      derived[[index]] <- cbind(scientific_base, model$summary$derived)
    } else {
      diagnostics[[index]] <- base
    }
    index <- index + 1L
  }
}

diagnostic_table <- do.call(rbind, diagnostics)
parameter_table <- do.call(rbind, parameters)
derived_table <- do.call(rbind, derived)
write.csv(
  diagnostic_table,
  "results/monte_carlo/diagnostics.csv",
  row.names = FALSE
)
write.csv(
  parameter_table,
  "results/monte_carlo/parameter_summary.csv",
  row.names = FALSE
)
write.csv(
  derived_table,
  "results/monte_carlo/derived_summary.csv",
  row.names = FALSE
)

group_columns <- c(
  "design", "initial_n", "profile_detection", "initial_profile",
  "process_model"
)
group_id <- interaction(diagnostic_table[group_columns], drop = TRUE)
aggregate_table <- do.call(rbind, lapply(
  split(diagnostic_table, group_id),
  function(group) {
    successful <- group[
      !group$failed & group$positive_definite,
      ,
      drop = FALSE
    ]
    data.frame(
      design = group$design[1L],
      initial_n = group$initial_n[1L],
      profile_detection = group$profile_detection[1L],
      initial_profile = group$initial_profile[1L],
      process_model = group$process_model[1L],
      attempted = nrow(group),
      successful = nrow(successful),
      positive_definite_rate = mean(successful$positive_definite),
      mean_parameter_coverage = mean(successful$parameter_coverage),
      mean_derived_coverage = mean(successful$derived_coverage),
      mean_parameter_abs_z = mean(successful$parameter_mean_abs_z),
      mean_derived_abs_z = mean(successful$derived_mean_abs_z),
      median_condition_number = median(successful$condition_number)
    )
  }
))
write.csv(
  aggregate_table,
  "results/monte_carlo/aggregate.csv",
  row.names = FALSE
)

cat("Monte Carlo fits:", nrow(diagnostic_table), "\n")
cat("Failures:", sum(diagnostic_table$failed), "\n")
cat("Positive-definite Hessians:", mean(
  diagnostic_table$positive_definite,
  na.rm = TRUE
), "\n\n")
print(aggregate(
  cbind(parameter_coverage, derived_coverage) ~ design + process_model,
  diagnostic_table[
    !diagnostic_table$failed & diagnostic_table$positive_definite,
  ],
  mean
))

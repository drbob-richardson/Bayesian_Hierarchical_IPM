source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
designs <- data.frame(
  design = c("one_by_twelve", "three_by_four", "six_by_two"),
  populations = c(1L, 3L, 6L),
  transitions = c(12L, 4L, 2L),
  stringsAsFactors = FALSE
)

scenario_grid <- merge(
  designs,
  expand.grid(
    initial_n = c(500L, 2500L),
    profile_detection = c(0.25, 1.00),
    initial_profile = c("baseline", "cohort_pulse"),
    replicate = seq_len(10L),
    stringsAsFactors = FALSE
  )
)
scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))

task_directory <- "results/monte_carlo/tasks"
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  scenario_grid,
  "results/monte_carlo/scenario_grid.csv",
  row.names = FALSE
)

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

run_scenario <- function(scenario) {
  output_path <- file.path(
    task_directory,
    sprintf("scenario_%04d.rds", scenario$scenario_id)
  )
  if (file.exists(output_path)) {
    return(output_path)
  }

  set.seed(20261000 + scenario$scenario_id)
  counts <- lapply(seq_len(scenario$populations), function(population_index) {
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

  model_results <- lapply(c("poisson", "gamma_poisson"), function(model) {
    tryCatch({
      prior <- default_inverse_ipm_priors("weak", model)
      fit <- fit_laplace_inverse_ipm(counts, mesh, prior, model)
      summary <- monte_carlo_fit_summary(
        fit,
        draws = 500L,
        seed = 20262000 + scenario$scenario_id +
          ifelse(model == "poisson", 0L, 10000L)
      )
      list(model = model, fit = fit, summary = summary, error = NULL)
    }, error = function(error) {
      list(
        model = model,
        fit = NULL,
        summary = NULL,
        error = conditionMessage(error)
      )
    })
  })
  names(model_results) <- c("poisson", "gamma_poisson")

  result <- list(
    scenario = scenario,
    annual_counts = lapply(counts, rowSums),
    models = model_results
  )
  saveRDS(result, output_path, compress = "xz")
  output_path
}

pending <- split(scenario_grid, seq_len(nrow(scenario_grid)))
available_cores <- detectCores(logical = FALSE)
if (!is.finite(available_cores)) {
  available_cores <- detectCores(logical = TRUE)
}
system_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN",
  intern = TRUE
)))
if (length(system_cores) == 1L && is.finite(system_cores)) {
  available_cores <- max(available_cores, system_cores, na.rm = TRUE)
}
if (!is.finite(available_cores)) {
  available_cores <- 2L
}
cores <- max(1L, min(6L, available_cores - 1L))
cat("Running", length(pending), "scenarios on", cores, "cores.\n")
paths <- mclapply(pending, run_scenario, mc.cores = cores, mc.preschedule = FALSE)
saveRDS(paths, "results/monte_carlo/task_paths.rds")
cat("Completed Monte Carlo task grid.\n")

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

library(parallel)

mesh <- seq(25, 160, by = 5)
replicates <- 30L
data_directory <- "results/auxiliary_confirmatory/data"
task_directory <- "results/auxiliary_confirmatory/tasks"
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)

make_dataset <- function(replicate) {
  path <- file.path(data_directory, sprintf("replicate_%03d.rds", replicate))
  if (file.exists(path)) {
    return(path)
  }
  set.seed(20268000 + replicate)
  counts <- lapply(seq_len(3L), function(population_index) {
    initial <- simulate_initial_population(n = 500, latent_quality_sd = 0)
    simulation <- simulate_population(
      years = 5L,
      initial_population = initial,
      vital_rates = default_fish_vital_rates(),
      environment = rep(0, 4L)
    )
    profiles <- sample_size_profiles(
      simulation$census,
      detection = function(size, year) rep(0.25, length(size)),
      measurement_sd = 0
    )
    make_profile_count_matrix(profiles, mesh)
  })
  saveRDS(counts, path, compress = "xz")
  path
}

regimes <- data.frame(
  regime = c(
    "weak",
    "correct_full_50",
    "correct_full_200",
    "correct_full_800",
    "biased_full_800",
    "biased_robust_800"
  ),
  sample_size = c(0L, 50L, 200L, 800L, 800L, 800L),
  accuracy = c("weak", "correct", "correct", "correct", "biased", "biased"),
  borrowing = c("none", "full", "full", "full", "full", "robust"),
  stringsAsFactors = FALSE
)
tasks <- merge(
  data.frame(replicate = seq_len(replicates)),
  regimes,
  all = TRUE
)
tasks$task_id <- seq_len(nrow(tasks))
write.csv(
  tasks,
  "results/auxiliary_confirmatory/task_grid.csv",
  row.names = FALSE
)

available_cores <- suppressWarnings(as.integer(system(
  "getconf _NPROCESSORS_ONLN",
  intern = TRUE
)))
if (!is.finite(available_cores)) {
  available_cores <- detectCores(logical = TRUE)
}
cores <- max(1L, min(6L, available_cores - 1L))

cat("Generating", replicates, "confirmatory datasets on", cores, "cores.\n")
data_paths <- mclapply(
  seq_len(replicates),
  make_dataset,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(data_paths, "results/auxiliary_confirmatory/data_paths.rds")

run_task <- function(task) {
  output_path <- file.path(
    task_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(output_path)) {
    return(output_path)
  }
  counts <- readRDS(file.path(
    data_directory,
    sprintf("replicate_%03d.rds", task$replicate)
  ))
  prior <- if (task$accuracy == "weak") {
    default_inverse_ipm_priors("weak", "gamma_poisson")
  } else {
    make_fecundity_study_prior(
      sample_size = task$sample_size,
      accuracy = task$accuracy,
      borrowing = task$borrowing,
      process_model = "gamma_poisson",
      robust_weight = 0.8,
      seed = 20269000 + task$replicate + 10L * task$sample_size +
        ifelse(task$accuracy == "biased", 100000L, 0L)
    )
  }

  result <- tryCatch({
    if (task$borrowing == "robust") {
      fit <- fit_robust_laplace_inverse_ipm(
        counts,
        mesh,
        prior,
        "gamma_poisson"
      )
      draws <- sample_robust_laplace_draws(
        fit,
        draws = 1000L,
        seed = 20270000 + task$task_id
      )
      summary <- summarize_inverse_ipm_draws(draws, mesh)
      robust_weight <- fit$posterior_external_weight
    } else {
      fit <- fit_laplace_inverse_ipm(counts, mesh, prior, "gamma_poisson")
      summary <- monte_carlo_fit_summary(
        fit,
        draws = 1000L,
        seed = 20270000 + task$task_id
      )
      robust_weight <- NA_real_
    }
    list(
      task = task,
      fit = fit,
      summary = summary,
      robust_weight = robust_weight,
      error = NULL
    )
  }, error = function(error) {
    list(
      task = task,
      fit = NULL,
      summary = NULL,
      robust_weight = NA_real_,
      error = conditionMessage(error)
    )
  })
  saveRDS(result, output_path, compress = "xz")
  output_path
}

cat("Running", nrow(tasks), "confirmatory fits on", cores, "cores.\n")
paths <- mclapply(
  split(tasks, seq_len(nrow(tasks))),
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
saveRDS(paths, "results/auxiliary_confirmatory/task_paths.rds")
cat("Completed confirmatory auxiliary-fecundity experiment.\n")

experiment <- Sys.getenv("PROFILE_ALT_EXPERIMENT", "main")
results_directory <- file.path("results/profile_process_alternatives", experiment)
derived <- read.csv(file.path(results_directory, "derived_results.csv"))

artifact_directory <- file.path("artifacts/profile_process_alternatives", experiment)
dir.create(artifact_directory, recursive = TRUE, showWarnings = FALSE)

scenario_order <- c(
  "standard", "dense_profiles", "longer_followup", "diverse_initial",
  "rich_design", "parent_excitation_sparse", "parent_excitation_dense"
)
scenario_labels <- c(
  "Standard", "Complete profiles", "Longer follow-up", "Diverse initial",
  "Rich design", "Parent excitation: sparse", "Parent excitation: dense"
)
method_order <- c(
  "recursive_poisson",
  "recursive_gamma_poisson",
  "one_step_poisson",
  "one_step_gamma_poisson",
  "one_step_shared_gamma_poisson",
  "one_step_finite_population",
  "one_step_finite_population_overdispersed",
  "one_step_gamma_poisson_shape_informed",
  "one_step_gamma_poisson_fixed_recruitment_slope"
)
method_labels <- c(
  "Recursive Poisson",
  "Recursive Gamma-Poisson",
  "One-step Poisson",
  "One-step Gamma-Poisson",
  "One-step shared-Gamma Cox",
  "Finite-population moments",
  "Finite-population + overdispersion",
  "One-step, informed fecundity shape",
  "One-step, fixed fecundity shape"
)
method_colors <- c(
  "#4C78A8", "#72B7B2", "#F58518", "#ECA82C", "#54A24B",
  "#B279A2", "#9D755D", "#BAB0AC", "#E45756"
)
scenario_labels <- setNames(scenario_labels, scenario_order)
method_labels <- setNames(method_labels, method_order)
method_colors <- setNames(method_colors, method_order)
scenario_order <- intersect(scenario_order, unique(derived$scenario))
method_order <- intersect(method_order, unique(derived$method))
scenario_labels <- unname(scenario_labels[scenario_order])
method_labels <- unname(method_labels[method_order])
method_colors <- unname(method_colors[method_order])
selected <- c(
  "survival_80",
  "growth_increment_80",
  "recruitment_80",
  "lambda"
)
selected_labels <- c(
  "Survival at 80 mm",
  "Growth increment at 80 mm",
  "Recruitment at 80 mm",
  "Population growth rate"
)

coverage <- aggregate(
  covers_truth ~ scenario + method + quantity,
  derived[derived$quantity %in% selected, ],
  mean,
  na.rm = TRUE
)

png(
  file.path(artifact_directory, "coverage_selected_quantities.png"),
  width = 2200,
  height = 1800,
  res = 180
)
par(mfrow = c(2, 2), mar = c(7, 5, 3, 1))
for (quantity_index in seq_along(selected)) {
  quantity <- selected[quantity_index]
  values <- matrix(
    NA_real_,
    nrow = length(method_order),
    ncol = length(scenario_order)
  )
  for (method_index in seq_along(method_order)) {
    for (scenario_index in seq_along(scenario_order)) {
      row <- coverage[
        coverage$method == method_order[method_index] &
          coverage$scenario == scenario_order[scenario_index] &
          coverage$quantity == quantity,
      ]
      if (nrow(row) == 1L) {
        values[method_index, scenario_index] <- row$covers_truth
      }
    }
  }
  barplot(
    values,
    beside = TRUE,
    col = method_colors,
    border = NA,
    names.arg = scenario_labels,
    las = 2,
    ylim = c(0, 1),
    ylab = "Empirical 95% interval coverage",
    main = selected_labels[quantity_index]
  )
  abline(h = 0.95, lty = 2, col = "gray35")
}
legend(
  "bottom",
  inset = -0.05,
  xpd = NA,
  horiz = TRUE,
  legend = method_labels,
  fill = method_colors,
  border = NA,
  bty = "n",
  cex = 0.8
)
dev.off()

derived$squared_error <- (derived$posterior_mean - derived$truth)^2
rmse <- aggregate(
  squared_error ~ scenario + method + quantity,
  derived[derived$quantity %in% selected, ],
  mean
)
rmse$rmse <- sqrt(rmse$squared_error)

png(
  file.path(artifact_directory, "rmse_selected_quantities.png"),
  width = 2200,
  height = 1800,
  res = 180
)
par(mfrow = c(2, 2), mar = c(7, 5, 3, 1))
for (quantity_index in seq_along(selected)) {
  quantity <- selected[quantity_index]
  values <- matrix(
    NA_real_,
    nrow = length(method_order),
    ncol = length(scenario_order)
  )
  for (method_index in seq_along(method_order)) {
    for (scenario_index in seq_along(scenario_order)) {
      row <- rmse[
        rmse$method == method_order[method_index] &
          rmse$scenario == scenario_order[scenario_index] &
          rmse$quantity == quantity,
      ]
      if (nrow(row) == 1L) {
        values[method_index, scenario_index] <- row$rmse
      }
    }
  }
  barplot(
    values,
    beside = TRUE,
    col = method_colors,
    border = NA,
    names.arg = scenario_labels,
    las = 2,
    ylab = "RMSE",
    main = selected_labels[quantity_index]
  )
}
legend(
  "bottom",
  inset = -0.05,
  xpd = NA,
  horiz = TRUE,
  legend = method_labels,
  fill = method_colors,
  border = NA,
  bty = "n",
  cex = 0.8
)
dev.off()

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")

fits <- readRDS("results/bayesian_identifiability_fits.rds")
parameter_summary <- read.csv("results/bayesian_identifiability_summary.csv")

derived <- do.call(rbind, lapply(names(fits), function(fit_name) {
  components <- strsplit(fit_name, "_", fixed = TRUE)[[1L]]
  transitions <- as.integer(tail(components, 1L))
  prior_regime <- paste(head(components, -1L), collapse = "_")
  summary <- posterior_derived_summary(fits[[fit_name]])
  summary$transitions <- transitions
  summary$prior_regime <- prior_regime
  summary
}))
write.csv(
  derived,
  "results/bayesian_identifiability_derived.csv",
  row.names = FALSE
)

derived_aggregate <- do.call(rbind, lapply(
  split(derived, list(derived$transitions, derived$prior_regime)),
  function(group) {
    data.frame(
      transitions = group$transitions[1L],
      prior_regime = group$prior_regime[1L],
      coverage = mean(group$covers_truth),
      mean_relative_error = mean(
        abs(group$posterior_mean - group$truth) / pmax(abs(group$truth), 0.1)
      )
    )
  }
))
write.csv(
  derived_aggregate,
  "results/bayesian_identifiability_derived_aggregate.csv",
  row.names = FALSE
)

cat("Parameter coverage by fit:\n")
print(aggregate(
  covers_truth ~ transitions + prior_regime,
  parameter_summary,
  mean
))
cat("\nDerived-quantity coverage and relative error by fit:\n")
print(derived_aggregate)

cat("\nTwelve-transition weak-prior derived quantities:\n")
print(derived[
  derived$transitions == 12 & derived$prior_regime == "weak",
  c(
    "quantity", "truth", "posterior_mean", "posterior_sd",
    "q025", "q975", "covers_truth"
  )
])

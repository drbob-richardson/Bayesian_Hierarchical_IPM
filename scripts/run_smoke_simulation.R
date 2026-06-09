source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")

set.seed(20260607)

initial <- simulate_initial_population(n = 500, latent_quality_sd = 0.25)
truth <- simulate_population(
  years = 6,
  initial_population = initial,
  environment = c(0.0, 0.5, -0.25, 0.25, 0.0)
)

profiles <- sample_size_profiles(
  truth$census,
  detection = function(size, year) {
    inv_logit(-0.6 + 0.012 * (size - 70))
  },
  measurement_sd = 1
)

mark_recapture <- sample_mark_recapture(
  truth$census,
  detection = function(size, year) rep(0.35, length(size)),
  measurement_sd = 1
)

profile_summary <- aggregate(id ~ year, profiles, length)
names(profile_summary)[2] <- "profile_n"
capture_summary <- aggregate(id ~ year, mark_recapture, length)
names(capture_summary)[2] <- "capture_n"

cat("Latent census counts by year:\n")
print(table(truth$census$year))
cat("\nObserved size-profile counts by year:\n")
print(profile_summary)
cat("\nObserved mark-recapture counts by year:\n")
print(capture_summary)

mesh <- seq(25, 160, length.out = 136)
ipm <- make_ipm_kernel(
  mesh = mesh,
  survival = function(size) inv_logit(-2.2 + 0.032 * size),
  growth_mean = function(size) pmax(size, size + 18 - 0.12 * size),
  growth_sd = function(size) 4.5,
  recruitment = function(size) 0.5 * exp(-7 + 0.075 * size),
  recruit_mean = function(size) 43 + 0.02 * size,
  recruit_sd = function(size) 3.5
)

initial_intensity <- hist(
  truth$census$size[truth$census$year == 0],
  breaks = c(mesh - ipm$delta / 2, tail(mesh, 1) + ipm$delta / 2),
  plot = FALSE
)$counts / ipm$delta

projected <- project_ipm(initial_intensity, ipm, transitions = 5)
cat("\nNumerical IPM projected abundance by year:\n")
print(round(rowSums(projected) * ipm$delta, 1))

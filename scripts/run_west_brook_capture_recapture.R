# West Brook brook trout: mark-recapture and integrated analyses, to benchmark
# the profile-only inverse IPM (scripts/run_west_brook_application.R).
#
# Detection is UNKNOWN for real data (unlike the simulated paired benchmark where
# capture probability was known). We therefore (1) estimate detection p with a
# constant-survival / constant-detection CJS over the full 20-occasion histories
# -- which separates survival from detection -- then (2) plug p-hat into a
# size-dependent survival+growth recapture fit, and (3) an integrated fit that
# adds the recapture likelihood to the profile likelihood. Growth from
# recaptures is detection-free and is the cleanest cross-check of the profile
# model.

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")
source("R/paired_benchmark.R")

results_directory <- "results/west_brook_application"
year_min <- as.integer(Sys.getenv("WB_YEAR_MIN", "2000"))
year_max <- as.integer(Sys.getenv("WB_YEAR_MAX", "2019"))

profile <- readRDS(file.path(results_directory, "profile_fit.rds"))
# Match the profile's observation window, mesh, and (recentred) prior.
mesh <- profile$mesh
evaluation_sizes <- c(80, 110, 140)

# Whole-network, tagged fish only; one annual record per fish = its earliest
# capture that year (by day-of-year), detected in ANY reach or season. This
# defines survival as "alive somewhere in the West Brook network next year",
# removing the main-stem emigration artifact.
raw <- read.csv("data/west_brook/cdWB_electro_DR.csv", stringsAsFactors = FALSE)
tg <- raw$species == "bkt" & raw$year >= year_min & raw$year <= year_max &
  !is.na(raw$observedLength) & !(is.na(raw$tag) | raw$tag == "")
b <- raw[tg, c("tag", "year", "yday", "observedLength")]
b$size <- b$observedLength
b <- b[order(b$tag, b$year, b$yday), ]
per_year <- b[!duplicated(paste(b$tag, b$year)), c("tag", "year", "size")]
cat(sprintf("Whole-network tagged brook trout: %d annual records, %d tags, %d years.\n",
            nrow(per_year), length(unique(per_year$tag)), length(unique(per_year$year))))
occasions <- sort(unique(per_year$year))
T <- length(occasions)
occasion_index <- setNames(seq_len(T), occasions)
per_year$occasion <- occasion_index[as.character(per_year$year)]

# ---- (1) Constant CJS to separate survival and detection -------------------
histories <- split(per_year$occasion, per_year$tag)
first <- vapply(histories, min, integer(1))
last <- vapply(histories, max, integer(1))
n_capt <- vapply(histories, length, integer(1))

cjs_negloglik <- function(par) {
  phi <- plogis(par[1]); p <- plogis(par[2])
  if (phi <= 0 || phi >= 1 || p <= 0 || p >= 1) return(Inf)
  # chi[m] = P(not seen again | m occasions remain after current), m = 0..T-1
  chi <- numeric(T)            # index by (m+1); chi for m=0 is 1
  chi[1] <- 1
  for (m in seq_len(T - 1L)) chi[m + 1L] <- (1 - phi) + phi * (1 - p) * chi[m]
  open_span <- last - first                 # survival intervals
  seen_open <- n_capt - 1L                   # detections after first capture
  notseen_open <- open_span - seen_open
  ll <- sum(open_span) * log(phi) +
    sum(seen_open) * log(p) + sum(notseen_open) * log(1 - p) +
    sum(log(chi[(T - last) + 1L]))
  -ll
}
cjs <- optim(c(0, 0), cjs_negloglik, method = "BFGS", hessian = TRUE)
phi_hat <- plogis(cjs$par[1]); p_hat <- plogis(cjs$par[2])
cat(sprintf("Constant CJS: mean annual survival phi = %.3f, detection p = %.3f\n",
            phi_hat, p_hat))

# ---- (2) Size-dependent survival + growth from consecutive recaptures -------
captures <- data.frame(
  year = per_year$year, id = per_year$tag,
  observed_size = per_year$size, capture_probability = p_hat
)
events <- prepare_capture_recapture_events(captures, final_year = max(occasions),
                                           recapture_detection = p_hat)
cat(sprintf("Recapture events: %d origins, %d recaptured (annual return rate %.3f)\n",
            nrow(events), sum(events$recaptured), mean(events$recaptured)))

cr_params <- c("survival_at_80", "survival_slope_20",
               "growth_increment_80", "growth_slope_20", "log_growth_sd")
weak <- profile$prior
cr_prior <- list(mean = weak$mean[cr_params], sd = weak$sd[cr_params],
                 alternative_starts = NULL)
cr_log_post <- function(parameter) {
  parameter <- setNames(parameter, cr_params)
  lp <- sum(dnorm(parameter, cr_prior$mean, cr_prior$sd, log = TRUE))
  if (!is.finite(lp)) return(-Inf)
  capture_recapture_log_likelihood(parameter, events) + lp
}
cr_fit <- fit_laplace_custom(cr_log_post, cr_prior, "capture_recapture")
cr_draws <- sample_custom_laplace_draws(cr_fit, cr_log_post, draws = 2000L, seed = 20350001)

# ---- (3) Integrated: profile likelihood + recapture likelihood -------------
full_prior <- weak
int_log_post <- function(parameter) {
  parameter <- setNames(parameter, names(full_prior$mean))
  profile_lp <- inverse_ipm_log_posterior(parameter, profile$counts, mesh,
                                          full_prior, "gamma_poisson")
  if (!is.finite(profile_lp)) return(-Inf)
  profile_lp + capture_recapture_log_likelihood(parameter, events)
}
int_fit <- fit_laplace_custom(int_log_post, full_prior, "integrated")
int_draws <- sample_custom_laplace_draws(int_fit, int_log_post, draws = 2000L, seed = 20350002)

# ---- Assemble a comparison table -------------------------------------------
summ <- function(values) c(median = median(values),
                           lo = quantile(values, 0.025, names = FALSE),
                           hi = quantile(values, 0.975, names = FALSE))
survival_at <- function(draws, size) {
  z <- (size - 80) / 20
  plogis(draws[, "survival_at_80"] + draws[, "survival_slope_20"] * z)
}
growth_at <- function(draws, size) {
  z <- (size - 80) / 20
  pmax(0, draws[, "growth_increment_80"] + draws[, "growth_slope_20"] * z)
}
recruit_at <- function(draws, size) {
  z <- (size - 80) / 20
  exp(draws[, "log_recruitment_at_80"] + draws[, "recruitment_slope_20"] * z)
}
lambda_of <- function(draws) apply(draws, 1, function(p) {
  max(Re(eigen(local({ ipm <- build_centered_ipm(setNames(p, colnames(draws)), mesh)
    ipm$delta * ipm$kernel }), only.values = TRUE)$values))
})

profile_draws <- profile$draws
rows <- list()
add <- function(label, analysis, vals) rows[[length(rows) + 1L]] <<-
  data.frame(quantity = label, analysis = analysis,
             median = vals[1], lo = vals[2], hi = vals[3])
for (s in evaluation_sizes) {
  add(sprintf("survival_%d", s), "profile",    summ(survival_at(profile_draws, s)))
  add(sprintf("survival_%d", s), "recapture",  summ(survival_at(cr_draws, s)))
  add(sprintf("survival_%d", s), "integrated", summ(survival_at(int_draws, s)))
}
add("growth_increment_80", "profile",    summ(growth_at(profile_draws, 80)))
add("growth_increment_80", "recapture",  summ(growth_at(cr_draws, 80)))
add("growth_increment_80", "integrated", summ(growth_at(int_draws, 80)))
add("recruitment_80", "profile",    summ(recruit_at(profile_draws, 80)))
add("recruitment_80", "integrated", summ(recruit_at(int_draws, 80)))
add("lambda", "profile",    summ(lambda_of(profile_draws)))
add("lambda", "integrated", summ(lambda_of(int_draws)))
comparison <- do.call(rbind, rows)

write.csv(comparison, file.path(results_directory, "analysis_comparison.csv"), row.names = FALSE)
saveRDS(list(cjs = list(phi = phi_hat, p = p_hat), events = events,
             cr_fit = cr_fit, int_fit = int_fit, comparison = comparison),
        file.path(results_directory, "capture_recapture_fits.rds"), compress = "xz")

cat("\n=== West Brook: profile-only vs mark-recapture vs integrated ===\n")
cat(sprintf("(detection p-hat = %.3f from constant CJS; %d recapture events)\n\n",
            p_hat, nrow(events)))
fmt <- function(r) sprintf("%.3f [%.3f, %.3f]", r$median, r$lo, r$hi)
for (q in unique(comparison$quantity)) {
  cat(sprintf("%-20s\n", q))
  sub <- comparison[comparison$quantity == q, ]
  for (i in seq_len(nrow(sub))) cat(sprintf("   %-11s %s\n", sub$analysis[i], fmt(sub[i, ])))
}

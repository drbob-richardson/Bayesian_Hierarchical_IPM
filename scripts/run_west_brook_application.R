# Real-data application 2: profile-only inverse IPM on the West Brook brook
# trout PIT-tag dataset (USGS ScienceBase, DOI 10.5066/P14PDHXM, CC0). This is
# the long-series counterpart to the short chub application: one population
# (brook trout, the whole West Brook network, season 2) observed over 20 annual
# occasions (2000-2019), i.e. 19 transitions. The population is the whole stream
# network (all reaches), not the main stem alone, because ~7% of fish move
# between the main stem and tributaries each year and a main-stem-only boundary
# counts those movers as deaths. Untagged captures (mostly sub-PIT-size
# young-of-year) ARE included here: they carry the recruitment signal in the
# size profile. Individual identities are used only in the separate
# mark-recapture / integrated script.

source("R/simulate_population.R")
source("R/observe_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")
source("R/monte_carlo.R")

process_model <- Sys.getenv("WB_PROCESS", "gamma_poisson")
species <- Sys.getenv("WB_SPECIES", "bkt")
season <- as.integer(Sys.getenv("WB_SEASON", "2"))
year_min <- as.integer(Sys.getenv("WB_YEAR_MIN", "2000"))
year_max <- as.integer(Sys.getenv("WB_YEAR_MAX", "2019"))
# Observation window starts at the reliably-detected size (~PIT-tag size).
# Young-of-year below this are dropped; recruitment is modelled as entry into
# the observable size class (recruit-size prior re-centred into the window).
mesh_min <- as.integer(Sys.getenv("WB_MESH_MIN", "70"))
recruit_center <- as.numeric(Sys.getenv("WB_RECRUIT_MEAN", "80"))
mesh <- seq(mesh_min, 225, by = 5)
results_directory <- "results/west_brook_application"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv("data/west_brook/cdWB_electro_DR.csv", stringsAsFactors = FALSE)
# Whole network (all reaches), one standardized annual season for the profile.
keep <- raw$species == species &
  raw$season == season & raw$year >= year_min & raw$year <= year_max &
  !is.na(raw$observedLength)
b <- raw[keep, ]
# Restrict to the mesh support (drops the ~1% of fish above 227.5 mm).
b <- b[b$observedLength >= min(mesh) - 2.5 & b$observedLength <= max(mesh) + 2.5, ]
# Dedupe tagged fish caught more than once within a season-year; keep all
# untagged captures (each is a distinct individual we cannot collapse).
b$size <- b$observedLength
tagged <- !(is.na(b$tag) | b$tag == "")
b_tag <- aggregate(size ~ tag + year, data = b[tagged, ], FUN = mean)
b_untag <- b[!tagged, c("year", "size")]
profile_rows <- rbind(b_tag[, c("year", "size")], b_untag)

cat(sprintf("West Brook %s (all reaches), season %d, %d-%d: %d profile records (%d tagged, %d untagged), %d years.\n",
            species, season, year_min, year_max, nrow(profile_rows),
            nrow(b_tag), nrow(b_untag), length(unique(profile_rows$year))))

counts <- make_profile_count_matrix(profile_rows, mesh)
transitions <- nrow(counts) - 1L
saveRDS(list(counts = counts, mesh = mesh, years = as.integer(rownames(counts))),
        file.path(results_directory, "profile_stream.rds"))

prior <- default_inverse_ipm_priors("weak", process_model)
# Recruitment = entry into the observable (>= mesh_min) class.
prior$mean["recruit_mean_80"] <- recruit_center
evaluation_sizes <- c(80, 110, 140)

cat("Fitting profile-only inverse IPM (", process_model, ",", transitions, "transitions)...\n")
fit <- fit_laplace_inverse_ipm(counts, mesh, prior, process_model)
draws <- sample_laplace_draws(fit, draws = 2000L, seed = 20340001)

param_summary <- data.frame(
  parameter = colnames(draws),
  map = fit$map[colnames(draws)],
  posterior_mean = colMeans(draws),
  posterior_sd = apply(draws, 2, sd),
  q025 = apply(draws, 2, quantile, 0.025),
  q500 = apply(draws, 2, quantile, 0.5),
  q975 = apply(draws, 2, quantile, 0.975),
  prior_sd = prior$sd[colnames(draws)],
  row.names = NULL
)
param_summary$contraction <- 1 - param_summary$posterior_sd / param_summary$prior_sd

derived_draws <- t(apply(draws, 1, function(p) {
  derived_inverse_ipm_quantities(setNames(p, colnames(draws)), mesh, evaluation_sizes)
}))
derived_summary <- data.frame(
  quantity = colnames(derived_draws),
  posterior_mean = colMeans(derived_draws),
  posterior_sd = apply(derived_draws, 2, sd),
  q025 = apply(derived_draws, 2, quantile, 0.025),
  q500 = apply(derived_draws, 2, quantile, 0.5),
  q975 = apply(derived_draws, 2, quantile, 0.975),
  row.names = NULL
)

correlation <- cov2cor(fit$covariance)
recruit_ridge <- correlation["log_recruitment_at_80", "recruitment_slope_20"]
information <- solve(fit$covariance)
info_eigen <- eigen((information + t(information)) / 2, symmetric = TRUE)
condition_number <- max(info_eigen$values) / max(min(info_eigen$values), 1e-12)

write.csv(param_summary, file.path(results_directory, "parameters.csv"), row.names = FALSE)
write.csv(derived_summary, file.path(results_directory, "derived.csv"), row.names = FALSE)
saveRDS(list(fit = fit, draws = draws, counts = counts, mesh = mesh, prior = prior,
             param_summary = param_summary, derived_summary = derived_summary,
             recruit_level_slope_correlation = recruit_ridge,
             condition_number = condition_number),
        file.path(results_directory, "profile_fit.rds"), compress = "xz")

cat("\n=== West Brook brook trout, profile-only fit ===\n")
cat("Transitions:", transitions, "| PD Hessian:", fit$positive_definite,
    "| info condition number:", format(condition_number, digits = 3), "\n")
cat("Posterior corr(recruitment level, slope):", round(recruit_ridge, 3), "\n")
cat("\nContraction (most to least identified):\n")
ps <- param_summary[order(-param_summary$contraction), c("parameter","posterior_mean","contraction")]
print(ps, row.names = FALSE, digits = 3)
cat("\nDerived vital rates (median [95%]):\n")
show <- derived_summary[derived_summary$quantity %in%
  c("survival_50","survival_80","survival_110","growth_increment_80",
    "recruitment_80","recruit_mean_80","lambda"), ]
for (i in seq_len(nrow(show))) {
  cat(sprintf("  %-20s %.3f [%.3f, %.3f]\n", show$quantity[i],
              show$q500[i], show$q025[i], show$q975[i]))
}

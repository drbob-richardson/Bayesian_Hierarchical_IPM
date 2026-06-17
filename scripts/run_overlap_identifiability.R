# Numerical companion to the decomposition-identifiability proposition.
# Claim: from profiles, the total per-capita contribution s(x)+r(x) of a source
# size is identified, but the split into survival-growth vs recruitment is
# identified only where the grown-size and recruit-size distributions have
# (near-)disjoint support. We vary the recruit-size location so that recruits
# move from clearly below the adult size distribution (separated; identifiable)
# into the bulk of it (overlapping; confounded), and measure the recruitment
# level/slope confounding and prior-to-posterior contraction as a function of
# the distributional overlap.

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/monte_carlo.R")

mesh <- seq(25, 170, by = 5)
recruit_means <- seq(35, 100, by = 5)
transitions <- 8L
results_directory <- "results/overlap_identifiability"
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

# Weak prior, widened/recentred on recruit mean so it is non-informative across
# the whole sweep (so the posterior reflects the likelihood, not the prior).
base_prior <- default_inverse_ipm_priors("weak", "poisson")
base_prior$mean["recruit_mean_80"] <- 65
base_prior$sd["recruit_mean_80"] <- 35

# Fixed initial profile (broad, ~6000 fish) so only the recruit location changes.
initial_intensity <- 6000 * dnorm(mesh, mean = 75, sd = 25)
initial_profile <- round(initial_intensity / sum(initial_intensity) * 6000)

overlap_coefficient <- function(truth) {
  ipm <- build_centered_ipm(truth, mesh)
  w <- abs(Re(eigen(ipm$delta * ipm$kernel)$vectors[, 1]))   # stable size distn
  w <- w / sum(w)
  rs <- exp(truth["log_recruit_sd"])
  b <- dnorm(mesh, truth["recruit_mean_80"], rs); b <- b / sum(b)
  sum(pmin(w, b))                                            # overlap coefficient in [0,1]
}

run_one <- function(rm) {
  truth <- inverse_ipm_truth()
  truth["recruit_mean_80"] <- rm
  truth["recruit_mean_slope_20"] <- 0          # constant recruit size for cleanliness
  set.seed(20370000 + rm)
  expected <- expected_profile_counts(truth, initial_profile, mesh, transitions)
  if (is.null(expected) || any(!is.finite(expected))) return(NULL)
  counts <- expected
  for (t in seq(2L, nrow(expected))) counts[t, ] <- rpois(ncol(expected), pmax(expected[t, ], 1e-9))
  counts[1L, ] <- initial_profile
  fit <- tryCatch(fit_laplace_inverse_ipm(list(counts), mesh, base_prior, "poisson"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  corr <- cov2cor(fit$covariance)["log_recruitment_at_80", "recruitment_slope_20"]
  post_sd <- sqrt(diag(fit$covariance))
  contraction <- 1 - post_sd / base_prior$sd[names(post_sd)]
  data.frame(
    recruit_mean = rm,
    overlap = overlap_coefficient(truth),
    recruit_level_slope_corr = corr,
    recruitment_contraction = contraction["log_recruitment_at_80"],
    survival_contraction = contraction["survival_at_80"],
    growth_contraction = contraction["growth_increment_80"]
  )
}

results <- do.call(rbind, lapply(recruit_means, run_one))
write.csv(results, file.path(results_directory, "overlap_results.csv"), row.names = FALSE)

# Figure: confounding and recruitment identifiability vs overlap.
png("paper/figures/overlap_identifiability.png", width = 1700, height = 800, res = 220)
par(mfrow = c(1, 2), mar = c(4.4, 4.4, 2.4, 1))
o <- order(results$overlap)
r <- results[o, ]
plot(r$overlap, abs(r$recruit_level_slope_corr), type = "b", pch = 19, col = "#d95f02",
     ylim = c(0, 1), xlab = "growth / recruit-size overlap", ylab = expression(group("|",cor(alpha[r],beta[r]),"|")),
     main = "(a) Recruitment confounding")
grid(col = "grey90")
plot(r$overlap, r$recruitment_contraction, type = "b", pch = 19, col = "#1b9e77",
     ylim = range(c(0, r$recruitment_contraction, r$survival_contraction)),
     xlab = "growth / recruit-size overlap", ylab = "prior-to-posterior contraction",
     main = "(b) Identifiability by rate")
lines(r$overlap, r$survival_contraction, type = "b", pch = 17, col = "#7570b3")
legend("topright", c("recruitment level", "survival"), col = c("#1b9e77", "#7570b3"),
       pch = c(19, 17), bty = "n", cex = 0.9)
grid(col = "grey90")
dev.off()

cat("=== Overlap-identifiability experiment ===\n")
print(results, row.names = FALSE, digits = 3)
cat("\nFigure: paper/figures/overlap_identifiability.png\n")

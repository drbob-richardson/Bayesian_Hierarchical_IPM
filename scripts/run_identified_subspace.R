# Identified-subspace analysis of the West Brook profile-only likelihood.
# We compute the observed Fisher information of the LOG-LIKELIHOOD at the mode
# (raw finite-difference Hessian -- NOT the regularized Laplace inverse, whose
# floor would destroy the flat directions), eigendecompose it, and read off
# which linear combinations of vital rates the profiles inform (stiff,
# large-eigenvalue directions) and which are prior-determined (flat,
# near-zero-eigenvalue directions). This makes the identified subspace explicit,
# the constructive form of Propositions 1-2.

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")
source("R/prior_information.R")

results_directory <- "results/west_brook_application"
profile <- readRDS(file.path(results_directory, "profile_fit.rds"))
counts <- profile$counts; mesh <- profile$mesh; prior <- profile$prior
pnames <- names(prior$mean); p <- length(pnames)
map <- profile$fit$map[pnames]

loglik <- function(theta) {
  theta <- setNames(theta, pnames)
  inverse_ipm_log_posterior(theta, counts, mesh, prior, "gamma_poisson") -
    inverse_ipm_log_prior(theta, prior)
}

# Finite-difference Hessian; per-parameter step on the Laplace scale.
h <- 1e-3 * pmax(sqrt(diag(profile$fit$covariance)), 1e-3)
f0 <- loglik(map)
H <- matrix(0, p, p)
for (i in seq_len(p)) {
  ei <- numeric(p); ei[i] <- h[i]
  H[i, i] <- (loglik(map + ei) - 2 * f0 + loglik(map - ei)) / h[i]^2
}
for (i in seq_len(p - 1L)) for (j in seq(i + 1L, p)) {
  ei <- numeric(p); ei[i] <- h[i]; ej <- numeric(p); ej[j] <- h[j]
  H[i, j] <- H[j, i] <- (loglik(map + ei + ej) - loglik(map + ei - ej) -
                           loglik(map - ei + ej) + loglik(map - ei - ej)) / (4 * h[i] * h[j])
}
info <- -(H + t(H)) / 2                 # observed data information
eig <- eigen(info, symmetric = TRUE)
vals <- pmax(eig$values, 0)
vecs <- eig$vectors
rownames(vecs) <- pnames

# Per-direction prior vs posterior variance and contraction.
prior_var_dir <- apply(vecs, 2, function(v) sum((v * prior$sd[pnames])^2))
post_prec_dir <- vals + 1 / prior_var_dir
post_var_dir <- 1 / post_prec_dir
contraction_dir <- 1 - sqrt(post_var_dir / prior_var_dir)

ord <- order(vals, decreasing = TRUE)
spectrum <- data.frame(
  direction = seq_len(p),
  eigenvalue = vals[ord],
  prior_sd_dir = sqrt(prior_var_dir[ord]),
  post_sd_dir = sqrt(post_var_dir[ord]),
  contraction = contraction_dir[ord]
)
loadings <- vecs[, ord, drop = FALSE]
colnames(loadings) <- paste0("dir", seq_len(p))

write.csv(spectrum, file.path(results_directory, "subspace_spectrum.csv"), row.names = FALSE)
write.csv(round(loadings, 3), file.path(results_directory, "subspace_loadings.csv"))

# Effective number of identified directions (how many directions the data
# meaningfully contract; threshold at contraction > 0.5).
n_identified <- sum(spectrum$contraction > 0.5)

cat("=== Identified-subspace analysis (West Brook profile-only) ===\n")
cat(sprintf("Effective # of data-identified directions (contraction>0.5): %d of %d\n",
            n_identified, p))
cat("\nSpectrum (stiff -> flat):\n")
print(round(spectrum, 3), row.names = FALSE)
cat("\nStiffest (most data-informed) direction loadings:\n")
print(round(sort(loadings[, 1], decreasing = TRUE), 2))
cat("\nFlattest (prior-determined) direction loadings:\n")
print(round(sort(loadings[, p], decreasing = TRUE), 2))

# Figure: spectrum + flattest-direction loadings.
png("paper/figures/identified_subspace.png", width = 1900, height = 800, res = 220)
par(mfrow = c(1, 2), mar = c(7.5, 4.5, 2.6, 1))
cols <- ifelse(spectrum$contraction > 0.5, "#1b9e77", "#d95f02")
barplot(spectrum$contraction, names.arg = seq_len(p), col = cols, ylim = c(0, 1),
        xlab = "information eigen-direction (stiff -> flat)",
        ylab = "prior-to-posterior contraction", main = "(a) Identified subspace")
abline(h = 0.5, lty = 2, col = "grey40")
legend("topright", c("data-informed", "prior-determined"), fill = c("#1b9e77", "#d95f02"), bty = "n", cex = 0.9)
flat <- loadings[, p]
flat <- flat[order(-abs(flat))]
par(mar = c(9, 4.5, 2.6, 1))
barplot(flat, las = 2, col = "#d95f02", ylab = "loading", cex.names = 0.7,
        main = "(b) Flattest direction")
dev.off()
cat("\nFigure: paper/figures/identified_subspace.png\n")

# Figures for the JABES manuscript: the four new result sets
# (confounding ridge, SBC calibration, observation robustness, West Brook
# application). Base-R graphics; writes PNGs to paper/figures/.

source("R/simulate_population.R")
source("R/numerical_ipm.R")
source("R/bayesian_inverse_ipm.R")

fig_dir <- "paper/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

mvn_draws <- function(map, covariance, n = 2000L, seed = 1L) {
  set.seed(seed)
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  transform <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), length(eig$values))
  d <- matrix(rnorm(n * length(map)), n) %*% t(transform)
  d <- sweep(d, 2, map, "+")
  colnames(d) <- names(map)
  d
}

# ===========================================================================
# Figure 1: recruitment confounding ridge (2-D profile log-likelihood).
# ===========================================================================
tryCatch({
  slice <- as.matrix(read.csv(file.path("results/confounding_slice/loglik_slice.csv"),
                              row.names = 1, check.names = FALSE))
  level <- as.numeric(rownames(slice)); slope <- as.numeric(colnames(slice))
  # Likelihood is extremely peaked (large, low-noise dataset); cap the relative
  # log-likelihood so the flat ridge near the optimum is visible.
  zfloor <- -40
  z <- pmax(slice - max(slice), zfloor)
  truth <- inverse_ipm_truth()
  png(file.path(fig_dir, "confounding_ridge.png"), width = 1500, height = 1300, res = 240)
  par(mar = c(4.5, 4.5, 2, 1))
  filled.contour(level, slope, z, zlim = c(zfloor, 0), nlevels = 22,
    color.palette = function(n) hcl.colors(n, "YlGnBu", rev = TRUE),
    xlab = expression(log~recruitment~intensity~(alpha[r])),
    ylab = expression(recruitment~slope~(beta[r])),
    key.title = title(main = expression(Delta*"log L"), cex.main = 0.8),
    plot.axes = {
      axis(1); axis(2)
      contour(level, slope, z, levels = c(-1, -3, -8, -20), add = TRUE, col = "grey25", labcex = 0.7)
      points(truth["log_recruitment_at_80"], truth["recruitment_slope_20"], pch = 4, lwd = 3, col = "red", cex = 1.6)
    })
  dev.off()
  cat("Figure 1 (confounding ridge) written.\n")
}, error = function(e) cat("Fig1 failed:", conditionMessage(e), "\n"))

# ===========================================================================
# Figure 2: SBC rank histograms, informative (calibrated) vs weak (overconfident).
# ===========================================================================
tryCatch({
  ri <- read.csv("results/sbc/informative_gamma_poisson/ranks.csv", check.names = FALSE)
  rw <- read.csv("results/sbc/weak_gamma_poisson/ranks.csv", check.names = FALSE)
  show <- c("survival_at_80", "growth_increment_80", "recruitment_slope_20", "log_process_precision")
  labs <- c("survival intercept", "growth increment", "recruitment slope", "process precision")
  L <- 99
  png(file.path(fig_dir, "sbc_ranks.png"), width = 2000, height = 1100, res = 230)
  par(mfrow = c(2, 4), mar = c(3.4, 3.4, 2.2, 0.8), mgp = c(2, 0.6, 0))
  draw_row <- function(R, tag) for (j in seq_along(show)) {
    h <- hist(R[[show[j]]], breaks = seq(-0.5, L + 0.5, length.out = 21), plot = FALSE)
    barplot(h$counts, space = 0, col = "grey75", border = "grey55",
            main = sprintf("%s\n%s", tag, labs[j]), cex.main = 0.9,
            xlab = "rank statistic", ylab = "count")
    exp_per_bin <- nrow(R) / 20
    abline(h = exp_per_bin, col = "red", lwd = 2, lty = 2)
  }
  draw_row(ri, "informative prior")
  draw_row(rw, "weak prior")
  dev.off()
  cat("Figure 2 (SBC ranks) written.\n")
}, error = function(e) cat("Fig2 failed:", conditionMessage(e), "\n"))

# ===========================================================================
# Figure 3: observation robustness -- coverage collapse under size-dependent detection.
# ===========================================================================
tryCatch({
  grid <- read.csv("results/detection_robustness/task_grid.csv")
  rows <- lapply(grid$task_id, function(id) {
    r <- readRDS(sprintf("results/detection_robustness/tasks/task_%04d.rds", id))
    if (is.null(r$summary)) return(NULL)
    d <- r$summary$derived
    data.frame(regime = r$task$regime, process = r$task$process_model,
               param = r$summary$diagnostics$parameter_coverage,
               surv80 = d$covers_truth[d$quantity == "survival_80"])
  })
  det <- do.call(rbind, rows)
  det <- det[det$process == "gamma_poisson", ]
  agg <- aggregate(cbind(param, surv80) ~ regime, det, mean)
  ord <- c("constant", "time_varying_effort", "size_dependent")
  agg <- agg[match(ord, agg$regime), ]
  png(file.path(fig_dir, "detection_robustness.png"), width = 1500, height = 1100, res = 240)
  par(mar = c(5, 4.5, 2, 1))
  m <- t(as.matrix(agg[, c("param", "surv80")]))
  bp <- barplot(m, beside = TRUE, names.arg = c("constant", "time-varying\neffort", "size-\ndependent"),
                col = c("grey65", "steelblue"), ylim = c(0, 1), ylab = "95% interval coverage",
                legend.text = c("all parameters (mean)", "survival at 80 mm"),
                args.legend = list(x = "topright", bty = "n"))
  abline(h = 0.95, lty = 2, col = "red"); box()
  dev.off()
  cat("Figure 3 (detection robustness) written.\n")
}, error = function(e) cat("Fig3 failed:", conditionMessage(e), "\n"))

# ===========================================================================
# Figure 4: West Brook -- profile vs mark-recapture vs integrated.
# ===========================================================================
tryCatch({
  prof <- readRDS("results/west_brook_application/profile_fit.rds")
  crint <- readRDS("results/west_brook_application/capture_recapture_fits.rds")
  comp <- read.csv("results/west_brook_application/analysis_comparison.csv")
  pd <- prof$draws
  cd <- mvn_draws(crint$cr_fit$map, crint$cr_fit$covariance, 2000L, 11)
  id <- mvn_draws(crint$int_fit$map, crint$int_fit$covariance, 2000L, 12)
  sizes <- seq(70, 165, by = 2)
  surv <- function(D, s) plogis(D[, "survival_at_80"] + D[, "survival_slope_20"] * (s - 80) / 20)
  grow <- function(D, s) pmax(0, D[, "growth_increment_80"] + D[, "growth_slope_20"] * (s - 80) / 20)
  band <- function(D, fn, col) {
    M <- sapply(sizes, function(s) quantile(fn(D, s), c(.025, .5, .975)))
    polygon(c(sizes, rev(sizes)), c(M[1, ], rev(M[3, ])), col = adjustcolor(col, 0.18), border = NA)
    lines(sizes, M[2, ], col = col, lwd = 2.3)
  }
  cols <- c(profile = "#1b9e77", recapture = "#d95f02", integrated = "#7570b3")
  png(file.path(fig_dir, "west_brook_comparison.png"), width = 2100, height = 800, res = 220)
  par(mfrow = c(1, 3), mar = c(4.3, 4.3, 2.4, 1))
  # A: survival
  plot(NA, xlim = range(sizes), ylim = c(0, 0.85), xlab = "size (mm)", ylab = "annual survival",
       main = "(a) Survival"); grid(col = "grey90")
  band(pd, surv, cols["profile"]); band(cd, surv, cols["recapture"]); band(id, surv, cols["integrated"])
  legend("topright", names(cols), col = cols, lwd = 2.3, bty = "n", cex = 0.9)
  # B: growth
  plot(NA, xlim = range(sizes), ylim = c(0, 45), xlab = "size (mm)", ylab = "annual growth increment (mm)",
       main = "(b) Growth"); grid(col = "grey90")
  band(pd, grow, cols["profile"]); band(cd, grow, cols["recapture"]); band(id, grow, cols["integrated"])
  # C: lambda intervals
  lam <- comp[comp$quantity == "lambda", ]
  plot(NA, xlim = c(0.5, 2.5), ylim = range(c(lam$lo, lam$hi, 0.95, 1.1)), xaxt = "n",
       xlab = "", ylab = expression(lambda), main = "(c) Population growth")
  axis(1, at = 1:2, labels = lam$analysis); abline(h = 1, lty = 3, col = "grey50"); grid(col = "grey90", ny = NULL, nx = NA)
  for (i in seq_len(nrow(lam))) {
    arrows(i, lam$lo[i], i, lam$hi[i], angle = 90, code = 3, length = 0.08, lwd = 2.3,
           col = cols[lam$analysis[i]])
    points(i, lam$median[i], pch = 19, col = cols[lam$analysis[i]], cex = 1.3)
  }
  dev.off()
  cat("Figure 4 (West Brook comparison) written.\n")
}, error = function(e) cat("Fig4 failed:", conditionMessage(e), "\n"))

# ===========================================================================
# Figure 5: West Brook observed annual size profiles (near-stationary).
# ===========================================================================
tryCatch({
  ps <- readRDS("results/west_brook_application/profile_stream.rds")
  counts <- ps$counts; mesh <- ps$mesh; years <- ps$years
  prop <- counts / rowSums(counts)
  png(file.path(fig_dir, "west_brook_profiles.png"), width = 1500, height = 1150, res = 240)
  par(mar = c(4.3, 4.3, 2, 1))
  image(mesh, years, t(prop), col = hcl.colors(32, "YlOrRd", rev = TRUE),
        xlab = "size (mm)", ylab = "year", main = "West Brook annual size profiles (proportion)")
  box()
  dev.off()
  cat("Figure 5 (West Brook profiles) written.\n")
}, error = function(e) cat("Fig5 failed:", conditionMessage(e), "\n"))

cat("Done. Figures in", fig_dir, "\n")

# Biological truths and transport-bias scenarios.

transportability_truths <- function() {
  baseline <- inverse_ipm_truth()

  low_survival <- baseline
  low_survival["survival_at_80"] <- qlogis(0.40)

  flatter_fecundity <- baseline
  flatter_fecundity["recruitment_slope_20"] <- 0.60

  list(
    baseline = baseline,
    low_survival = low_survival,
    flatter_fecundity = flatter_fecundity
  )
}

fecundity_transport_bias <- function(multiplier) {
  stopifnot(
    length(multiplier) == 1L,
    is.finite(multiplier),
    multiplier >= 0
  )
  multiplier * c(
    log_recruitment_at_80 = 0.50,
    recruitment_slope_20 = -0.30
  )
}

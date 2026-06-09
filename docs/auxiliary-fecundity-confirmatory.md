# Confirmatory Auxiliary-Fecundity Experiment

## Why This Strengthens the Result

The earlier prior-strength experiment used abstract coefficient standard
deviations and only five profile datasets. This confirmatory experiment asks a
more directly interpretable design question:

> How many fish in an independent fecundity study are needed to resolve the
> recruitment confounding in repeated size profiles?

Prior information is generated from an independent study that observes
offspring counts across a balanced range of parent sizes. Its joint
recruitment-level and recruitment-slope posterior becomes the prior for the
profile-only inverse IPM. This retains the correlation between the two
fecundity coefficients.

## Design

- 30 newly simulated and independently generated profile datasets
- Three populations with four transitions each
- Initial population of 500 fish per population
- 25% profile detection
- Gamma-Poisson profile likelihood
- Independent fecundity studies containing 50, 200, or 800 fish
- Correctly transportable and systematically biased external studies
- Full and robust mixture borrowing
- Paired comparisons against the same profile dataset

All 180 analyses completed. Two Laplace Hessians from one difficult dataset
were non-positive-definite and were excluded from the corresponding biased
prior summaries.

## Confirmatory Results

| Analysis | Mean recruitment-curve RMSE | RMSE ratio versus profiles only | Recruitment-at-80 coverage |
|---|---:|---:|---:|
| Profiles only | 0.280 | 1.000 | 0.533 |
| Correct study, 50 fish | 0.190 | 0.679 | 0.433 |
| Correct study, 200 fish | 0.123 | 0.440 | 0.600 |
| Correct study, 800 fish | 0.084 | 0.299 | 0.867 |
| Biased study, 800 fish, full borrowing | 0.102 | 0.366 | 0.000 |
| Biased study, 800 fish, robust borrowing | 0.238 | 0.852 | 0.276 |

Bootstrap 95% intervals for the paired RMSE ratios were:

| Analysis | Paired RMSE-ratio interval |
|---|---:|
| Correct study, 50 fish | 0.502 to 0.921 |
| Correct study, 200 fish | 0.306 to 0.654 |
| Correct study, 800 fish | 0.206 to 0.444 |
| Biased study, full borrowing | 0.260 to 0.543 |
| Biased study, robust borrowing | 0.739 to 0.925 |

## Interpretation

The useful result is now stronger than “informative priors help.”

1. Auxiliary fecundity data improve recruitment recovery monotonically as the
   external study grows.
2. Two hundred independently studied fish reduce recruitment-curve error by
   more than half on average.
3. Eight hundred correctly studied fish reduce error by about 70% and raise
   recruitment-at-80 coverage from 53% to 87%.
4. A precise but non-transportable study can improve point-estimate RMSE while
   producing zero interval coverage. Accuracy alone would therefore give a
   dangerously favorable assessment.
5. Robust mixture borrowing reduces the damage from biased transport, but the
   profile data often cannot diagnose the conflict strongly enough. It also
   sacrifices much of the potential efficiency gain.

The biased-study result is particularly important. A strong external study can
pull estimates toward a curve that happens to predict the profile dynamics
well while being wrong for the individual-level vital-rate estimand.

## Robust Borrowing

The robust prior is an 80:20 mixture of the external-study posterior and the
diffuse fecundity prior. Because the posterior is genuinely multimodal, it is
approximated as a two-component Laplace mixture rather than by a single
posterior mode.

With the biased 800-fish study, the mean posterior probability assigned to the
external-study component was still about 0.67. Repeated profiles alone
therefore provided insufficient evidence to reliably reject the
non-transportable study.

This supports reporting the no-borrowing analysis alongside any informative
analysis. Robust borrowing is useful protection, but it is not a substitute
for assessing transportability scientifically.

## Important Limitations

- The external study directly samples effective offspring production and is
  idealized relative to real fecundity studies.
- Parent sizes are deliberately balanced across the observed range.
- Sex ratio, maturation, offspring detection, and reproductive observation
  error are not yet represented.
- Only the promising `3 x 4` profile design was used.
- Results use Laplace approximations; the robust mixture uses component-wise
  Laplace evidence approximations.
- Thirty datasets provide useful confirmation but not highly precise coverage
  estimates.

## Best Next Experiments

The easiest high-value extensions are:

1. vary transport bias continuously to estimate a borrowing safety boundary;
2. replace idealized offspring counts with realistic female-only or
   young-of-year observations;
3. compare auxiliary fecundity sampling with additional profile years at equal
   cost;
4. audit selected robust-mixture fits with full component-wise MCMC; and
5. increase the confirmatory endpoints to at least 100 datasets before final
   paper claims.

## Reproduction

```r
source("scripts/run_auxiliary_fecundity_experiment.R")
source("scripts/summarize_auxiliary_fecundity_experiment.R")
source("scripts/plot_auxiliary_fecundity_experiment.R")

source("scripts/run_confirmatory_auxiliary_fecundity.R")
source("scripts/summarize_confirmatory_auxiliary_fecundity.R")
source("scripts/plot_confirmatory_auxiliary_fecundity.R")
```

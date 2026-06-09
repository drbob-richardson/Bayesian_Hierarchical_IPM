# Transportability Boundary Across Biological Truths

## Question

The confirmatory auxiliary-fecundity experiment showed that a correctly
transportable study can resolve recruitment confounding, while a biased study
can create severe false certainty. This experiment asks:

> How much external-study bias can be tolerated before borrowing is less safe
> than analyzing the profiles alone, and is that boundary stable across
> different population biology?

## Design

Three individual-level biological truths were considered:

| Truth | Main difference |
|---|---|
| Baseline | Original generating vital rates |
| Lower survival | Survival at 80 mm reduced from 0.589 to 0.400 |
| Flatter fecundity | Recruitment slope per 20 mm reduced from 1.50 to 0.60 |

For each truth:

- ten fresh datasets were independently generated;
- each dataset contained three populations with four annual transitions;
- each population began with 500 fish;
- 25% of each annual profile was observed;
- an external fecundity study observed 800 fish; and
- profiles-only, full borrowing, and robust mixture borrowing were compared.

Transport bias was applied jointly to recruitment level and size slope:

```text
log recruitment at 80 mm: +0.50 * bias multiplier
recruitment slope per 20 mm: -0.30 * bias multiplier
```

At a bias multiplier of 0.25, the external population therefore differs by
only:

- 27% higher recruitment at 50 mm;
- 13% higher recruitment at 80 mm; and
- 1% higher recruitment at 110 mm.

The tested bias multipliers were 0, 0.25, 0.50, 0.75, 1.00, and 1.50.

The experiment included 390 fits. There were no fit failures and 98.2% of
Laplace Hessians were positive definite. One difficult baseline dataset had a
non-positive-definite profiles-only fit and corresponding robust-component
fits; those unmatched pairs were excluded from boundary summaries.

## Main Result

There is no universal transportability boundary.

Under full borrowing, even the smallest tested mismatch reduced
recruitment-at-80 coverage below the profiles-only baseline for all three
biological truths. Full borrowing continued to improve point-estimate RMSE in
many biased scenarios, but its intervals became severely miscalibrated.

Pooled across biological truths, full borrowing showed the distinction
especially clearly:

| Bias multiplier | RMSE ratio versus profiles only | Recruitment-at-80 coverage |
|---:|---:|---:|
| 0.00 | 0.356 | 0.862 |
| 0.25 | 0.396 | 0.483 |
| 0.50 | 0.390 | 0.103 |
| 0.75 | 0.442 | 0.000 |
| 1.00 | 0.532 | 0.000 |
| 1.50 | 0.760 | 0.000 |

Thus a transported study can continue improving point-estimate accuracy after
its uncertainty intervals have become unusable.

Under the deliberately conservative safety definition:

1. the paired RMSE-ratio 95% interval must remain below one; and
2. coverage must not fall below profiles-only coverage.

the largest tested safe biases were:

| Biological truth | Full borrowing | Robust borrowing |
|---|---:|---:|
| Baseline | 0.00 | 0.00 |
| Lower survival | 0.00 | 0.50 |
| Flatter fecundity | 0.00 | 0.75 |

These boundaries are exploratory because each point currently uses only nine
or ten paired datasets.

## Accuracy Versus Calibration

The experiment reinforces that RMSE alone is an unsafe performance criterion.

For baseline biology, full borrowing at bias multiplier 1.0 retained a mean
recruitment-curve RMSE ratio of 0.49 relative to profiles alone, but
recruitment-at-80 coverage was zero.

For lower-survival biology, full borrowing at bias multiplier 1.5 retained an
RMSE ratio of 0.58, but coverage was again zero.

For flatter-fecundity biology, full borrowing became worse even for point
estimation: its RMSE ratio rose above one at bias multiplier 1.0 and reached
2.0 at multiplier 1.5.

## What Robust Borrowing Can And Cannot Do

Robust borrowing used an 80:20 mixture of the external-study posterior and the
diffuse fecundity prior. It protected point-estimate accuracy over a wider bias
range and generally prevented the complete collapse seen under full
borrowing.

However, its ability to diagnose non-transportability was strongly dependent
on the true biology:

- For flatter fecundity, posterior external-study weight fell from about 0.91
  with no bias to effectively zero by bias multiplier 1.5.
- For baseline and lower-survival populations, external-study weight remained
  around 0.55 to 0.70 even under the largest tested bias.

Repeated profiles can therefore reject an incompatible external study only
when the disagreement changes profile dynamics along a sufficiently
likelihood-informed direction. Robust borrowing cannot guarantee safety when
the transport bias lies along the original confounding ridge.

## Paper-Level Interpretation

This supports a stronger and more nuanced claim:

> The value and safety of external vital-rate information depend on whether
> transport bias aligns with profile-likelihood-informed directions.
> Robust borrowing helps when profiles can detect the conflict, but cannot
> protect against bias along weakly identified demographic directions.

This links prior transportability directly to inverse-IPM identifiability.
Rather than asking whether external information is generally suitable, the
method should diagnose whether plausible between-population differences occur
along directions that the observed profiles can test.

## Limitations And Required Confirmation

- Each boundary point currently uses only ten datasets.
- Only one direction of transport bias was tested.
- The external study remains an idealized effective-offspring-count study.
- The robust mixture weight was fixed at 0.80.
- Results use component-wise Laplace approximations.
- Coverage differences of one or two datasets should not be overinterpreted.

Before final paper claims, the most important additions are:

1. increase replication around bias multipliers 0, 0.25, 0.50, and 0.75;
2. test transport bias directions aligned with the local sensitivity or
   likelihood-information eigenvectors;
3. audit selected full and robust fits using MCMC; and
4. repeat with a realistic auxiliary reproductive observation model.

## Reproduction

```r
source("scripts/run_transportability_boundary.R")
source("scripts/summarize_transportability_boundary.R")
source("scripts/plot_transportability_boundary.R")
```

# External Vital-Rate Information: First Experiment

## Question

Size-profile inverse IPMs are unusual because substantial information about
vital rates often exists before the profiles are analyzed. Previous
mark-recapture studies, life-history knowledge, related populations, and small
local auxiliary samples can all inform survival, growth, maturation,
fecundity, and recruit-size distributions.

The relevant question is therefore not whether priors help. It is:

> Which external information breaks the profile-likelihood confounding, how
> strong must it be, and how sensitive are conclusions to imperfect
> transportability?

## Focused Experiments

The first experiment fitted six prior-information regimes to 15 independently
simulated, data-poor datasets: five datasets for each of the `1 x 12`,
`3 x 4`, and `6 x 2` transition designs. Each dataset had an initial population
of 500 and 25% profile detection.

The regimes were:

- diffuse priors;
- external survival and growth information;
- external recruit-size distribution information;
- external fecundity information;
- external information on all vital rates; and
- systematically biased external information on all vital rates.

The external priors were centered on the individual-level generating rates and
were twice as uncertain as the earlier oracle priors. For example, their
approximate 95% intervals at 80 mm were:

| Quantity | External-prior 95% interval |
|---|---:|
| Survival | 0.31 to 0.82 |
| Growth increment | 0.6 to 16.2 mm |
| Recruitment | 0.057 to 0.596 recruits per fish |
| Recruit mean size | 36.8 to 52.4 mm |

The second experiment used five `3 x 4` datasets and varied full external-prior
uncertainty from one-half to four times the earlier oracle uncertainty. It
compared correctly centered and systematically biased external information.

Both experiments used Gamma-Poisson recovery and Laplace posterior
approximations. All 135 fitted posterior Hessians were positive definite.

## Main Results

Moderate external information changed most posterior results surprisingly
little. Across all quantities:

| Prior information | Mean 95% coverage | Mean relative absolute error |
|---|---:|---:|
| Diffuse | 0.747 | 0.179 |
| Survival and growth | 0.769 | 0.175 |
| Recruit-size distribution | 0.751 | 0.179 |
| Fecundity | 0.742 | 0.161 |
| All vital rates | 0.747 | 0.155 |
| All vital rates, biased | 0.769 | 0.165 |

This is not evidence that prior information is unimportant. The profile
likelihood is highly concentrated, including along directions where it targets
a population-effective pseudo-truth rather than the individual-level
generating rates. Moderately informative priors are therefore often
overwhelmed.

Information was most useful when it directly targeted the confounded
mechanism. Moderate fecundity information reduced recruitment-at-80 posterior
error by about 11%, while recruit-size information did not materially improve
recruitment intensity.

## How Strong Must External Information Be?

For the `3 x 4` design, diffuse-prior recruitment-at-80 error averaged 0.538
relative to the truth and 95% coverage was 0.60.

| Prior uncertainty multiplier | Correct-prior error / diffuse error | Correct coverage | Biased-prior error / diffuse error | Biased coverage |
|---:|---:|---:|---:|---:|
| 4.0 | 0.86 | 0.60 | 1.06 | 0.60 |
| 2.0 | 0.71 | 0.60 | 1.17 | 0.60 |
| 1.0 | 0.42 | 0.60 | 1.46 | 0.40 |
| 0.5 | 0.14 | 0.80 | 1.80 | 0.00 |

The strongest correctly centered prior had a recruitment-at-80 95% interval
of approximately 0.137 to 0.247 around a truth of 0.184. It largely resolved
the recruitment ridge. The equivalently strong biased prior tightly favored a
related-population value and failed badly.

Population growth `lambda` barely benefited from stronger priors. It was
already tightly constrained by profile dynamics, and stronger priors on
individual-level rates did not reliably move it closer to the individual-level
truth.

## Statistical Interpretation

External vital-rate information can make the inverse IPM useful, but it changes
the interpretation of the analysis:

1. Weakly identified vital rates may be learned primarily from external
   information, with profiles updating only likelihood-informed combinations.
2. Agreement between posterior rates and biological expectations is not
   evidence that profiles identified those rates.
3. External studies may estimate individual-level rates, while profile
   dynamics identify population-effective transition rates.
4. Strong borrowing is valuable only when transportability across populations,
   years, gears, and estimands is credible.

The methods paper should report prior-to-posterior contraction and posterior
sensitivity by demographic mechanism, not only a single weak-versus-strong
prior comparison.

## Recommended Prior Model

A mature version should replace fixed independent coefficient priors with a
hierarchical external-evidence model:

- express expert or literature information on vital-rate curves at meaningful
  sizes rather than on regression coefficients;
- estimate between-study and between-population heterogeneity;
- distinguish individual-level and population-effective rates;
- include robust mixture or heavy-tailed borrowing so local profiles can
  reject incompatible external information;
- perform prior predictive checks on population trajectories and size
  profiles; and
- report a no-borrowing analysis and prior-likelihood conflict diagnostics.

This creates a strong paper theme: inverse IPMs as a framework for combining
large anonymous size-profile datasets with smaller but biologically specific
sources of demographic information.

## Reproduction

```r
source("scripts/run_prior_information_experiment.R")
source("scripts/summarize_prior_information_experiment.R")
source("scripts/plot_prior_information_experiment.R")

source("scripts/run_prior_strength_gradient.R")
source("scripts/summarize_prior_strength_gradient.R")
source("scripts/plot_prior_strength_gradient.R")
source("scripts/plot_prior_predictive_information.R")
```

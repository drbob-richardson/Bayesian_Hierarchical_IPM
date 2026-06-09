# Broad Monte Carlo Exploration: First Grid

## Purpose

This experiment broadens the initial single-dataset analyses and asks how
inverse-IPM recovery changes with:

- allocation of a fixed twelve-transition budget;
- population size;
- fraction of fish observed in each annual profile;
- initial size-profile shape; and
- Poisson versus Gamma-Poisson profile likelihoods.

The main grid uses a Laplace posterior approximation for speed. Six
representative datasets were subsequently fitted with full MCMC to audit that
approximation.

## Scenario Grid

The grid contains:

| Dimension | Levels |
|---|---|
| Transition design | `1 x 12`, `3 x 4`, `6 x 2` |
| Initial population per trajectory | 500, 2500 |
| Fraction observed per annual profile | 0.25, 1.00 |
| Initial profile | baseline, small-fish cohort pulse |
| Independent replicates per scenario | 10 |
| Recovery likelihood | Poisson, Gamma-Poisson |

This produces:

- 240 independently simulated datasets;
- 480 Bayesian inverse-IPM fits;
- zero fit failures; and
- 98.96% positive-definite posterior Hessians.

Five Poisson fits with non-positive-definite Hessians were retained as
approximation-failure records but excluded from scientific coverage summaries.
Their Gaussian Laplace draws produced obviously invalid transformed tails.

## Laplace Approximation Audit

Six representative datasets spanning all trajectory designs and sparse versus
dense observations were refitted using four-chain adaptive-Metropolis MCMC.

Audit results:

| Diagnostic | Result |
|---|---:|
| Median absolute Laplace-MCMC mean difference | 0.106 MCMC posterior SD |
| Median Laplace/MCMC posterior-SD ratio | 0.988 |
| Parameter-coverage agreement | 100% |
| Maximum MCMC split R-hat | 1.063 |
| Minimum approximate MCMC ESS | 66 |

The Laplace approximation is adequate for broad screening in this grid.
Selected sparse or weakly identified scenarios still require full MCMC for
final inference.

## Overall Recovery

### Coverage by Transition Design

| Design | Poisson parameters | Poisson derived | Gamma-Poisson parameters | Gamma-Poisson derived |
|---|---:|---:|---:|---:|
| `1 x 12` | 0.734 | 0.701 | 0.764 | 0.718 |
| `3 x 4` | 0.795 | 0.775 | 0.818 | 0.818 |
| `6 x 2` | 0.795 | 0.813 | 0.816 | 0.836 |

Replicated trajectories outperform one long trajectory on average. However,
none of the designs achieve nominal 95% coverage because several quantities,
especially recruitment and long-run population growth, remain biased or
overconfident.

### Coverage by Data Volume

For Poisson recovery:

| Initial N | Fraction observed | Parameter coverage | Derived coverage |
|---:|---:|---:|---:|
| 500 | 0.25 | 0.755 | 0.756 |
| 2500 | 0.25 | 0.768 | 0.759 |
| 500 | 1.00 | 0.780 | 0.758 |
| 2500 | 1.00 | 0.795 | 0.777 |

More fish and more complete profiles help only modestly. The primary limitation
is not uncertainty about annual profile shapes. It is the decomposition of
profile change into demographic mechanisms.

## What Is Recoverable?

Across the full grid, Gamma-Poisson 95% interval coverage was approximately:

| Quantity | Coverage |
|---|---:|
| Recruit-size SD | 0.958 |
| Recruit mean at 110 mm | 0.929 |
| Recruit mean at 80 mm | 0.908 |
| Growth increment at 50 mm | 0.892 |
| Survival at 50 mm | 0.879 |
| Survival at 80 mm | 0.829 |
| Growth SD | 0.754 |
| Population growth `lambda` | 0.617 |
| Recruitment at 80 mm | 0.571 |

The size distribution of recruits is recoverable. Recruitment intensity is
not. Survival and mean growth are moderately recoverable, while variance
parameters and `lambda` remain less reliable.

## Recruitment Confounding

Recruitment level and its size slope have parameter coverage near 0.57 under
Gamma-Poisson and near 0.51–0.52 under Poisson. These parameters can compensate for
each other over the parent-size range that contributes most offspring.

This confounding propagates into low coverage for:

- recruitment across all evaluated sizes; and
- long-run population growth.

The next model-design work should focus on breaking this confounding rather
than simply increasing profile sample size.

## Replication Versus Long Follow-Up

With the same total number of transitions, multiple populations generally
improve coverage. Replication averages over the idiosyncratic demographic
history of a single population.

However, `6 x 2` does not identify every mechanism. Two transitions per
population still leave substantial recruitment-level/slope confounding. A
balanced design such as `3 x 4` appears promising, but a formal cost-based
design analysis is still needed.

## Cohort Pulse

Applying the same small-fish cohort pulse to every population did not uniformly
improve recovery. It particularly reduced coverage for the `3 x 4` design.

This is not evidence that excitation is unhelpful. It indicates that useful
excitation must create variation that distinguishes competing demographic
mechanisms. Repeating the same perturbation across every trajectory can simply
concentrate observations in a less informative part of the size domain.

The next excitation experiment should use deliberately diverse initial
profiles or observed environmental interventions across replicate populations.

## Gamma-Poisson Likelihood

Gamma-Poisson recovery modestly improves average coverage, especially for
recruitment and replicated designs. It does not resolve the main structural
confounding.

Independent overdispersion is useful for uncertainty calibration but is not a
replacement for demographic information.

## Important Limitations

- This grid uses one fixed set of individual-level vital-rate functions.
- Observation probability is constant across size and time.
- Populations are closed.
- No environmental variation or persistent individual heterogeneity is
  included.
- Ten replicates per scenario are sufficient for broad patterns, not precise
  coverage estimates.
- The target remains individual-level vital rates; a population-effective
  target may show different behavior.

## Recommended Next Grid

The next targeted experiment should focus on breaking recruitment confounding:

1. diverse initial size profiles across replicate populations;
2. known reproductive exclusion below maturation size;
3. small auxiliary fecundity or young-of-year samples;
4. mark-recapture benchmarks for survival and growth;
5. observed environmental variation affecting only selected vital rates; and
6. cost-based comparison of additional years versus additional populations.

## Reproduction

```r
source("scripts/run_monte_carlo_grid.R")
source("scripts/summarize_monte_carlo_grid.R")
source("scripts/audit_monte_carlo_mcmc.R")
source("scripts/plot_monte_carlo_grid.R")
```

# Research Roadmap

## Proposed Statistical Contribution

The paper should not claim that repeated size profiles generally identify
individual vital rates. Its contribution should be a framework for determining
what is identified, diagnosing confounding, and designing the additional data
or priors needed for useful inference.

A working title is:

> What Can Repeated Size Profiles Identify? Bayesian Inference and
> Confounding in Inverse Integral Projection Models

## Main Experiment

Generate populations from individual biological events rather than from the
recovery IPM:

1. Each individual has identity, size, age, sex, and optional latent quality.
2. Survival is generated as an individual Bernoulli event.
3. Conditional growth is generated as an individual size transition.
4. Reproduction is generated as an individual count process.
5. Recruit size is generated separately from parental size.
6. Optional movement, density dependence, environmental variation, and
   individual heterogeneity create controlled model misspecification.
7. Observation models independently produce repeated size profiles and
   mark-recapture histories.

The profile-only inverse IPM ignores identities. Mark-recapture data are used as
an external benchmark, or included only in explicitly integrated analyses.

## Estimands Must Be Declared First

There are four distinct inferential targets:

1. **Individual vital rates:** survival, growth, and reproductive output as
   functions of individual size.
2. **Population-scale effective rates:** redistribution rates after averaging
   over unobserved heterogeneity, movement, and detection.
3. **Kernel functionals:** population growth, stable size distribution,
   reproductive value, and elasticity.
4. **Prediction targets:** next-year and multi-year size profiles.

Good prediction does not imply that individual vital rates were recovered.
Simulation results must report these targets separately.

## Kernel Parameterization

Use a biologically decomposed positive kernel:

```text
K(x', x) = P(x', x) + F(x', x)
P(x', x) = s(x) g(x' | x)
F(x', x) = r(x) b(x' | x)
```

Require both conditional size distributions to integrate to one:

```text
integral g(x' | x) dx' = 1
integral b(x' | x) dx' = 1
```

Then:

```text
s(x) = integral P(x', x) dx'
r(x) = integral F(x', x) dx'
```

This removes arbitrary scaling between recruitment intensity and recruit-size
distribution. It does not, by itself, identify the decomposition of `K` into
`P` and `F`. That decomposition must be learned through temporal information,
biological restrictions, priors, or auxiliary data.

Do not estimate an unconstrained total kernel and attempt to label its pieces
as vital rates afterward.

## Positive Numerical Method

The recommended baseline is numerical quadrature on a fine size mesh:

1. Model survival and maturation through link functions.
2. Model realized recruitment with a log link.
3. Use normalized positive distributions for growth and recruit size.
4. Construct `K` directly from these positive pieces.
5. Propagate latent intensity using positive matrix multiplication.

This is transparent, stable, and avoids Fourier-basis negativity.

For flexible vital-rate curves:

- use nonnegative B-spline or M-spline bases with constrained coefficients;
- use I-splines when monotonicity is biologically justified;
- represent growth and recruit distributions as normalized positive mixtures;
- use a log-spline or log-Gaussian process only for process discrepancy, not as
  a substitute for the demographic kernel.

Mesh sensitivity and boundary-eviction corrections must be reported.

## Observation Model

The observation model should be explicit:

```text
observed intensity_t(x)
  = effort_t * detection_t(x) * latent population intensity_t(x)
```

Candidate likelihoods depend on survey design:

- point-process likelihood for unbinned locations or sizes;
- Poisson or negative-binomial likelihood for independent bin counts;
- binomial thinning when a finite population and detection probability are
  explicitly represented;
- repeated-pass removal or capture models when available.

Sampling effort, size-dependent detection, and measurement error can otherwise
be mistaken for demographic change.

## Simulation Design

Vary the following factors systematically:

| Dimension | Example levels |
|---|---|
| Number of annual transitions | 2, 4, 8, 16 |
| Individuals observed per year | 50, 200, 1000 |
| Number of replicate populations | 1, 4, 16 |
| Profile resolution | exact size, 2 mm, 10 mm bins |
| Detection | constant, size-dependent, time-varying |
| Sampling effort | known constant, known variable, unknown variable |
| Initial state | stable, perturbed, cohort pulse |
| Environment | constant, observed varying, latent varying |
| Vital-rate shape | parametric truth, spline truth, threshold truth |
| Open population processes | none, immigration, emigration |
| Individual heterogeneity | none, persistent latent quality |
| Recovery model | correct, mildly wrong, structurally wrong |

Variation and perturbation are especially important. A population repeatedly
observed near its stable size distribution may contain little information about
the operator that generated it. Cohort pulses, environmental covariates, and
replicate populations can provide the "excitation" needed to distinguish vital
rates.

## Prior Experiments

Treat prior information as part of the design, not as an afterthought.

For each vital-rate parameter or function:

1. Fit with weak, realistic, strong-correct, and strong-wrong priors.
2. Measure prior-to-posterior contraction.
3. Measure sensitivity to parameterization and prior scale.
4. Run simulation-based calibration under the assumed model.
5. Assess frequentist bias and coverage under model misspecification.
6. Identify parameters whose posterior is effectively the prior.

The paper should distinguish:

- structural nonidentifiability;
- practical weak identification;
- apparent identification caused by an informative but incorrect prior.

The first confirmatory auxiliary-fecundity experiment indicates that external
data should be expressed as an explicit study design whenever possible, rather
than only as abstract prior standard deviations. Study sample size,
transportability, and robust borrowing can then be evaluated directly. See
`docs/auxiliary-fecundity-confirmatory.md`.

The first multi-truth transportability experiment further indicates that no
single amount of external-study bias is universally safe. Robust borrowing can
reject conflict only when that conflict alters profile dynamics along a
likelihood-informed direction. See
`docs/transportability-boundary-first-run.md`.

## Competing Analyses

Each simulated dataset should support at least four comparisons:

1. **Profile-only inverse IPM**
2. **Mark-recapture analysis**
3. **Integrated IPM using both data streams**
4. **Oracle analysis using complete simulated individual histories**

Mark-recapture alone will not identify all reproductive processes, so compare
methods only on estimands each method can target. The simulator's known truth
remains the ultimate benchmark.

## Performance Metrics

Report more than parameter bias:

- bias, RMSE, interval width, and interval coverage for vital rates;
- integrated error over each vital-rate curve;
- error in the total kernel and in its survival/reproduction components;
- recovery of population growth and stable size distribution;
- one-step and multi-step predictive performance;
- posterior correlation and confounding manifolds;
- prior-to-posterior information gain;
- computational cost and failure rate.

## Identifiability Diagnostics

Develop diagnostics that can be used on real datasets:

- singular values of the local sensitivity or expected-information operator;
- likelihood-informed parameter directions;
- posterior correlation and ridge visualization;
- profile likelihoods or conditional posterior slices;
- prior-versus-posterior overlap;
- data-cloning or replicated-likelihood behavior;
- predictive equivalence of biologically different parameter sets.

An especially useful output would state which linear or nonlinear combinations
of vital rates are data-informed, even when individual parameters are not.

## Real-Data Program

Use real datasets in two roles:

1. **External validation:** fit profiles while withholding mark-recapture data,
   then compare inferred survival and growth with mark-recapture estimates.
2. **Realistic stress test:** use the observed sample sizes, effort variation,
   and time-series lengths to define simulation scenarios.

The southern leatherside chub data are valuable because both data streams
exist, but the short time series should be presented as a difficult case rather
than definitive proof of recovery. A second dataset with a longer series and
consistent effort would materially strengthen the paper.

## Recommended Paper Sequence

1. Formalize the positive profile-only model and estimands.
2. Establish simple analytical or numerical confounding examples.
3. Present broad simulation experiments.
4. Introduce identifiability and prior-sensitivity diagnostics.
5. Validate against withheld mark-recapture data.
6. Apply to the chub and a longer size-profile dataset.
7. Release reusable software and reproducible workflows.

## Main Risks

- Confusing effective population-scale rates with individual vital rates.
- Treating many fish within a year as many demographic transitions.
- Using priors to conceal rather than characterize nonidentifiability.
- Ignoring changing effort, detection, movement, or measurement error.
- Validating only under a recovery model identical to the generator.
- Claiming decomposition of a total kernel without information that separates
  survival-growth from reproduction-recruitment.

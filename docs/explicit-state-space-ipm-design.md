# Explicit Finite-Population State-Space IPM

## Transition Model

Let \(\mathbf{N}_{pt}\) be the integer-valued latent abundance by size class.
The IPM vital rates define separate survival-growth and fertility matrices,

\[
Q_\theta = \Delta P_\theta,
\qquad
F_\theta = \Delta R_\theta.
\]

For each source size class \(j\), survivors move among destination classes or
die:

\[
(S_{1j},\ldots,S_{mj},D_j)
\sim
\operatorname{Multinomial}\left[
N_{tj};
Q_{1j},\ldots,Q_{mj},1-\sum_iQ_{ij}
\right].
\]

Sex is not observed in the profiles. The transition therefore samples the
number of females in each source class,

\[
M_{tj}\sim\operatorname{Binomial}(N_{tj},0.5),
\]

and females produce at twice the population-average fertility rate:

\[
R_{ij}\mid M_{tj}
\sim
\operatorname{Poisson}(2M_{tj}F_{ij}).
\]

The next latent population is

\[
N_{t+1,i}=\sum_jS_{ij}+\sum_jR_{ij}.
\]

Marginally, this is a compound/binomial-Poisson recruitment process: the
number of contributing female parents is random and each contributing parent
produces Poisson offspring. It matches the individual-based fish simulator
more closely than a strict compound-Poisson clutch model with a Poisson number
of reproductive events. A strict clutch model can be added as a subsequent
recruitment-process alternative.

Observed profiles retain the Poisson-process observation layer:

\[
Y_{ti}\mid N_{ti},q\sim\operatorname{Poisson}(qN_{ti}).
\]

## Inference

An auxiliary particle filter integrates over the latent integer states. Its
look-ahead weights use the exact transition mean and marginal variance, while
the propagated particles always follow the explicit multinomial and
compound-Poisson transition. Particle marginal Metropolis-Hastings estimates
the ten vital-rate parameters.

An event-level importance proposal that tilted survivor destinations and
recruitment toward the upcoming profile was also tested. It retained a valid
likelihood correction but increased weight variance on the pilot data, so the
current default uses the exact transition after auxiliary ancestor selection.

## Initial Filtering Boundary

The transition simulator reproduced its theoretical marginal means within
0.8% and variances within 2.7% in a 10,000-transition check.

With 500 particles and the eight-population parent-excitation pilot, six
population-specific likelihood estimates had Monte Carlo standard deviations
between 0.21 and 0.54. Population 7 had standard deviation 13.16 and
population 8 had standard deviation 1.95. The two largest-parent populations
generate explosive recruitment, producing severe particle collapse.

The first estimation pilot therefore uses the six stable populations, whose
initial size distributions still span centers from approximately 60 to
117 mm. This establishes whether the explicit transition can estimate vital
rates within its current computationally feasible region. The explosive
recruitment cases define a target for improved conditional proposals,
correlated pseudo-marginal methods, or compiled inference.

## Run

```sh
EXPLICIT_STATE_POPULATIONS=6 \
EXPLICIT_STATE_PARTICLES=500 \
EXPLICIT_STATE_ITERATIONS=2000 \
EXPLICIT_STATE_WARMUP=1000 \
Rscript scripts/run_explicit_state_space_ipm_pilot.R
```

The correctly specified control experiment generates profiles from the same
explicit state process and can be run with:

```sh
EXPLICIT_MATCHED_PARTICLES=500 \
EXPLICIT_MATCHED_ITERATIONS=2000 \
EXPLICIT_MATCHED_WARMUP=1000 \
Rscript scripts/run_explicit_state_space_matched_pilot.R
```

Multiple matched-model chains are run with:

```sh
EXPLICIT_MATCHED_CHAINS=4 \
EXPLICIT_MATCHED_PARTICLES=500 \
EXPLICIT_MATCHED_ITERATIONS=2000 \
EXPLICIT_MATCHED_WARMUP=1000 \
Rscript scripts/run_explicit_state_space_matched_multichain.R
```

## Pilot Results

For the individual-based fish simulation, the explicit-state likelihood was
computationally stable for the first six parent-excitation populations. With
500 particles, the combined log-likelihood estimate had Monte Carlo standard
deviation 0.844. A short 600-iteration PMMH chain had acceptance 0.193 and
covered 7 of 10 vital-rate parameters, but its effective sample sizes were
only 8 to 21 and it was not treated as a converged result.

The correctly specified control generated profiles from the explicit
multinomial and compound-Poisson transition itself. Four PMMH chains used 500
particles and 1,500 iterations each. Results were:

- acceptance rates from 0.084 to 0.108;
- all 10 vital-rate parameter intervals covered truth;
- all selected derived vital-rate curves and population growth covered truth;
- minimum effective sample size 19;
- maximum split \(\widehat{R}=1.383\).

The poor \(\widehat{R}\) and effective sample size show that the current PMMH
implementation is not yet suitable for final inference. Growth increment and
recruit-size parameters mixed most slowly. Nevertheless, recovery under the
matched control demonstrates that the explicit transition is statistically
estimable. The remaining problem is computational efficiency.

The next implementation priorities are:

1. blockwise or gradient-informed parameter proposals;
2. correlated pseudo-marginal likelihood estimates;
3. an observation-conditioned transition bridge for explosive recruitment;
4. compiled transition and filtering code;
5. a binomial-detection observation alternative matching the individual
   simulator.

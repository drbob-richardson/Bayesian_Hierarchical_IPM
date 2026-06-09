# Particle-Filtered State-Space IPM

## Model

Let \(\boldsymbol{\Lambda}_{pt}\) be the latent population intensity over size
bins for population \(p\) at census \(t\). The vital-rate parameters
\(\theta\) generate the positive IPM transition matrix

\[
A_\theta = \Delta K_\theta
=
\Delta(P_\theta + F_\theta).
\]

The first estimable state-space prototype uses

\[
\begin{aligned}
\operatorname{E}[
  \boldsymbol{\Lambda}_{p,t+1}
  \mid \boldsymbol{\Lambda}_{pt},\theta
]
&= A_\theta\boldsymbol{\Lambda}_{pt},\\
\Lambda_{p,t+1,j}\mid\boldsymbol{\Lambda}_{pt},\theta,\psi
&\sim
\operatorname{Gamma}\left(
  \frac{[A_\theta\boldsymbol{\Lambda}_{pt}]_j}{\psi},
  \text{scale}=\psi
\right),\\
Y_{ptj}\mid\Lambda_{ptj},q
&\sim
\operatorname{Poisson}(q\Lambda_{ptj}).
\end{aligned}
\]

Thus the IPM-derived transition matrix controls the conditional mean of the
latent demographic process. The state Fano factor \(\psi\) controls process
variance, and the observed profile is a Poisson point-process discretization
conditional on the latent intensity.

The baseline latent state is conditioned on the first profile using
\(\boldsymbol{\Lambda}_{p0}=\mathbf{Y}_{p0}/q\). Later latent states are
integrated out with a bootstrap particle filter.

## Estimation

The particle filter provides an unbiased estimate of the likelihood
conditional on the baseline state. Particle marginal Metropolis-Hastings
(PMMH) uses this estimate to sample the posterior distribution of all ten
vital-rate parameters and \(\log\psi\).

Because the Gamma transition and Poisson observation are conjugate, the
implementation uses a fully adapted particle filter. Conditional on a
particle's predicted mean, the profile has an exact Gamma-Poisson predictive
distribution and the updated latent intensity has an exact Gamma conditional
distribution. This greatly reduces degeneracy relative to a bootstrap particle
filter.

This is a genuine state-space model. Unlike the one-step likelihood, observed
profiles after the baseline are not substituted for latent population states.
Information is propagated through the IPM transition across years, while each
new profile updates the filtering distribution.

## Current Approximation and Planned Extension

The Gamma transition is positive and computationally convenient, but it treats
destination-bin process innovations as conditionally independent. A later
particle filter can replace it with the finite-population transition:

- multinomial survival-growth movements derived from \(P_\theta\);
- compound-Poisson recruitment derived from \(F_\theta\);
- Poisson or binomial profile observation.

That extension changes the state simulator inside the particle filter without
changing the overall inferential architecture.

## Pilot

Run the initial pilot with:

```sh
STATE_SPACE_PARTICLES=300 \
STATE_SPACE_ITERATIONS=2000 \
STATE_SPACE_WARMUP=1000 \
Rscript scripts/run_state_space_ipm_pilot.R
```

Run multiple independent PMMH chains with:

```sh
STATE_SPACE_CHAINS=4 \
STATE_SPACE_PARTICLES=1000 \
STATE_SPACE_ITERATIONS=3000 \
STATE_SPACE_WARMUP=1500 \
Rscript scripts/run_state_space_ipm_multichain.R
```

The first 1,000-particle, 3,000-iteration single-chain pilot recovered 9 of 10
vital-rate parameters within their 95% posterior intervals. The recruitment
slope narrowly missed, with an upper interval bound of 1.491 versus truth
1.500. Derived survival at 50, 80, and 110 mm, growth increment at 80 mm,
recruitment at 80 and 110 mm, recruit mean size, and population growth rate
all covered truth. This is evidence that the state-space formulation can
estimate the vital-rate functions, but multi-chain and repeated-simulation
validation remain necessary.

The four-chain pilot used 1,000 particles and 3,000 iterations per chain.
Acceptance rates ranged from 0.103 to 0.180. All 10 vital-rate parameters and
all selected derived quantities covered their simulation truths. The minimum
effective sample size was approximately 54 and the maximum split
\(\widehat{R}\) was 1.103. The two largest \(\widehat{R}\) values belonged to
the recruit-size mean and slope, showing that further sampler tuning or longer
chains remain necessary.

This pilot is intentionally favorable: eight populations provide strong
variation in parent-size composition, detection is known, the initial latent
state is conditioned on the baseline profile, and the Gamma state transition
is an approximation to the individual-based generator. Its purpose is to
establish that a matrix-transition state-space IPM can estimate the vital-rate
functions before expanding the simulation study.

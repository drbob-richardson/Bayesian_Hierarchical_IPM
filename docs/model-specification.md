# Initial Model Specification

## Biological Generator

At census time `t`, individual `i` has size `x_it`, age `a_it`, sex `q_i`, and
optional latent quality `u_i`.

Survival:

```text
S_it ~ Bernoulli(s(x_it, a_it, q_i, u_i, E_t, N_t))
```

Growth conditional on survival:

```text
x_i,t+1 ~ G(. | x_it, a_it, q_i, u_i, E_t, N_t)
```

Realized recruits entering the next census:

```text
R_it ~ Count(r(x_it, a_it, q_i, u_i, E_t, N_t))
x_recruit,t+1 ~ B(. | x_it, E_t)
```

This generator is individual based. It need not imply that population
intensities exactly follow the fitted inverse IPM.

## Latent Population Recovery Model

Let `mu_t(x)` be latent population intensity over size. The demographic process
is:

```text
mu_t+1(x') = integral K_t(x', x) mu_t(x) dx
K_t(x', x) = s_t(x) g_t(x' | x) + r_t(x) b_t(x' | x)
```

The baseline implementation uses midpoint quadrature:

```text
mu_t+1 = delta_x * K_t %*% mu_t
```

All elements remain nonnegative when the kernel and starting intensity are
nonnegative.

Conditional normal transition distributions are integrated exactly over each
target-size bin using normal CDF differences. Evaluating densities only at bin
centers can materially inflate or deflate variance when bins are wide relative
to the biological transition variance.

## Observation Model

For survey effort `e_t` and size-dependent detection `q_t(x)`:

```text
lambda_t(x) = e_t q_t(x) mu_t(x)
Y_t ~ point process with intensity lambda_t(x)
```

For binned observations:

```text
Y_tj ~ Poisson(integral_bin_j lambda_t(x) dx)
```

Alternatives such as negative-binomial counts or finite-population binomial
sampling should be compared when appropriate.

## Recovering Vital Rates

If the survival-growth component `P` and reproduction-recruitment component
`F` are separately identified:

```text
s(x) = integral P(x', x) dx'
g(x' | x) = P(x', x) / s(x)
r(x) = integral F(x', x) dx'
b(x' | x) = F(x', x) / r(x)
```

If only total `K = P + F` is identified, these vital rates cannot generally be
recovered uniquely. Biological constraints can help:

- recruits occupy a restricted size range;
- surviving individuals cannot shrink beyond a plausible amount;
- maturation constrains which sizes can reproduce;
- external growth or survival data constrain selected components.

The consequences of each constraint must be evaluated through simulation.

## Open Decisions

- Census timing and whether reproduction uses size at `t` or `t+1`
- Demographic versus environmental process error
- Poisson, negative-binomial, or finite-population observation model
- Parametric, spline, or mixture vital-rate functions
- Computational backend: Stan, NIMBLE, TMB, or a combination
- Formal identifiability diagnostic and its uncertainty

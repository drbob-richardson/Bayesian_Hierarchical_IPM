# Chub Application: Profile-Only Inverse IPM (First Run)

Real-data application for the identifiability paper. Two southern leatherside
chub populations (`data/chub.csv`: Location, year, Size; 4,412 fish) are fitted
with the profile-only positive inverse IPM, Gamma–Poisson observation, weak
prior. The two populations are **not pooled** — biologists were explicit that
they cannot be combined — so each is fitted separately with its own vital-rate
parameters, as a single annual-profile trajectory.

`scripts/run_chub_application.R` → `results/chub_application/`.

| Population | Years | Transitions | n / year | Size range |
|---|---|---|---|---|
| Salina | 2003–2006 | 3 | 715, 911, 623, 541 | 39–144 mm |
| Lost Creek | 2003–2005 | 2 | 645, 669, 308 | 35–137 mm |

## Result: the simulation identifiability pattern reproduces on real data

Both fits had positive-definite Hessians. Posterior medians [95%]:

| Quantity | Salina | Lost Creek |
|---|---|---|
| Survival @ 80 mm | 0.786 [0.603, 0.900] | 0.293 [0.243, 0.348] |
| Growth increment @ 80 mm | 4.46 [1.67, 7.35] | 6.50 [4.88, 8.02] |
| Recruitment @ 80 mm | 0.113 [0.036, 0.365] | 0.404 [0.262, 0.627] |
| Recruit mean size @ 80 mm | 48.5 [45.4, 51.4] | 43.8 [42.2, 45.4] |
| λ (dominant eigenvalue) | 0.965 [0.822, 1.449] | 0.613 [0.533, 0.703] |

- **Identified:** recruit *size*, survival, and growth contract strongly from
  the prior. Recruit-mean-at-80 contraction 0.84 (Salina) / 0.92 (Lost Creek).
- **Confounded:** recruitment *intensity slope* is nearly unmoved from the
  prior (Salina contraction 0.185), and λ intervals are wide — for Salina the
  λ interval spans decline through strong growth, i.e. population growth is not
  identified from the profile series.
- **Recruitment level/slope ridge is present in the real data:** posterior
  correlation −0.708 (Salina, 3 transitions) and −0.863 (Lost Creek, 2
  transitions). The ridge is stronger for the shorter series, matching the
  simulation finding that short series retain recruitment confounding.
- Information condition numbers 370 (Salina) and 274 (Lost Creek).

## The populations are demographically distinct

Survival-at-80 differs by a factor of ~2.7 (0.79 vs 0.29) and λ by ~0.35
(0.97 vs 0.61) between populations. This is direct empirical support for the
biologists' position that the two cannot be pooled.

## Caveats (must accompany the section)

- **Size-dependent detection.** Chub are sampled by electrofishing. The
  observation-robustness run (`docs/supplementary-confirmatory-runs.md`)
  showed size-dependent detection collapses survival-at-80 coverage to 0. The
  survival and λ estimates here assume constant, size-independent detection
  (absorbed by the linear projection) and must be read with that risk.
- **Short series, closed population.** 2–3 transitions; no immigration/
  emigration, density dependence, or environmental variation is modeled.
- **No mark–recapture stream in this file.** The withheld-CR external
  validation requires individual capture histories, which `chub.csv` does not
  contain.

## Framing for the paper

This is a *diagnostic* application, not a recovery claim. It demonstrates that
the identifiability framework — contraction, the recruitment ridge, information
conditioning — behaves on real short-series data exactly as the simulations
predict: survival, growth, and recruit size are learnable; recruitment
intensity and population growth are not.

## Run

```sh
Rscript scripts/run_chub_application.R                 # gamma_poisson (default)
CHUB_PROCESS=poisson Rscript scripts/run_chub_application.R
```

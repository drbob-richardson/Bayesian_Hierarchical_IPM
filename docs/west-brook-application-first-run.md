# West Brook Application: Profile-only vs Mark-Recapture vs Integrated

Long-series real-data application for the identifiability paper, using the West
Brook (MA) brook trout PIT-tag dataset (USGS ScienceBase, DOI
10.5066/P14PDHXM, CC0; `data/west_brook/cdWB_electro_DR.csv`). This is the
long-series, two-stream counterpart to the short profiles-only chub
application.

Scripts: `scripts/run_west_brook_application.R` (profile stream + fit),
`scripts/run_west_brook_capture_recapture.R` (mark-recapture + integrated).
Output: `results/west_brook_application/`.

## Population and observation scoping (decisions and why)

- **Species/population:** brook trout, the **whole West Brook network** (all 4
  reaches), 2000-2019. A main-stem-only boundary was rejected after the data
  showed ~7% of fish (316/4570) move to tributaries each year, so main-stem
  "survival" was really site fidelity (apparent survival collapsed to 0.11).
  Whole-network apparent survival is 0.28 with detection 0.97.
- **Profile stream:** standardized season-2 cross-section, **sizes >= 70 mm**.
  Sub-70 mm young-of-year are dropped because they are below the PIT-tag size
  and electrofishing detects them only sparsely and variably -- including them
  injected the size-dependent-detection failure mode (see
  `docs/supplementary-confirmatory-runs.md`) and degenerated the fit (zero
  growth, reversed survival). Recruitment is therefore modelled as **entry into
  the observable (>= 70 mm) class** (recruit-size prior re-centred into the
  window).
- **Mark-recapture stream:** tagged fish only, one annual record per fish
  (earliest capture each year), detected in any reach/season. Detection p-hat =
  0.974 from a constant-survival/constant-detection CJS; this p-hat is plugged
  into the size-dependent survival+growth fit. Growth from recaptures is
  detection-free.

## Result: a plausible profile-only fit is not a recovered one

20 annual occasions (19 transitions); 20,077 recapture events. Posterior median
[95%]:

| Quantity | Profile-only | Mark-recapture | Integrated |
|---|---|---|---|
| Survival @ 80 mm | 0.43 [0.29, 0.60] | 0.31 [0.31, 0.32] | 0.31 [0.31, 0.32] |
| Survival @ 110 mm | 0.45 [0.35, 0.55] | 0.20 [0.20, 0.21] | 0.20 [0.19, 0.21] |
| Survival @ 140 mm | 0.46 [0.34, 0.58] | 0.12 [0.12, 0.13] | 0.12 [0.11, 0.12] |
| Growth increment @ 80 mm | 12.6 [9.3, 15.6] | 32.3 [31.7, 32.8] | 32.1 [31.5, 32.7] |
| lambda | 0.99 [0.97, 1.11] | -- | 1.016 [1.00, 1.05] |

- The profile-only fit is internally plausible (positive growth, moderate
  survival, lambda ~ 1, PD Hessian) but **wrong on both quantities the
  mark-recapture can check**:
  - **Growth underestimated**: 12.6 vs 32.3 mm/yr. A near-stationary size
    profile hides the true rate at which individuals move through it -- the
    stable-distribution blindness central to the paper.
  - **Size-dependent survival missed**: profile survival is roughly flat
    (~0.44) while mark-recapture shows a clear decline (0.31 -> 0.12), because
    profiles cannot see that large fish leave the network.
- **Recruitment intensity remains unidentified** by profiles (level/slope
  posterior correlation -0.995), as in the simulations and the chub.
- The **integrated** analysis adopts the mark-recapture survival/growth and
  tightens lambda to 1.016 [1.00, 1.05].

## Caveats (for the text)

- Mark-recapture "survival" is **apparent** survival: it folds in the
  size-dependent emigration of large fish out of the network, which is why it
  declines with size. It is the right benchmark for the *size pattern* the
  profile misses, not an absolute true-survival claim.
- The growth contrast rests on the recaptured (survivor, non-emigrant,
  detectable) sample.
- Letcher et al. (J. Anim. Ecol. 2015 integrated CMR; 2023 CJS growth) are the
  published external benchmark.

## Relation to the chub application

The chub (`docs/chub-application-first-run.md`) is a short (2-3 transition)
profiles-only case where confounding is severe and little is identified. West
Brook is the long (19-transition) two-stream case where the profile-only fit
*looks* successful yet is contradicted by mark-recapture -- a stronger and more
cautionary demonstration of the paper's thesis, and the natural headline
application.

## Run

```sh
Rscript scripts/run_west_brook_application.R          # profile-only (>=70mm)
Rscript scripts/run_west_brook_capture_recapture.R    # CR + integrated
```

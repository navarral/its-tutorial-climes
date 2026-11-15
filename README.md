# Interrupted time series for a policy: preterm births in Lombardy

This tutorial teaches interrupted time series (ITS) analysis in R with `brms`.
It uses simulated data for the preterm rate in 12 Lombardy provinces.
Each data point covers one province and one 4-week period.
No prior knowledge of ITS or Bayesian statistics is needed.

## What you will learn

- How ITS estimates a level change and a slope change after a policy start.
- How to simulate preterm births, with seasons and autocorrelation.
- How to fit an ITS model with `brms` and fast ADVI instead of MCMC.
- How to read results with the 95% credible interval, pd and ROPE.
- How to calculate power by simulation.
  - Window length before and after the policy start.
  - Births per period.
  - Size of the level change and the slope change.

## Files

| File | Content |
|------|---------|
| `its_tutorial_lombardy.R` | The full tutorial, with comments in 12 parts. |

## Requirements

- R 4.1 or later.
- A working C++ toolchain for Stan (see the [brms installation guide](https://paulbuerkner.com/brms/)).
- These R packages:

```r
install.packages(c("brms", "tidyverse", "bayestestR"))
```

## How to run

1. Download `its_tutorial_lombardy.R`.
2. Open it in RStudio.
3. Run the script from top to bottom, part by part.

The first model fit compiles Stan code and takes about one minute.
The full script runs in about 10 minutes.

## Tutorial outline

| Part | Topic |
|------|-------|
| 0-2 | Words, the ITS idea, and the Bayesian idea |
| 3-5 | Set up, simulate data, and plot the data |
| 6-7 | Design choices and the Bayesian model |
| 8 | Results: rate ratios, ROPE, counterfactual, residual check |
| 9 | Power by simulation for the level change and the slope change |
| 10 | Power check with `brms` |
| 11-12 | How to choose a design, and exercises |

## Important notes

- The data are simulated. Province sizes come from approximate birth counts.
- ADVI is an approximation and can understate uncertainty. Use MCMC for real analyses.
- The power grids use a fast `glm` model. Part 10 checks one design with `brms`.

## Sources

- [brms](https://paulbuerkner.com/brms/)
- [bayestestR articles](https://easystats.github.io/bayestestR/articles/)
- [bayestestR guidelines](https://easystats.github.io/bayestestR/articles/guidelines.html)

## License

GNU GENERAL PUBLIC LICENSE Version 3

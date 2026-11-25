ST920 | Childhood malaria in The Gambia
======================================
Dense Gaussian process inference using Stan and SPDE-INLA
Revised rectangular prediction domain

1. Getting started
------------------
Keep all seven R scripts in one folder, replacing 02_setup.R with the
latest rectangular-domain version. Open that folder as the RStudio
working directory. Install missing packages using the commented commands
in 01_dependencies.R. CmdStan and INLA are required for the full analysis.
The rectangular domain no longer requires GEOS >= 3.11; the old comment
about concave hulls in 01_dependencies.R can be ignored.

2. Script organisation
----------------------
01_dependencies.R    Load packages and provide installation instructions.
02_setup.R           Prepare village counts, priors, boundaries, grid and maps.
03_fit_stan.R        Fit or reload the dense GP binomial model.
04_fit_inla.R        Construct the mesh and fit the SPDE-INLA model.
05_analyse_stan.R    Diagnostics, prior/posterior plots and GP predictions.
06_analyse_inla.R    Summarise INLA predictions and produce maps.
07_compare_models.R  Compare predictions on the common grid.

3. Running the analysis
-----------------------
At the start of each R session:
  source("01_dependencies.R")
  source("02_setup.R")

Inspect images/gambia_prediction_domain.pdf before proceeding.
Copy an existing gambia_stan_fit.rds into the working directory and
results/stan_model_data.rds into results/ if available.

  source("03_fit_stan.R")
  source("04_fit_inla.R")
  source("05_analyse_stan.R")
  source("06_analyse_inla.R")
  source("07_compare_models.R")

The default refit_stan <- FALSE reuses a saved compatible fit. If no fit
exists, script 03 will sample again. For legacy fits without metadata,
verify the model, priors, observed data and village ordering yourself.
Changing only the prediction domain does not require refitting Stan.
Rebuild the INLA mesh and fit, and regenerate all predictions and comparisons.
Do not reuse prediction results from the old concave-hull grid.

4. Prediction domain and settings
---------------------------------
Edit the settings at the top of 02_setup.R:
  grid_spacing <- 5   # kilometres between prediction locations
  map_buffer <- 5     # kilometres beyond the full map
  M <- 1000L          # posterior predictive realisations
  threshold <- 0.30   # illustrative prevalence exceedance threshold

The grid covers the full extent of the villages and original geoR boundary
features, with a margin and limits rounded outwards to grid multiples.
All grid points are retained. NA-separated boundary segments, including
river features, are preserved for plotting. No external CRS is assumed.
The mesh is not constrained by the national boundary or river.
Both methods use the same prediction grid.

This is a rectangular study region, not a national land polygon. Predictions
may extend outside The Gambia and over water; interpret these as model-based
extrapolations, not validated estimates of malaria prevalence in those areas.

5. Stan prediction troubleshooting
----------------------------------
A larger grid increases memory and computation for joint GP simulations.
For a quick test, set M <- 10L after sourcing setup and before script 05.
Restore M <- 1000L for the final analysis and rerun script 05.

If using the original 05_analyse_stan.R, retain the previously supplied
numeric coercion fixes. At the start of matern_cov(), add:
  D <- as.matrix(D)
  sigma <- as.numeric(sigma)
  rho <- as.numeric(rho)
Extract sigma, rho, beta and the sampled field with as.numeric() in the
prediction loop to avoid non-conformable array errors.

After the prediction loop completes, verify all(is.finite(eta_draws)) is
TRUE, then run prob_draws <- plogis(eta_draws) before summarising prevalence.
Do not summarise incomplete predictive draws.

6. Outputs
----------
Figures are saved in images/. Intermediate prediction and fit objects are
saved in results/; the Stan fit and final Stan results also have RDS files
in the working directory. Preserve the fitted Stan model when updating
prediction settings. Regenerate downstream results on the new grid.

The revised README was checked against the latest setup and existing
scripts. The complete analysis has not been executed in R here.

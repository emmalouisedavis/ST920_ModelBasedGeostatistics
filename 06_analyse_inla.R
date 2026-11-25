# ST920 | Childhood malaria in The Gambia
# 06_analyse_inla.R

# Run 01 and 02 first; requires the saved INLA fit from 04.
inla_saved <- readRDS("results/gambia_inla_fit.rds")
stopifnot(identical(inla_saved$spec, analysis_spec))
fit <- inla_saved$fit
inla_grid <- grid

# Extract prevalence-scale summaries at the same grid used for Stan.
  ii <- inla_saved$ii
  fs <- fit$summary.fitted.values[ii, ]
  inla_grid$inla_mean <- fs$mean
  inla_grid$inla_sd <- fs$sd
  inla_grid$inla_lower <- fs[["0.025quant"]]
  inla_grid$inla_upper <- fs[["0.975quant"]]
  # The logit is monotone: P(p > t) = P(eta > logit(t)).
  inla_grid$inla_exceed <- vapply(ii, function(i) {
    1 - inla.pmarginal(qlogis(threshold), fit$marginals.linear.predictor[[i]])
  }, numeric(1))

p_inla_mean <- make_map(inla_grid, "inla_mean", "INLA posterior mean", "Prevalence", c(0, 1))
  p_inla_sd <- make_map(inla_grid, "inla_sd", "INLA posterior SD", "SD")
  p_inla_exceed <- make_map(inla_grid, "inla_exceed", "INLA exceedance probability",
                            "Probability", c(0, 1))
  ggsave("images/gambia_inla_prediction.pdf",
         p_inla_mean / p_inla_sd / p_inla_exceed, width = 10, height = 6)

saveRDS(list(grid = inla_grid, spec = analysis_spec,
             runtime_inla = inla_saved$runtime_inla),
        "results/gambia_inla_results.rds")

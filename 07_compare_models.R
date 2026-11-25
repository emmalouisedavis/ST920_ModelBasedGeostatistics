# ST920 | Childhood malaria in The Gambia
# 07_compare_models.R

# Run 01 and 02 first; requires the saved analyses from 05 and 06.
stan_results <- readRDS("results/gambia_stan_results.rds")
inla_results <- readRDS("results/gambia_inla_results.rds")
stopifnot(identical(stan_results$spec, analysis_spec),
          identical(inla_results$spec, analysis_spec))
comparison <- stan_results$grid
stopifnot(identical(comparison[, c("x", "y")],
                    inla_results$grid[, c("x", "y")]))
comparison$inla_mean <- inla_results$grid$inla_mean
comparison$inla_difference <- comparison$inla_mean - comparison$mean

# Compare prevalence on the same grid and colour scale.
p_stan <- make_map(comparison, "mean", "Dense GP: Stan posterior mean",
                   "Prevalence", c(0, 1))
p_inla <- make_map(comparison, "inla_mean", "SPDE-INLA posterior mean",
                   "Prevalence", c(0, 1))
p_diff <- signed_map(comparison, "inla_difference", "INLA minus Stan")
ggsave("images/gambia_mcmc_inla.pdf", p_stan / p_inla / p_diff,
       width = 10, height = 6)

# Differences include mesh approximation, INLA approximation and Monte Carlo error.
print(summary(comparison$inla_difference))
print(sqrt(mean(comparison$inla_difference^2)))

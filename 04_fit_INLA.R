# ST920 | Childhood malaria in The Gambia
# 04_fit_INLA.R

# Run 01 and 02 first. This fit does not depend on Stan.
mesh <- inla.mesh.2d(loc = coords, loc.domain = pred_coords,
                     max.edge = c(15, 50), cutoff = 3, offset = c(20, 60))
pdf("images/gambia_mesh.pdf", width = 10, height = 4.5)
plot(mesh, asp = 1)
for (seg in split(border, border$segment)) {
  lines(seg$x, seg$y, col = "grey35", lwd = 0.8)
}
points(coords, col = "red", pch = 16, cex = 0.6)
dev.off()

spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                            prior.range = c(rho0, prior_tail),
                            prior.sigma = c(sigma0, prior_tail))
# Projection matrices evaluate the mesh basis at observations/predictions.
A_obs <- inla.spde.make.A(mesh, loc = coords)
A_pred <- inla.spde.make.A(mesh, loc = pred_coords)
field_index <- inla.spde.make.index("field", n.spde = spde$n.spde)

# Stack observed responses and missing responses at prediction locations.
stack_obs <- inla.stack(
  data = list(Y = dat$Y, Ntrials = dat$N), A = list(A_obs, 1),
  effects = list(field_index, data.frame(intercept = rep(1, n))),
  tag = "observations"
)
stack_pred <- inla.stack(
  data = list(Y = rep(NA_integer_, n_pred), Ntrials = rep(1L, n_pred)),
  A = list(A_pred, 1),
  effects = list(field_index, data.frame(intercept = rep(1, n_pred))),
  tag = "prediction"
)
stack_full <- inla.stack(stack_obs, stack_pred)
stack_data <- inla.stack.data(stack_full)
formula <- Y ~ -1 + intercept + f(field, model = spde)

runtime_inla <- system.time({
  fit <- inla(
    formula, family = "binomial", data = stack_data,
    Ntrials = stack_data$Ntrials, control.family = list(link = "logit"),
    control.fixed = list(mean = 0, prec = 1 / beta_sd^2),
    control.predictor = list(A = inla.stack.A(stack_full), compute = TRUE, link = 1),
    control.compute = list(config = TRUE, return.marginals.predictor = TRUE)
  )
})
print(fit$summary.fixed)
print(fit$summary.hyperpar)
print(runtime_inla)

# Save the fit and prediction indices for independent analysis.
saveRDS(list(fit = fit, mesh = mesh,
             ii = inla.stack.index(stack_full, "prediction")$data,
             runtime_inla = runtime_inla, spec = analysis_spec),
        "results/gambia_inla_fit.rds")

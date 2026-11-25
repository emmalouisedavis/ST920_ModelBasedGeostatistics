# ST920 | Childhood malaria in The Gambia
# 05_analyse_stan.R

# Run 01 and 02 first; requires the saved Stan fit from 03.
fit_stan <- readRDS("results/gambia_stan_fit.rds")

# Check both covariance parameters and all latent spatial effects.
params <- c("beta", "sigma", "rho")

all_diag <- fit_stan$summary(
  variables = c(params, "S"),
  "mean",
  "median",
  "sd",
  ~ posterior::quantile2(
    .x, probs = c(0.025, 0.975)
  ),
  "rhat",
  "ess_bulk",
  "ess_tail"
)

param_summary <- all_diag |>
  filter(variable %in% params)
print(param_summary)

field_diagnostics <- all_diag |>
  filter(grepl("^S\\[", variable))

max(field_diagnostics$rhat)
min(field_diagnostics$ess_bulk)
min(field_diagnostics$ess_tail)

diag_hmc <- fit_stan$diagnostic_summary()
print(diag_hmc)

# Check for any problematic spatial effects (empty tibble is ideal)
field_diagnostics |>
  filter(
    rhat > 1.01 |
      ess_bulk < 400 |
      ess_tail < 400
  )

bad <- all_diag |> filter(!is.finite(rhat)|rhat>1.01|ess_bulk<400|ess_tail<400)
print(bad)

# Diagnostic figures

draws_params <- fit_stan$draws(variables = params)
ggsave("images/gambia_stan_trace.pdf",
       mcmc_trace(draws_params, pars = params) + theme_bw(),
       width = 10, height = 6)
ggsave("images/gambia_stan_rank.pdf",
       mcmc_rank_overlay(draws_params, pars = params), width = 10, height = 4)
ggsave("images/gambia_stan_pairs.pdf",
       mcmc_pairs(draws_params, pars = params), width = 8, height = 8)

# These thresholds are diagnostic guidelines, not proof of convergence.
if (nrow(bad) > 0 || any(diag_hmc$num_divergent > 0) ||
    any(diag_hmc$num_max_treedepth > 0) || any(diag_hmc$ebfmi < 0.3)) {
  stop("Review sampling diagnostics before prediction.")
}

# Extract joint draws, retaining the original village ordering.
post <- as.matrix(as_draws_matrix(fit_stan$draws(variables = c(params, "S"))))
field_columns <- paste0("S[", seq_len(n), "]")
stopifnot(all(field_columns %in% colnames(post)))

# Baseline prevalence is evaluated at S = 0, not averaged over the country.
print(quantile(plogis(post[, "beta"]), c(0.025, 0.5, 0.975)))
print(c(posterior_prob_short_range = mean(post[, "rho"] < rho0)))

# Compare prior and posterior distributions ----------------------------

prior_plot <- function(samples, prior, limits, title, xlab) {
  
  # Posterior samples
  posterior_df <- data.frame(
    x = as.numeric(samples)
  )
  
  # Evaluate the prior density over a regular grid
  prior_df <- data.frame(
    x = seq(limits[1], limits[2], length.out = 1000)
  )
  
  prior_df$density <- prior(prior_df$x)
  
  # Overlay prior and posterior densities
  ggplot() +
    geom_density(
      data = posterior_df,
      aes(x = x, colour = "Posterior", linetype = "Posterior"),
      linewidth = 1
    ) +
    geom_line(
      data = prior_df,
      aes(x = x, y = density, colour = "Prior", linetype = "Prior"),
      linewidth = 1
    ) +
    coord_cartesian(xlim = limits) +
    scale_colour_manual(
      values = c(Prior = "grey40", Posterior = "steelblue")
    ) +
    scale_linetype_manual(
      values = c(Prior = "dashed", Posterior = "solid")
    ) +
    labs(
      title = title,
      x = xlab,
      y = "Density",
      colour = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}

p_beta <- prior_plot(post[, "beta"], function(x) dnorm(x, 0, beta_sd),
                     c(-10, 10), "Intercept", expression(beta[0]))
p_sigma <- prior_plot(post[, "sigma"], function(x) dexp(x, lambda_sigma),
                      c(0, 4), "Spatial standard deviation", expression(sigma))
p_rho <- prior_plot(post[, "rho"],
                    function(x) lambda_rho / x^2 * exp(-lambda_rho / x),
                    c(1, 200), "Spatial range", expression(rho ~ "(km)"))
p_priors <- (p_beta | p_sigma | p_rho) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("images/gambia_prior_posterior.pdf", p_priors, width = 11, height = 3.5)

# Predict at unsampled locations ---------------------------------------

# Matern covariance (nu = 1), matching the Stan model.
matern_cov <- function(D, sigma, rho) {
  
  D <- as.matrix(D)
  sigma <- as.numeric(sigma)
  rho <- as.numeric(rho)
  
  stopifnot(
    length(sigma) == 1L,
    length(rho) == 1L,
    is.finite(sigma),
    is.finite(rho),
    sigma > 0,
    rho > 0
  )
  
  h <- sqrt(8) / rho * D
  
  K <- matrix(sigma^2, nrow(D), ncol(D))
  
  pos <- h > 0
  
  K[pos] <- sigma^2 * h[pos] * besselK(h[pos], nu = 1)
  
  return(K)
}

cross_dist <- function(A, B) {
  sqrt(outer(A[, 1], B[, 1], "-")^2 + outer(A[, 2], B[, 2], "-")^2)
}
Dso <- cross_dist(pred_coords, coords)
Dss <- as.matrix(dist(pred_coords))

# For each posterior draw, condition the GP at prediction locations on the
# sampled field at observed villages. Simulate jointly to preserve dependence.
# Cholesky solves avoid explicit matrix inversion.
set.seed(456)
M_used <- min(M, nrow(post))
stopifnot(M_used >= 3)
idx <- sample(seq_len(nrow(post)), M_used)
eta_draws <- matrix(NA_real_, n_pred, M_used)
runtime_prediction <- system.time({
  for (m in seq_len(M_used)) {
    j <- idx[m]
    sigma <- as.numeric(post[j, "sigma"])
    rho   <- as.numeric(post[j, "rho"])
    Koo <- matern_cov(D, sigma, rho)
    diag(Koo) <- diag(Koo) + jitter
    Kso <- matern_cov(Dso, sigma, rho)
    Kss <- matern_cov(Dss, sigma, rho)
    
    R <- chol(Koo)
    a <- backsolve(
      R,
      forwardsolve(t(R), as.numeric(post[j, field_columns]))
    )
    mu_star <- drop(Kso %*% a)
    V <- forwardsolve(t(R), t(Kso))
    Cstar <- Kss - crossprod(V)
    Cstar <- (Cstar + t(Cstar)) / 2
    diag(Cstar) <- diag(Cstar) + jitter
    
    S_star <- mu_star + drop(t(chol(Cstar)) %*% rnorm(n_pred))
    eta_draws[, m] <- as.numeric(post[j, "beta"]) + S_star
    if (m %% 100 == 0) message("Prediction draw ", m, " / ", M_used)
  }
})
print(runtime_prediction)

# Transform every draw to the prevalence scale before summarising.
prob_draws <- plogis(eta_draws)
saveRDS(list(grid = grid, eta_draws = eta_draws, prob_draws = prob_draws),
        "results/gambia_stan_predictions.rds")

# Summarise and visualise spatial predictions --------------------------

# Pointwise posterior summaries of underlying prevalence, not future counts.
stan_grid <- grid
stan_grid$mean <- rowMeans(prob_draws)
stan_grid$sd <- apply(prob_draws, 1, sd)
stan_grid$lower <- apply(prob_draws, 1, quantile, 0.025)
stan_grid$upper <- apply(prob_draws, 1, quantile, 0.975)
stan_grid$width <- stan_grid$upper - stan_grid$lower
stan_grid$exceed <- rowMeans(prob_draws > threshold)

# Show three possible surfaces alongside their posterior mean.
draw_ids <- 1:3
plot_df <- bind_rows(
  lapply(draw_ids, function(j) {
    mutate(stan_grid, prevalence = prob_draws[, j], panel = paste("Posterior draw", j))
  }),
  mutate(stan_grid, prevalence = mean, panel = "Posterior mean")
)
plot_df$panel <- factor(plot_df$panel,
                        levels = c(paste("Posterior draw", draw_ids), "Posterior mean"))
p_draws <- raster_map(plot_df, "prevalence", NULL) +
  facet_wrap(~panel, ncol = 2) +
  scale_fill_viridis_c(limits = c(0, 1), labels = label_percent(), name = "Prevalence") +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
ggsave("images/gambia_draws.pdf", p_draws, width = 11, height = 4)

# Maps of posterior mean, uncertainty and probability of elevated prevalence.
p_mean <- make_map(stan_grid, "mean", "Posterior mean prevalence", "Prevalence", c(0, 1))
p_sd <- make_map(stan_grid, "sd", "Posterior standard deviation", "SD")
p_exceed <- make_map(stan_grid, "exceed", paste("Probability that prevalence exceeds",
                                                percent(threshold)), "Probability", c(0, 1))
ggsave("images/gambia_prediction.pdf", p_mean / p_sd / p_exceed,
       width = 10, height = 6)

# Nonlinear transformation: E[logistic(eta)] is not logistic(E[eta]).
stan_grid$transform_first <- stan_grid$mean
stan_grid$average_first <- plogis(rowMeans(eta_draws))
stan_grid$difference <- stan_grid$transform_first - stan_grid$average_first
ggsave("images/gambia_nonlinear.pdf",
       signed_map(stan_grid, "difference", "Transform first minus average first"),
       width = 9, height = 3)
saveRDS(list(data = dat, border = border, coords = coords, grid = stan_grid,
             post = post, eta_draws = eta_draws, prob_draws = prob_draws,
             runtime_prediction = runtime_prediction, spec = analysis_spec), "results/gambia_stan_results.rds")

# ST920 | Childhood malaria in The Gambia
# 03_fit_stan.R

# Run 01 and 02 first. Reuse an existing fit unless refit_stan is TRUE.
stopifnot(beta_sd == 5, jitter == 1e-8) # Fixed constants in this Stan model.

# Non-centred representation: S = L z, z ~ N(0, I), LL' = covariance.
# Stan writes its model to a separate file automatically; no manual editing.
stan_code <- '
data {
  int<lower=1> n;
  array[n] int<lower=0> Y;
  array[n] int<lower=0> N;
  matrix[n,n] D;
  real<lower=0> lambda_sigma;
  real<lower=0> lambda_rho;
}
parameters {
  real beta;
  real<lower=0> sigma;
  real<lower=0> rho;
  vector[n] z;
}
transformed parameters {
  vector[n] S;
  {
    real kappa = sqrt(8.0)/rho;
    matrix[n,n] K;
    for (i in 1:n) {
      K[i,i] = square(sigma)+1e-8;
      if (i<n) {
        for (j in (i+1):n) {
          real h = kappa*D[i,j];
          if (h>0) K[i,j] = square(sigma)*h*modified_bessel_second_kind(1,h);
          else K[i,j] = square(sigma);
          K[j,i] = K[i,j];
        }
      }
    }
    S = cholesky_decompose(K)*z;
  }
}
model {
  beta ~ normal(0,5);
  sigma ~ exponential(lambda_sigma);
  rho ~ inv_gamma(1,lambda_rho);
  z ~ std_normal();
  Y ~ binomial_logit(N,beta+S);
}
generated quantities {
  vector[n] p = inv_logit(beta+S);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = binomial_logit_lpmf(Y[i]|N[i],beta+S[i]);
}
'

writeLines(stan_code, "gambia_gp.stan")
stan_data <- list(n = n, Y = as.integer(dat$Y), N = as.integer(dat$N),
                  D = D, lambda_sigma = lambda_sigma, lambda_rho = lambda_rho)

# Save the fit so that plotting and prediction do not require another HMC run.
# An older saved fit has no provenance record: verify its data ordering,
# coordinate units and priors yourself before reusing it.
fit_file <- "results/gambia_stan_fit.rds"
metadata_file <- "results/stan_model_data.rds"
metadata <- list(data = stan_data, model = stan_code)
if (!refit_stan && file.exists(fit_file)) {
  if (file.exists(metadata_file)) {
    stopifnot(identical(readRDS(metadata_file), metadata))
  } else {
    warning("Reusing an older fit: verify model, priors and village ordering.")
  }
  fit_stan <- readRDS(fit_file)
} else {
  mod <- cmdstan_model("gambia_gp.stan")
  runtime_stan <- system.time({
    fit_stan <- mod$sample(
      data = stan_data, seed = 123, chains = 4,
      parallel_chains = parallel_chains,
      iter_warmup = 1000, iter_sampling = 1000,
      adapt_delta = 0.95, max_treedepth = 12, refresh = 100
    )
  })
  print(runtime_stan)
  fit_stan$save_object(file = fit_file)
  saveRDS(metadata, metadata_file)
}
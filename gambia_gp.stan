
data {
  int<lower=1> n;

  array[n] int<lower=0> Y;
  array[n] int<lower=0> N;

  matrix[n, n] D;

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
    real kappa = sqrt(8.0) / rho;

    matrix[n, n] K;
    matrix[n, n] L;

    // Matern covariance, nu = 1
    for (i in 1:n) {

      K[i, i] = square(sigma) + 1e-8;

      if (i < n) {
        for (j in (i + 1):n) {

          real h = kappa * D[i, j];

          if (h > 0) {
            K[i, j] =
              square(sigma) * h *
              modified_bessel_second_kind(1, h);
          } else {
            K[i, j] = square(sigma);
          }

          K[j, i] = K[i, j];
        }
      }
    }

    L = cholesky_decompose(K);

    // Non-centred parameterisation
    S = L * z;
  }
}

model {
  // Prior distributions
  beta ~ normal(0, 5);

  sigma ~ exponential(lambda_sigma);

  rho ~ inv_gamma(1, lambda_rho);

  z ~ std_normal();

  // Binomial likelihood on logit scale
  Y ~ binomial_logit(N, beta + S);
}

generated quantities {
  vector[n] p;
  vector[n] log_lik;

  p = inv_logit(beta + S);

  for (i in 1:n) {
    log_lik[i] =
      binomial_logit_lpmf(
        Y[i] | N[i], beta + S[i]
      );
  }
}


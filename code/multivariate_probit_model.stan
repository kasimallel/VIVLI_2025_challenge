data {
  int N, P, D;                 // N = number of patients
  // P = number of covariates (e.g. gender, agegroup, department)
  // D = number of resistance outcomes (7 antibiotics)
  
  matrix[N, P] X;              // Design matrix: covariates (already processed in R)
  array[N, D] int<lower=0, upper=1> Y; // Resistance outcomes for each patient and drug
}
parameters {
  vector[D] alpha;            // Intercepts for each antibiotic/resistance outcome
  matrix[P, D] beta;  // A different coefficient vector for each outcome
  cholesky_factor_corr[D] L; // L captures the correlation structure among the 7 resistances
  matrix[N,D] Z_raw;         // Z_raw represents patient-specific latent variables for each resistance
}
transformed parameters {
  matrix[N,D] Z;
  for (n in 1:N) {
    Z[n] = alpha' + X[n] * beta + Z_raw[n] * L'; // The covariate effects (X * beta) are added to correlated noise (Z_raw * L') for each outcome.
    }
}

model {
  // Priors
  alpha ~ normal(0,1); // change this at one point to sperate priors potentially (maybe hierarchical prior, which would just mean putting an mu and sigma on alpha and draw from that)
  to_vector(beta) ~ normal(0,1);
  L ~ lkj_corr_cholesky(2);
  to_vector(Z_raw) ~ normal(0,1);
  // Likelihood
  for (n in 1:N)
    for (d in 1:D)
      Y[n,d] ~ bernoulli(Phi(Z[n,d]));
}

generated quantities {
  array[N, D] int Y_rep;
  matrix[N, D] log_lik;
  matrix[N, D] pred_prob;

  for (n in 1:N) {
    for (d in 1:D) {
      real p = Phi(Z[n, d]);
      Y_rep[n, d] = bernoulli_rng(p);
      log_lik[n, d] = bernoulli_lpmf(Y[n, d] | p);
      pred_prob[n, d] = p;
    }
  }
}

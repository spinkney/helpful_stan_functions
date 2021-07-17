functions {
 #include multi_skewed_t.stan
}
data {
  int<lower=1> p;
  int<lower=1> n;
  matrix[p, n] y;
}
parameters {
  vector[p] lambda;
  vector<lower=0>[p] sigma;
  vector[p] mu;
  real<lower=0> nu;
  cholesky_factor_corr[p] L;
}
model {
  mu ~ normal(0, 5);
  nu ~ exponential(1);
  sigma ~ gamma(2, 0.1);
  lambda ~ normal(0, 10);
  L ~ lkj_corr_cholesky(1);
  
  {
    matrix[p, n] y_norm;
    
    for (i in 1:p)
      y_norm[i] = (y[i] - mu[i]);
    
    y_norm ~ multi_skewed_t(diag_pre_multiply(sigma, L), sigma, lambda, nu);
  }
}
generated quantities {
  matrix[p, p] Sigma_raw = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma, L));
  matrix[p, p] Sigma = (nu / (nu - 2)) * Sigma_raw - quad_form_diag(mu * mu', sigma);
}
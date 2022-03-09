functions {
  #include centered_gaussian_copula.stanfunctions
}
data {
  int<lower=0> N; // number of observations
  array[3] int<lower=0> J; // total number of outcomes
  int<lower=0> K; // number of covariates for concatenated design matrices (X1, .., XK)
  matrix[N, J[1]] Yn; // normal outcomes
  array[N, J[2]] int<lower=0, upper=1> Yb; // bernoulli outcomes
  array[N, J[3]] int<lower=0> Yp; // poisson outcomes
  matrix[N, K] X; // concatenated design matrices (X1, ..., XK)
  array[3] int<lower=0> Kj; // J-dim integer array giving how many covariates per outcome
  int<lower=0> special;
}
transformed data {
  int J_all = sum(J);
  // Create separate design matrices for each outcome type
  matrix[N, Kj[1]] Xn = X[ : , 1:Kj[1]];
  matrix[N, Kj[2]] Xb = X[ : , (Kj[1] + 1):(Kj[1] + Kj[2])];
  matrix[N, Kj[3]] Xp = X[ : , (Kj[2] + 1):(Kj[2] + Kj[3])];
  vector[J_all] mu_zero = rep_vector(0, J_all);
}
parameters {
  vector[Kj[1]] beta_n; // Vector of normal regression coefficients
  vector[Kj[2]] beta_b; // Vector of bernoulli regression coefficients
  vector[Kj[3]] beta_p; // Vector of poisson regression coefficients
  cholesky_factor_corr[J_all] L; // Cholesky decomposition of JxJ correlation matrix
  vector<lower=0>[J[1]] sigmasq; // Jn-dim vector of variances (may be 0-dim)
  matrix<lower=0, upper=1>[N, J[2]] uraw_b; // latent variables for bernoulli outcomes
  matrix<lower=0, upper=1>[N, J[3]] uraw_p; // latent variables for Poisson outcomes
}
model {
  // initialize variables
  vector[J[1]] sigma = sqrt(sigmasq); // stdev
  // Calculate the means for each regression
  matrix[N, J[1]] mu_n = Xn * rep_matrix(beta_n, J[1]);
  matrix[N, J[2]] mu_b = Xb * rep_matrix(beta_b, J[2]);
  matrix[N, J[3]] mu_p = Xp * rep_matrix(beta_p, J[3]);

  sigmasq ~ inv_gamma(1.0e-4, 1.0e-4);
  // priors for regression coefficients
  beta_n ~ normal(0.0, 10.0);
  beta_b ~ normal(0.0, 10.0);
  beta_p ~ normal(0.0, 10.0);
  { normal_marginal(Yn, mu_n, sigma),
    bernoulli_marginal(Yb, mu_b, uraw_b),
    poisson_marginal(Yp, mu_p, uraw_p) } ~ centered_gaussian_copula_cholesky(L);
}
generated quantities {
  corr_matrix[J_all] Gamma = multiply_lower_tri_self_transpose(L);
}

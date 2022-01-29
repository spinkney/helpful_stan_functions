functions {
  #include mixed_copula.stanfunctions
}
data {
  int<lower=0> N; // number of observations
  array[3] int<lower=0> J; // total number of outcomes
  int<lower=0> K; // number of covariates for concatenated design matrices (X1, .., XK)
  array[N, J[1]] real Yn; // normal outcomes
  array[N, J[2]] int<lower=0, upper=1> Yb; // bernoulli outcomes
  array[N, J[3]] int<lower=0> Yp; // poisson outcomes
  matrix[N, K] X; // concatenated design matrices (X1, ..., XK)
  array[3] int<lower=0> Kj; // J-dim integer array giving how many covariates per outcome
  int<lower=0> special;
}
transformed data {
  int J_all = sum(J);
  array[2] matrix[J_all, N] bounds_ind = get_bounds_indicator(Yn, Yb, Yp, N,
                                                              J[1], J[2],
                                                              J[3], J_all);
  vector[J_all] mu_zero = rep_vector(0, J_all);
}
parameters {
  vector[K] beta; // long vector of all regression coefficients
  cholesky_factor_corr[J_all] L; // Cholesky decomposition of JxJ correlation matrix
  array[J[1]] real<lower=0> sigmasq; // Jn-dim array of variances (may be 0-dim)
  array[J[2]] vector<lower=0, upper=1>[N] uraw_b; // latent variables for bernoulli outcomes
  array[J[3]] vector<lower=0, upper=1>[N] uraw_p; // latent variables for Poisson outcomes
}
model {
  // initialize variables
  array[J[1]] real sigma; // stdev
  matrix[N, J_all] mu_glm = get_means(beta, X, J, Kj);
  
  // get standard deviations and declare prior on sigmasq if applicable
  if (J[1] > 0) {
    sigma = sqrt(sigmasq);
    sigmasq ~ inv_gamma(1.0e-4, 1.0e-4);
    
    for (j in 1 : J[1]) 
      to_vector(Yn[ : , j]) ~ normal(to_vector(mu_glm[ : , j]), sigma[j]);
  }
  // prior for beta
  beta ~ normal(0.0, 10.0);
  
  // log target increment
    if (special == 1) {
        target += mixed_cop_sp_lp(Yn, Yb, Yp,
                   uraw_b, uraw_p,
                   mu_glm,
                   sigma,
                   L,
                   bounds_ind,
                   mu_zero,
                   J, N);
    } else target += mixed_cop_lp(Yn, Yb, Yp,
                   uraw_b, uraw_p,
                   mu_glm,
                   sigma,
                   L,
                   bounds_ind,
                   mu_zero,
                   J, N);
}
generated quantities {
  corr_matrix[J_all] Gamma = multiply_lower_tri_self_transpose(L);
}

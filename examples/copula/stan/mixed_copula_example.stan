functions {
  #include mixed_copula.stan
}

data {
  int<lower=0> N;                 // number of observations
  int<lower=0> J[3];                 // total number of outcomes
  int<lower=0> K;                 // number of covariates for concatenated design matrices (X1, .., XK)
  real Yn[N, J[1]];                 // normal outcomes
  int<lower=0,upper=1> Yb[N, J[2]]; // bernoulli outcomes
  int<lower=0> Yp[N, J[3]];         // poisson outcomes
  matrix[N, K] X;                  // concatenated design matrices (X1, ..., XK)
  int<lower=0> Kj[3];             // J-dim integer array giving how many covariates per outcome
  int<lower=0> special;
}
transformed data {
  int J_all = sum(J);
  matrix[J_all, N] bounds_ind[2] = get_bounds_indicator(Yn, Yb, Yp, N, J[1], J[2], J[3], J_all);
  vector[J_all] mu_zero = rep_vector(0, J_all);
}
parameters {
  vector[K] beta;                             // long vector of all regression coefficients
  cholesky_factor_corr[J_all] L;                  // Cholesky decomposition of JxJ correlation matrix
  real<lower=0> sigmasq[J[1]];                  // Jn-dim array of variances (may be 0-dim)
  vector<lower=0,upper=1>[N] uraw_b[J[2]];         // latent variables for bernoulli outcomes
  vector<lower=0,upper=1>[N] uraw_p[J[3]];         // latent variables for Poisson outcomes
}
model {
  // initialize variables
  real sigma[J[1]];       // stdev
  matrix[N, J_all] mu_glm = get_means(beta, X, J, Kj);

  // get standard deviations and declare prior on sigmasq if applicable
  if ( J[1] > 0 ) {
    sigma   = sqrt(sigmasq);
    sigmasq ~ inv_gamma(1.0e-4, 1.0e-4);

    for (j in 1:J[1]) to_vector(Yn[, j]) ~ normal( to_vector(mu_glm[, j]), sigma[j] );
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

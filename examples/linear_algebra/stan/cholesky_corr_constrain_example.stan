functions {
  #include corr_constrain_lp.stan
  #include cholesky_corr_constrain_lp.stan
}
 data {
  int K;
  real eta;
}
transformed data {
  int k_choose_2 = (K * (K - 1)) / 2;
}
parameters {
  vector<lower=0>[k_choose_2] y_raw;
}
transformed parameters {
  matrix[K, K] y = cholesky_corr_constrain_lp(y_raw, K);
}
model {
  y ~ lkj_corr_cholesky(eta);
}
generated quantities {
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(y);
}

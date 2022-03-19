functions {
  // #include normal_copula.stanfunctions
}
data {
  int<lower=0> N;
  int K;
  matrix[N, K] x;
}
parameters {
  array[K - 1] real mu;
  array[K + 1] real<lower=0> sigma;
  cholesky_factor_corr[K] L;
}
model {
  target += normal_lpdf(x[ : , 1] | mu[1], sigma[1]);
  target += gumbel_lpdf(x[ : , 2] | mu[2], sigma[2]);
  target += lognormal_lpdf(x[ : , 3] | mu[3], sigma[3]);
  target += weibull_lpdf(x[ : , 4] | sigma[4], sigma[5]);
  
  {
    matrix[K, N] y;
    for (n in 1 : N) {
      y[1, n] = inv_Phi(normal_cdf(x[n, 1] | mu[1], sigma[1]));
      y[2, n] = inv_Phi(gumbel_cdf(x[n, 2] | mu[2], sigma[2]));
      y[3, n] = inv_Phi(lognormal_cdf(x[n, 3] | mu[3], sigma[3]));
      y[4, n] = inv_Phi(weibull_cdf(x[n, 4] | sigma[4], sigma[5]));
    }
     y ~ multi_normal_cholesky_copula(L);
  }
}
generated quantities {
  matrix[K, K] Sigma = multiply_lower_tri_self_transpose(L);
}

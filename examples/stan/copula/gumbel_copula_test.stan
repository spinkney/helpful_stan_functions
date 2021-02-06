functions {
  #include gumbel_copula.stan
}
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  int is_vector;
}
parameters {
  real mu[2];
  real<lower=0> sigma[2];
  real<lower=0, upper=1> tau;
}
transformed parameters {
  real<lower=1> theta = 1 / (1 - tau);
}
model {
  target += lognormal_lpdf(x | mu[1], sigma[1]);
  target += lognormal_lpdf(y | mu[2], sigma[2]);
  
  if(is_vector == 0){
    for (n in 1:N)
    target += gumbel_copula(lognormal_cdf(x[n], mu[1], sigma[1]),
                            lognormal_cdf(y[n], mu[2], sigma[2]), theta);
  } else {
      target += gumbel_copula_vector(lognormal_cdf(x, mu[1], sigma[1]),
                                     lognormal_cdf(y, mu[2], sigma[2]), theta);
  }

}

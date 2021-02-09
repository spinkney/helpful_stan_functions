functions {
  #include clayton_copula.stan
}
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real mu[2];
  real<lower=0> sigma[2];
  real<lower=0> theta;
}
model {
  target += lognormal_lpdf(x | mu[1], sigma[1]);
  target += lognormal_lpdf(y | mu[2], sigma[2]);
  
  for (n in 1:N)
    target += clayton_copula(lognormal_cdf(x[n], mu[1], sigma[1]),
                            lognormal_cdf(y[n], mu[2], sigma[2]), theta);
}

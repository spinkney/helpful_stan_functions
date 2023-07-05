functions {
  #include frank_copula.stanfunctions
}
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real theta;
  real mu;
  real<lower=0> sigma;
  real<lower=0> alpha;
}
model {
  vector[N] u;
  vector[N] v;
  
  target += lognormal_lpdf(x | mu, sigma);
  target += exponential_lpdf(y | alpha);
  
  for (n in 1:N) {
    u[n] = lognormal_cdf(x[n] | mu, sigma);
    v[n] = exponential_cdf(y[n] | alpha);
  }
  
  target += frank_copula_lpdf(u | v, theta);
}

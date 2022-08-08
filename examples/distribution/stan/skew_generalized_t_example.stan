functions {
  #include skew_generalized_t.stanfunctions
}
data {
  int N;
  vector[N] x;
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=-1, upper=1> lambda;
  real<lower=0> q;
  real<lower=2.0/q> p;
}
model {
  mu ~ normal(0, 4);
  sigma ~ exponential(1);
  lambda ~ std_normal();
  p ~ normal(2, 1);
  q ~ normal(7, 2);
  x ~ skew_generalized_t(mu, sigma, lambda, p, q);
}

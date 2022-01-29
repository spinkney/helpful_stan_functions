functions {
  #include unit_johnson_su.stanfunctions
}
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  y ~ unit_johnson(mu, sigma);
}

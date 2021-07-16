functions {
 #include multi_wallenius_hypergeometric.stan
}
data {
  int<lower=0> N;
  int<lower=0> C;
  int y[N, C + 1];
  vector[C] m;
  real tol;
}
transformed data {
  real x_r[0];
}
parameters {
  simplex[C] probs;
}
model {
  for (i in 1:N) 
    y[i] ~ multi_wallenius(m, probs, x_r, tol);
}

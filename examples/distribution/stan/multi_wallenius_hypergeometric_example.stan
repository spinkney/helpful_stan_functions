functions {
  #include multi_wallenius_hypergeometric.stanfunctions
}
data {
  int<lower=0> N;
  int<lower=0> C;
  array[N, C + 1] int y;
  vector[C] m;
  real tol;
}
transformed data {
  array[0] real x_r;
}
parameters {
  simplex[C] probs;
}
model {
  for (i in 1:N) 
    y[i] ~ multi_wallenius(m, probs, x_r, tol);
}

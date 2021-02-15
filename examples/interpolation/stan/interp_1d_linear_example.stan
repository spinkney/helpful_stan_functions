functions {
#include interp_1d_linear.stan
}

data {
  int<lower=2> N_in;
  int<lower=1> N_out;
  int<lower=1> D;
  vector[D] y[N_in];
  real x[N_in];
  real x_out[N_out];
}

parameters {
  real theta; // dummy parameter, because at least one param needed to
  // sample model with cmdstanr
}

model {
  theta ~ std_normal();
}
generated quantities {
  vector[D] y_out[N_out] = interp_1d_linear(y, x, x_out);
}

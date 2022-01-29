functions {
  #include interp_1d_linear.stanfunctions
}
data {
  int<lower=2> N_in;
  int<lower=1> N_out;
  int<lower=1> D;
  array[N_in] vector[D] y;
  array[N_in] real x;
  array[N_out] real x_out;
}
parameters {
  real theta; // dummy parameter, because at least one param needed to
  // sample model with cmdstanr
}
model {
  theta ~ std_normal();
}
generated quantities {
  array[N_out] vector[D] y_out = interp_1d_linear(y, x, x_out);
}

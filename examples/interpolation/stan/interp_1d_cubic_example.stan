functions {
  vector derivative_fun(real t, vector y, data array[] int a0, vector theta) {
    return (cos(theta[1] * y));
  }
  #include interp_1d_cubic.stanfunctions
}
data {
  int<lower=2> N_in;
  int<lower=1> N_out;
  int<lower=1> D;
  array[N_in] vector[D] y;
  array[N_out] vector[D] y_data;
  array[N_in] real x;
  array[N_out] real x_out;
}
transformed data {
  array[0] int a0;
}
parameters {
  vector<lower=0>[1] theta;
}
transformed parameters {
  array[N_out] vector[D] y_out = interp_1d_cubic(y, x, x_out, a0, theta);
}
model {
  theta ~ std_normal();
  for (n in 1 : N_out) {
    y_data[n] ~ normal(y_out[n], 0.2);
  }
}

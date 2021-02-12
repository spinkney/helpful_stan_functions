functions {
  vector derivative_fun(real t, vector y, vector theta, 
      data real[] x_r, data int[] x_i) {
    return(cos(theta[1]*y));
  }
#include interp_1d_cubic.stan
}

data {
  int<lower=2> N_in;
  int<lower=1> N_out;
  int<lower=1> D;
  vector[D] y[N_in];
  vector[D] y_data[N_out];
  real x[N_in];
  real x_out[N_out];
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  vector<lower=0>[1] theta;
}

transformed parameters {
  vector[D] y_out[N_out] = interp_1d_cubic(y, x, x_out, theta, x_r, x_i);
}

model {
  theta ~ std_normal();
  for(n in 1:N_out) {
    y_data[n] ~ normal(y_out[n], 0.2);
  }
}

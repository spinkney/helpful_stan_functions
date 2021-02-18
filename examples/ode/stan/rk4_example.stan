functions {
  // Lotka-Volterra system
  vector derivative_fun(real t, vector y, vector theta, 
      data real[] x_r, data int[] x_i) {
    vector[2] dydt;
    dydt[1] = theta[1]*y[1] - y[1]*y[2];
    dydt[2] = y[1]*y[2] - theta[2]*y[2];
    return dydt;
  }
#include interp_1d_cubic.stan
#include odeint_rk4.stan
}

data {
  int<lower=1> N;
  real t_data[N]; // must be increasing
  vector[2] y_data[N];
  vector[2] y0;
  int<lower=1> num_steps;
  real<lower=0> h; // step size
  real t_grid[num_steps + 1];
}

transformed data {
  real x_r[0];
  int x_i[0];
  real t0 = t_grid[1];
}

parameters {
  vector<lower=0>[2] theta;
}

transformed parameters {
  vector[2] y_grid[num_steps + 1] = odeint_rk4(t0, y0, h, num_steps,
      theta, x_r, x_i);
}

model {
  theta ~ normal(1.0, 0.2);
}

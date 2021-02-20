functions {
  // Lotka-Volterra system
  vector derivative_fun(real t, vector y, data int[] a0, vector theta) {
    vector[2] dydt;
    dydt[1] = theta[1]*y[1] - y[1]*y[2];
    dydt[2] = y[1]*y[2] - theta[2]*y[2];
    return dydt;
  }
#include interp_1d_cubic.stan
#include odeint_euler.stan
#include odeint_midpoint.stan
#include odeint_rk4.stan
}

data {
  int<lower=1> N;
  real t_data[N]; // must be increasing
  vector[2] y_data[N];
  vector[2] y0;
  real t0;
  int<lower=1> num_steps; // number of steps
  real<lower=0> h; // step size
}

transformed data {
  int a0[0];
  int G = num_steps + 1;
  real t_grid[G];
  t_grid[1] = t0;
  for (j in 1:num_steps) {
    t_grid[j+1] = t0 + j*h;
  }
}

parameters {
  vector<lower=0>[2] theta;
  real<lower=0> sigma;
}

transformed parameters {
  
  // Solve ODE at a grid of points (using different solvers)
  vector[2] y_grid_rk4[G] = odeint_rk4(t0, y0, h, num_steps, a0, theta);
  vector[2] y_grid_mid[G] = odeint_midpoint(t0, y0, h, num_steps, a0, theta);
  vector[2] y_grid_eul[G] = odeint_euler(t0, y0, h, num_steps, a0, theta);
      
  // Interpolate solution to data time points
  vector[2] y_rk4[N] = interp_1d_cubic(y_grid_rk4, t_grid, t_data, a0, theta);
  vector[2] y_mid[N] = interp_1d_cubic(y_grid_mid, t_grid, t_data, a0, theta);
  vector[2] y_eul[N] = interp_1d_cubic(y_grid_eul, t_grid, t_data, a0, theta);
    
  // Compute solution using built-in BDF solver as reference
  vector[2] y_bdf[N] = ode_bdf(derivative_fun, y0, t0, t_data, a0, theta);
}

model {
  theta ~ normal(1.0, 0.2);
  sigma ~ inv_gamma(5, 5);
  for(j in 1:N){
    y_data[j] ~ normal(y_rk4[j], sigma);
  }
}

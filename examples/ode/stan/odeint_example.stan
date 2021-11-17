functions {
  // Lotka-Volterra system
  vector derivative_fun(real t, vector y, data array[] int a0, vector theta) {
    vector[2] dydt;
    dydt[1] = theta[1] * y[1] - y[1] * y[2];
    dydt[2] = y[1] * y[2] - theta[2] * y[2];
    return dydt;
  }
  #include interp_1d_cubic.stanfunctions
  #include odeint_euler.stanfunctions
  #include odeint_midpoint.stanfunctions
  #include odeint_rk4.stanfunctions
}
data {
  int<lower=1> N;
  array[N] real t_data; // must be increasing
  array[N] vector[2] y_data;
  vector[2] y0;
  real t0;
  int<lower=1> num_steps; // number of steps
  real<lower=0> h; // step size
}
transformed data {
  array[0] int a0;
  int G = num_steps + 1;
  array[G] real t_grid;
  t_grid[1] = t0;
  for (j in 1 : num_steps) {
    t_grid[j + 1] = t0 + j * h;
  }
}
parameters {
  vector<lower=0>[2] theta;
  real<lower=0> sigma;
}
transformed parameters {
  // Solve ODE at a grid of points (using different solvers)
  array[G] vector[2] y_grid_rk4 = odeint_rk4(t0, y0, h, num_steps, a0, theta);
  array[G] vector[2] y_grid_mid = odeint_midpoint(t0, y0, h, num_steps, a0, theta);
  array[G] vector[2] y_grid_eul = odeint_euler(t0, y0, h, num_steps, a0, theta);
  
  // Interpolate solution to data time points
  array[N] vector[2] y_rk4 = interp_1d_cubic(y_grid_rk4, t_grid, t_data, a0, theta);
  array[N] vector[2] y_mid = interp_1d_cubic(y_grid_mid, t_grid, t_data, a0, theta);
  array[N] vector[2] y_eul = interp_1d_cubic(y_grid_eul, t_grid, t_data, a0, theta);
  
  // Compute solution using built-in BDF solver as reference
  array[N] vector[2] y_bdf = ode_bdf(derivative_fun, y0, t0, t_data, a0, theta);
}
model {
  theta ~ normal(1.0, 0.2);
  sigma ~ inv_gamma(5, 5);
  for (j in 1 : N) {
    y_data[j] ~ normal(y_rk4[j], sigma);
  }
}

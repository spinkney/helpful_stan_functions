library(cmdstanr)

# Compile model
fp <- file.path("../stan/rk4_example.stan")
inc1 <- file.path("../../../functions/interpolation/")
inc2 <- file.path("../../../functions/ode/")
inc_paths <- c(inc1, inc2)
model <- cmdstan_model(fp, include_paths = inc_paths, force_recompile = T)

# Create test data
t_data <- seq(0.2, 20, by = 0.2)
N <- length(t_data)
y_data <- array(0, c(N, 2))
y0 <- c(1,2)
t_grid <- c(0, t_data)
stan_data <- list(
  N = N, t_data = t_data, y_data = y_data, y0 = y0, t_grid = t_grid,
  h = 0.2, num_steps = length(t_data)
)

# Run model
fit <- model$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
)

# Plot
N_grid <- stan_data$num_steps + 1
y_out <- fit$draws(variables = "y_grid")
plot(t_data, y_data[,1], "n", xlim = c(0, 20), ylim = c(0, 4))
for (j in 1:3) {
  yj <- as.vector(y_out[j,1,])
  y1 <- yj[1:N_grid]
  y2 <- yj[(N_grid + 1):(2*N_grid)]
  lines(t_grid, y1, col = "firebrick3")
  lines(t_grid, y2, col = "steelblue")
}

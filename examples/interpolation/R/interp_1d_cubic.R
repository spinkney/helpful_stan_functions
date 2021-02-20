library(cmdstanr)

# Compile model
fp <- file.path("../stan/interp_1d_cubic_example.stan")
inc_paths <- file.path("../../../functions/interpolation/")
model <- cmdstan_model(fp, include_paths = inc_paths, force_recompile = T)

# Create test data
x <- c(0, 8, 10)
y <- sin(0.8*x)
y <- cbind(y)
x_out <- seq(0.2, 10, by = 0.23)
y_data <- sin(0.8*x_out) + 0.2*rnorm(n = length(x_out))
y_data <- cbind(y_data)
N_out <- length(x_out)
stan_data <- list(
  x = x, y = y, x_out = x_out, D = ncol(y), y_data = y_data,
  N_in = length(x), N_out = N_out
)

# Run model
fit <- model$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
)

# Plot data (black dots) and interpolation (blue line) and points
# where to interpolate from (red dots)
par(mfrow = c(1,1))
plot(x_out, y_data, "n", xlim = c(0, 10))
points(x, y, pch = 16, col = "red")
points(x_out, y_data, pch = 20)
for (j in 1:10) {
  y_out <- as.vector(fit$draws(variables = "y_out")[j,1,])
  lines(x_out, y_out, col = "steelblue")
}

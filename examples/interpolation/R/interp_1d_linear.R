library(cmdstanr)

# Compile model
fp <- file.path("../stan/interp_1d_linear_example.stan")
inc_paths <- file.path("../../../functions/interpolation/")
model <- cmdstan_model(fp, include_paths = inc_paths, force_recompile = T)

# Create test data
x <- c(3, 5.5, 6.9, 7.1, 10)
y1 <- c(4, 5, -1, 1.5, 0)
y2 <- sin(x)
y <- cbind(y1, y2)
x_out <- seq(0.1, 10, by = 0.6)
N_out <- length(x_out)
stan_data <- list(
  x = x, y = y, x_out = x_out, D = ncol(y),
  N_in = length(x), N_out = N_out
)

# Run model
fit <- model$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 5,
  iter_sampling = 5,
)

# Extract and plot interpolation result (black = original, red = interpolation)
y_out <- fit$draws(variables = "y_out")[1,1,]
y1_out <- y_out[1:N_out]
y2_out <- y_out[(N_out + 1):(2*N_out)]

par(mfrow = c(2,1))
plot(x, y1, 'o', pch  = 16)
points(x_out, y1_out, col = "red", pch = 20)

plot(x, y2, 'o', pch  = 16)
points(x_out, y2_out, col = "red", pch = 20)

library(cmdstanr)

# Compile model
fp <- file.path("../stan/odeint_example.stan")
inc1 <- file.path("../../../functions/interpolation/")
inc2 <- file.path("../../../functions/ode/")
inc_paths <- c(inc1, inc2)
model <- cmdstan_model(fp, include_paths = inc_paths, force_recompile = T)

# Create test data
t_data <- seq(0.1, 15, by = 0.1)
N <- length(t_data)
y_data <- array(0, c(N, 2))
y0 <- c(1, 2)
t_grid <- c(0, t_data)
stan_data <- list(
  N = N, t_data = t_data, y_data = y_data, y0 = y0, t_grid = t_grid,
  h = 0.1, num_steps = length(t_data)
)

# Run model
fit <- model$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
)

get_output <- function(fit, name, chain_idx, draw_idx) {
  y_out <- fit$draws(variables = name)
  yj <- as.matrix(y_out[draw_idx, chain_idx, ])
  N_grid <- length(yj) / 2
  y1 <- yj[1:N_grid]
  y2 <- yj[(N_grid + 1):(2 * N_grid)]
  return(cbind(y1, y2))
}

# Plot one draw
j <- 100
y_rk4 <- get_output(fit, "y_rk4", 1, j)
y_eul <- get_output(fit, "y_eul", 1, j)
y_mid <- get_output(fit, "y_mdp", 1, j)

par(mfrow = c(2, 1))
cols <- c("steelblue", "orange", "firebrick3")
for (j in 1:2) {
  plot(t_data, y_data[, j], "n",
    xlim = c(0, 15), ylim = c(0, 4),
    xlab = "t", ylab = paste0("y", j)
  )
  lines(t_grid, y_eul[, j], col = cols[1] , lwd = 2)
  lines(t_grid, y_mid[, j], col = cols[2], lwd = 2)
  lines(t_grid, y_rk4[, j], col = cols[3], lty = 2, lwd = 2)
}
legend(
  x = 0.5, y =4, lty = c(1, 1, 2), lwd = 2, col = cols,
  legend = c("euler", "midpoint", "rk4")
)

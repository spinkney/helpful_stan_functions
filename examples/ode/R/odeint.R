library(cmdstanr)

# Compile model
fp <- file.path("../stan/odeint_example.stan")
inc1 <- file.path("../../../functions/interpolation/")
inc2 <- file.path("../../../functions/ode/")
inc_paths <- c(inc1, inc2)
model <- cmdstan_model(fp, include_paths = inc_paths, force_recompile = T)

# Load test data
dat <- readRDS("test_data.rds")
P <- 80
t_data <- dat$t_data[1:P]
y_data <- dat$y_data[1:P,]
N <- length(t_data)

# Set initial state
y0 <- c(1,2)
t0 <- 0
h <- 0.37
num_steps <- ceiling(max(t_data)/h)
stan_data <- list(
  N = N, t_data = t_data, y_data = y_data, y0 = y0, t0 = t0,
  h = h, num_steps = num_steps
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
t_grid <- seq(0, num_steps*h, by = h)
j <- 1
y_grid_rk4 <- get_output(fit, "y_grid_rk4", 1, j)
y_grid_eul <- get_output(fit, "y_grid_eul", 1, j)
y_grid_mid <- get_output(fit, "y_grid_mid", 1, j)
y_rk4 <- get_output(fit, "y_rk4", 1, j)
y_eul <- get_output(fit, "y_eul", 1, j)
y_mid <- get_output(fit, "y_mid", 1, j)
y_bdf <- get_output(fit, "y_bdf", 1, j)

par(mfrow = c(2, 1))
par(mar = c(2.25,4,1,1))
cols <- c("black", "steelblue", "orange", "firebrick3")
for (j in 1:2) {
  plot(t_data, y_data[, j], col = "gray", pch = 20,
    xlim = c(0, max(t_data)), ylim = c(0, 3.6),
    xlab = "t", ylab = paste0("y", j)
  )
  points(t_grid, y_grid_eul[, j], col = cols[2], pch = 20)
  points(t_grid, y_grid_mid[, j], col = cols[3], pch = 20)
  points(t_grid, y_grid_rk4[, j], col = cols[4], pch = 20)
  lines(t_data, y_bdf[, j], col = cols[1])
  lines(t_data, y_eul[, j], col = cols[2])
  lines(t_data, y_mid[, j], col = cols[3])
  lines(t_data, y_rk4[, j], col = cols[4], lty = 2)
  if (j == 1) {
    legend(
      x = 0.5, y = 3.5, lty = c(1, 1, 1, 2), lwd = 2, col = cols,
      legend = c("bdf", "euler", "midpoint", "rk4")
    )
  }
}


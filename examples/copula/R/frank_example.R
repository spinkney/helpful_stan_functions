library(cmdstanr)
library(copula)

fp <- file.path("./examples/copula/stan/frank_copula_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/copula")

# Set the copula parameter
theta <- -10

# Create the Frank copula object
copula <- frankCopula(param = theta)

# Generate random variables from the copula
n <- 200 # Number of data points
u <- rCopula(n, copula)

# Apply inverse transformation to get simulated values for x and y
x <- qlnorm(u[, 1])
y <- qexp(u[, 2])

mod_out <- mod$sample(
  data = list(N = n,
              x = x,
              y = y),
  parallel_chains = 4
)

library(cmdstanr)
fp <- file.path("examples/inverse_cdf/stan/lognormal_qf_example.stan")
mod <- cmdstan_model(fp,  
                     include_paths = "./functions/qf")

log_norm_loc <- function (mu, sigma) {
  log(mu^2 / sqrt(mu^2 + sigma^2))
}

log_norm_scale <- function (mu, sigma) {
  log( 1 + ( sigma^2 / mu^2 ) )
}

mu <- 15
sigma <- 10

mod_trunc <- mod$sample(
  data = list(N = 500,
              mu = log_norm_loc(mu, sigma),
              sigma = log_norm_scale(mu, sigma),
              lb = 0,
              ub = 20),
  chains = 2,
  parallel_chains = 2,
  fixed_param = TRUE
)

hist(mod_trunc$draws("y_out"))

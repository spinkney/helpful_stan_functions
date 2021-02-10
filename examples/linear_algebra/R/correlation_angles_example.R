library(cmdstanr)
fp <- file.path("./examples/stan/linear_algebra/to_psd.stan")
mod <- cmdstan_model(fp) #, include_paths = "./functions/copula")

T <- matrix(c(0.8, -0.9, -0.9,
              -0.9, 1.1, 0.3,
              -0.9, 0.4, 0.9), 3, 3)


mod_out <- mod$optimize(
  data = list(N = 3,
              R = T)
)

mod_out <- mod$sample(
  data = list(N = 3,
              R = T),
  chains = 2,
  adapt_delta = 0.9,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200
)

mat <- matrix(mod_out$summary("R_hat")$mean, 3, 3)
mat
chol(mat)
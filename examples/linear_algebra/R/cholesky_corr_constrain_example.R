library(cmdstanr)

fp <- file.path("./examples/linear_algebra/stan/cholesky_corr_constrain_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/linear_algebra")

fp <- file.path("../positive_corr.stan")
mod <- cmdstan_model(fp, force_recompile = T)

K <- 10
mod_out <- mod$sample(
  data = list(K = K,
              eta = 100),
  chains = 2,
  adapt_delta = 0.8,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200
)

matrix(mod_out$summary("Omega")$mean, K, K)

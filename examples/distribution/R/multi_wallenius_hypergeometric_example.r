### multi-wallenius
fp <- file.path("./examples/distribution/stan/multi_wallenius_hypergeometric_example.stan")
mod <- cmdstan_model(fp, force_recompile = T)

N <- 10
m <- c(500, 1000, 800)
n <- 55

odds <- c(0.2, 0.7, 0.1)
y <- rMWNCHypergeo(N, m, n, odds)


meanMWNCHypergeo(m, n, odds, precision=1E-7)


mod_out <- mod$sample(
  data = list(N = N,
              C = length(m),
              y = cbind(n, t(y)),
              m = m,
              tol = 0.01),
  chains = 2,
  init = 1,
  adapt_delta = 0.8,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200
)

mod_out

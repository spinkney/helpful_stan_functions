library(sgt)
library(cmdstanr)

# SINGLE VARIABLE ESTIMATION:
### generate random variable
set.seed(7900)
n <- 1000
x <- rsgt(n, mu = 2, sigma = 2, lambda = -0.25, p = 1.7, q = 7)
### Get starting values and estimate the parameter values
start <- list(mu = 0, sigma = 1, lambda = 0, p = 2, q = 10)
result <- sgt.mle(X.f = ~ x, start = start, method = "nlminb")
print(result)
print(summary(result))

fp <- file.path("examples/distribution/stan/skew_generalized_t_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/distribution")

mod_skew_gen_t <- mod$sample(
  data = list(N = n,
              x = x),
  chains = 2,
  parallel_chains = 2
)
mod_skew_gen_t$summary(c("mu", "sigma", "lambda", "p", "q"))

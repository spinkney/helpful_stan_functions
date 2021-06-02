library(cmdstanr)

fp <- file.path("./examples/distribution/stan/unit_johnson_su_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/distribution", 
                     force_recompile = T)

y <- c(0.04, 0.02, 0.06, 0.12, 0.14, 0.08, 0.22, 0.12, 0.08, 0.26, 0.24, 0.04, 0.14, 0.16, 0.08, 0.26, 0.32, 0.28, 0.14, 0.16, 0.24,
       0.22, 0.12, 0.18, 0.24, 0.32, 0.16, 0.14, 0.08, 0.16, 0.24, 0.16, 0.32, 0.18, 0.24, 0.22, 0.16, 0.12, 0.24, 0.06, 0.02, 0.18,
       0.22, 0.14, 0.06, 0.04, 0.14, 0.26, 0.18, 0.16)

mod$sample(
  data = list(
    N = length(y),
    y = y
  ),
  seed = 12312,
  chains = 2,
  adapt_delta = 0.8,
  parallel_chains = 2
)
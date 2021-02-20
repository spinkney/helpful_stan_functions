library(cmdstanr)
fp <- file.path("./test/copula/stan/bivariate_normal_copula_cdf.stan")
mod <- cmdstan_model(fp, include_paths = c("./functions/unit_test",
                                           "./functions/copula"))

my_samples <- mod$sample(chains = 1,
                         iter_sampling = 1,
                         fixed_param = T)

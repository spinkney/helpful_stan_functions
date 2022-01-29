library(cmdstanr)
library(copula)
library(ExtDist)
library(VGAM)

fp <- file.path("./examples/copula/stan/normal_copula_example_test.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/copula", 
                     force_recompile = T)

G3 <- normalCopula(-.3, dim=2)
gMvd2 <- mvdc(G3, c("norm","Gumbel"), param = 
                list(list(mean = 3.5, sd = 0.5), 
                     list(location = 0, scale = 1)))
set.seed(11)

N <- 500
x <- rMvdc(N, gMvd2)

mod_out <- mod$sample(
  data = list(N = N,
              x = x[, 1],
              y = x[, 2],
              is_vector = 0),
  chains = 2,
  adapt_delta = 0.8,
  parallel_chains = 2,
  iter_warmup = 1000,
  iter_sampling = 1000
)



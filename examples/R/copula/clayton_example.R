library(cmdstanr)
library(copula)

fp <- file.path("./examples/stan/copula/clayton_copula_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/copula")

clayton_mod <- claytonCopula(1.37, dim=2)
clayton_out <-  mvdc(clayton_mod, c("lnorm","lnorm"), param = 
                list(list(meanlog = 3.5, sdlog = 0.5), 
                     list(meanlog = 2.75, sdlog = 0.75)))
set.seed(11)

N <- 100
x <- rMvdc(N, clayton_out)

mod_out <- mod$sample(
  data = list(N = N,
              x = x[, 1],
              y = x[, 2]),
  chains = 2,
  adapt_delta = 0.8,
  parallel_chains = 2,
  iter_warmup = 1000,
  iter_sampling = 1000
)

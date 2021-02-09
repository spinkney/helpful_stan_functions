# library(rstan)
# expose_stan_functions("gumbel_copula.stan")
# gumbel_copula_pdf <- function(uv, theta) 
#   exp(gumbel_copula(uv[1], uv[2], theta))
# 
# library(cubature)
# adaptIntegrate(gumbel_copula_pdf, lowerLimit = c(0,0), upperLimit = c(1,1), 
#                theta = sqrt(2), doChecking = TRUE)


library(cmdstanr)
library(copula)

fp <- file.path("./examples/stan/copula/gumbel_copula_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/copula")

G3 <- gumbelCopula(1.37, dim=2)
gMvd2 <- mvdc(G3, c("lnorm","lnorm"), param = 
                list(list(meanlog = 3.5, sdlog = 0.5), 
                     list(meanlog = 2.75, sdlog = 0.75)))
set.seed(11)

N <- 100
x <- rMvdc(N, gMvd2)

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

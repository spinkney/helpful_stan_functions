library(cmdstanr)
library(copula)
library(ExtDist)
library(VGAM)

fp <- file.path("./examples/stan/copula/multi_normal_copula_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/copula", 
                     force_recompile = T)

G3 <- normalCopula(param=c(0.4,0.2,-0.8, -0.3, 0.15, -0.5), dim = 4, dispstr = "un")
gMvd2 <- mvdc(G3, c("norm","Gumbel", "lnorm", "weibull"), param = 
                list(list(mean = 3.5, sd = 0.5), 
                     list(location = 0, scale = 1),
                     list(meanlog = 2.75, sdlog = 0.75),
                     list(shape = 1.0, scale = 1.5)))
set.seed(11)

N <- 500
x <- rMvdc(N, gMvd2)

mod_out <- mod$sample(
  data = list(N = N,
              K = 4,
              x = x),
  chains = 2,
  adapt_delta = 0.8,
  parallel_chains = 2,
  iter_warmup = 1000,
  iter_sampling = 1000
)

mod_out$summary("Sigma")
p2P(c(0.4,0.2,-0.8, -0.3, 0.15, -0.5))

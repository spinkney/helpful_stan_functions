library(cmdstanr)
# setwd("~\helpful_stan_functions")
fp <- file.path("./examples/copula/stan/centered_gaussian_copula_example.stan")
mod <- cmdstan_model(fp, include_paths = c("./functions",
                                           "./functions/copula",
                                           "./functions/distribution",
                                           "../functions"))

## obtain non-extreme random correlation matrix
Gamma = matrix(1, 3, 3)
while( any( abs(Gamma[lower.tri(Gamma)]) > 0.6 ) ) {
  Gamma = randcorr::randcorr(3)
}
n     = 500
x     = rnorm(n, mean = 5, sd = 5/3)
X     = cbind(1, x)
beta  = c(1, -.25)
eta   = X %*% beta
mu1   = eta
mu2   = binomial()$linkinv(eta)
mu3   = exp(eta)
Z     = mvtnorm::rmvnorm(n, sigma = Gamma)
U     = pnorm(Z)
yn    = qnorm(U[, 1])
yb    = qbinom(U[, 2], 1, mu2)
yp    = qpois(U[, 3], mu3)
data  = data.frame(yn, yb, yp, x)


formula.list = list(
  yb ~ x,
  yn ~ x,
  yp ~ x
)
family.list = list(
  binomial(),
  gaussian(),
  poisson()
)

#' Joint analysis of multivariate GLMs of mixed types
#'
#' Uses rstan to estimate a multivariate GLM with correlated
#' responses via Gaussian copula
#'
#' @param formula.list a list of formulas
#' @param family.list  a list of families giving the density of each formula. Canonical link functions always used.
#' @param data         data.frame giving all variables in `formula.list`
#' @param ...          arguments to pass onto `rstan::sampling`
#'
#' @return             rstan object fitting the posterior
get_data = function(formula.list, family.list, data) {

  N  = nrow(data)
  J  = length(formula.list)

  ## reorder variables if necessary
  n.indx    = sapply(family.list, function(x) x$family) == 'gaussian'
  b.indx    = sapply(family.list, function(x) x$family) == 'binomial'
  p.indx    = sapply(family.list, function(x) x$family) == 'poisson'
  new.order = c(which(n.indx), which(b.indx), which(p.indx))

  ## gtotal number of endpoints and matrix of each outcome type
  Jn = sum(n.indx)
  Jb = sum(b.indx)
  Jp = sum(p.indx)
  Y  = sapply(formula.list, function(x) data[, all.vars(x)[1]])
  Yn = Y[, n.indx, drop = F]
  Yb = Y[, b.indx, drop = F]
  Yp = Y[, p.indx, drop = F]

  ## store character vector giving dist name and number so results make sense
  distnames = c( paste0('gaussian[', 1:Jn, ']'),
                 paste0('binomial[', 1:Jb, ']'),
                 paste0('poisson[', 1:Jp, ']')
  )


  if ( !all( new.order == 1:J ) ) {
    message("Rearranging order of dependent variables to normal, binomial, poisson")
    formula.list = formula.list[new.order]
    family.list  = family.list[new.order]
  }

  ## get design matrices and number of covariates per regression
  Xlist = lapply(formula.list, model.matrix)
  Xbig  = do.call(cbind, Xlist)
  K     = ncol(Xbig)
  Kj    = sapply(Xlist, ncol)

  ## obtain matrix giving start and end values of x matrix per endpoint
  Xindx  = matrix(nrow = J, ncol = 2)
  xstart = 1
  for ( j in 1:J ) {
    xend = xstart + Kj[j] - 1
    Xindx[j, ] = c(xstart, xend)
    xstart = xend + 1
  }

  ## construct stan data
  standat = list(
    'N'     = N,
    'J'     = J,
    'Jn'    = Jn,
    'Jb'    = Jb,
    'Jp'    = Jp,
    'K'     = K,
    'Yn'    = Yn,
    'Yb'    = Yb,
    'Yp'    = Yp,
    'X'     = Xbig,
    'Xindx' = Xindx,
    'Kj'    = Kj
  )

  ## fit stan model
  # fit           = rstan::sampling(stanmod, data = standat, ...)
  # fit = stanmod$sample(data = standat,
  #                      parallel_chains = 4,
  #                      iter_warmup = 1000,
  #                      iter_sampling = 1000)
  # ## rename betas for clarity about which equation each comes from, e.g., beta[j][0:(Kj-1)]
  # beta.old.indx = which(grepl('beta', names(fit)))
  # beta.new      = character()
  # for( j in 1:J ) {
  #   beta.new = c(beta.new, paste0(distnames[j], '[', seq(0,Kj[j] - 1), ']') )
  # }
  # names(fit)[beta.old.indx] = beta.new
  #
  ## return stan object
  return(standat)
}



stan_data <- get_data(formula.list[1:3], family.list[1:3], data)
stan_data[["J"]] <- rep(1, 3)

# run with "optimizations" which turn out to be slower for this case
stan_data[["special"]] <- 1
# or not
# stan_data[["special"]] <- 0
mod_out <- mod$sample(data = stan_data,
                                    seed = 908453,
                                    init = 0 ,
                                    chains = 2,
                                    parallel_chains = 2,
                                    iter_warmup = 400,
                                    iter_sampling = 400)


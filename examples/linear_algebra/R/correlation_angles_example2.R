library(cmdstanr)
fp <- file.path("./examples/linear_algebra/stan/correlation_angles_constrain_example.stan")
mod <- cmdstan_model(fp, include_paths = "./functions/linear_algebra", force_recompile = T)

library(rstan)
expose_stan_functions(fp)

library(StanHeaders)
stanFunction("multiply_lower_tri_self_transpose", x = L)


N <- 3
test_mat <- matrix(c(0.8, -0.9, -0.9,
              -0.9, 1.1, 0.3,
              -0.9, 0.4, 0.9), 3, 3)

T <- matrix(c(1, -0.9, -0.9,
              -0.9, 1, 0.3,
              -0.9, 0.4, 1), 3, 3)


test_mat <- matrix(c(1.0000,   0,         0,         0,   -0.9360,
                     0,    1.0000,   -0.5500,   -0.3645,   -0.5300,
                     0,   -0.5500,    1.0000,   0,    0.0875,
                     0,   -0.3645,    0,    1.0000,    0.4557,
                     -0.9360,   -0.5300,    0.0875,    0.4557,    1.0000),
                   byrow = T, 5, 5)

where_zero <- (test_mat[upper.tri(test_mat)] == 0) * 1
sum(where_zero) * 1

fit <- stan(file = fp, 
            data = list(N = nrow(test_mat),
                        R = test_mat,
                        is_symmetric = 1,
                        Z = sum(where_zero),
                        K = ( nrow(test_mat) * (nrow(test_mat) - 1 )) / 2,
                        where_zero = where_zero ),
            chains = 1)

K <- 5
theta_test <- build_angle_mat(where_zero, runif(6, 0, pi), K)
L = zero_constrain(theta_test, K)
tcrossprod(L)

angle_raw <-  runif(6)
N <- length(where_zero);
mat <- matrix(0, N, N)
count <- 1;
raw_count <- 1;

for (k in 2:K){
  mat[k, k] = 0;
  for (j in 1:k - 1) {
    if (where_zero[raw_count] != 1) {
      mat[k, j] = angle_raw[count];
      count <= count +  1;
    }
    raw_count = raw_count + 1;
  }
}

count <- 0

for (k in 2:K){
  for (j in 1:(k - 1)) {
    count <- count + 1
    print(count)
  }}



mod_out <- mod$sample(
  data = list(N = nrow(test_mat),
              R = test_mat,
              is_symmetric = 1,
              Z = sum(where_zero),
              K = ( nrow(test_mat) * (nrow(test_mat) - 1 )) / 2,
              where_zero = where_zero ),
  chains = 2,
  #init = 0,
  #seed = 23421,
  #adapt_delta = 0.8,
  #max_treedepth = 10,
  parallel_chains = 4,
  iter_warmup = 400,
  iter_sampling = 400
)

N <- nrow(test_mat)

chol(test_mat)

mod_out <- mod$optimize(
  data = list(N = nrow(test_mat),
                      R = test_mat,
                      is_symmetric = 1,
                      Z = sum(where_zero),
                      K = ( nrow(test_mat) * (nrow(test_mat) - 1 )) / 2,
                      where_zero = where_zero ))

round(matrix(mod_out$summary("R_out")$estimate,nrow(test_mat), nrow(test_mat)), 3)

mod_out <- mod$sample(
  data = list(N = nrow(test_mat),
              R = test_mat,
              is_symmetric = 1),
  chains = 2,
  seed = 23421,
  adapt_delta = 0.9,
  max_treedepth = 10,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200
)
mod_out$summary("R_out")
round(matrix(mod_out$summary("R_out")$mean, N, N), 3)

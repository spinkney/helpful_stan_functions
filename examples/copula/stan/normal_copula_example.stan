functions {
  #include normal_copula.stanfunctions
}
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  int is_vector;
}
parameters {
  array[2] real mu;
  array[2] real<lower=0> sigma;
  real<lower=-1, upper=1> rho;
}
model {
  target += normal_lpdf(x | mu[1], sigma[1]);
  target += gumbel_lpdf(y | mu[2], sigma[2]);
  
  if(is_vector == 0){
    for (n in 1:N)
         target += normal_copula(inv_Phi(normal_cdf(x[n] | mu[1], sigma[1])),
                                 inv_Phi(gumbel_cdf(y[n] | mu[2], sigma[2])), rho);
  
  } else {
     vector[N] x_p;
     vector[N] y_p;
  
      for (n in 1:N){
        x_p[n] = inv_Phi(normal_cdf(x[n] | mu[1], sigma[1]));
        y_p[n] = inv_Phi(gumbel_cdf(y[n] | mu[2], sigma[2]));
      }
  
    target += normal_copula_vector(x_p, y_p, rho);
  }
}

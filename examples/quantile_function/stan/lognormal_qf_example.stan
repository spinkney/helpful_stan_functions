functions {
  #include lognormal_qf.stanfunctions
  
  real lognormal_trunc_rng(real mu, real sigma, real lb, real ub){
    real p_ub = lognormal_cdf(ub, mu, sigma);
    real p_lb = lognormal_cdf(lb, mu, sigma);
    real p = uniform_rng(p_lb, p_ub);
    
    return lognormal_qf(p, mu, sigma);
  }
}
data {
  int N;
  real mu;
  real sigma;
  real<lower=0> lb;
  real<lower=lb> ub;
}
generated quantities {
  array[N] real y_out;

  for(n in 1:N) 
    y_out[n] = lognormal_trunc_rng(mu, sigma, lb, ub);
}

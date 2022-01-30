functions {
 #include correlation_angles.stan
 #include triangular.stan
}
data {
  int N;
  matrix[N, N] R;
  int<lower=0, upper=1> is_symmetric;
}
transformed data {
  int K = N * (N - 1) / 2;
  row_vector[N] m;
  
  for (n in 1:N)
    m[n] = N - n + 1;
}
parameters {
  vector<lower = 0, upper = pi()>[K - (N - 2)] theta_free;
} 
transformed parameters {
  vector[K] theta = append_row(rep_vector(0.5 * pi(), N - 2), theta_free);
}
model {
   matrix[N, N] L = angle2chol(angle_vec2angle_mat(theta, N));
   matrix[N, N] R_hat = multiply_lower_tri_self_transpose(L);
   
  // log_absd_jacs 
  // sin(theta) is always > 0 since theta in (0, pi)
  // because of this the diagonal of L is always > 0 also
   target += 2 * sum(log(sin(theta)));            // angle2chol
   target += N * log(2) + m * diagonal(log(L));   // cholesky tcrossprod
   
   if (is_symmetric == 1)
    lower_elements(R, K) ~ normal(lower_elements(R_hat, K), 0.001);
    else to_vector(R) ~ normal(to_vector(R_hat), 0.001);
}
generated quantities {
   matrix[N, N] R_out = multiply_lower_tri_self_transpose(angle2chol(angle_vec2angle_mat(theta, N)));
}
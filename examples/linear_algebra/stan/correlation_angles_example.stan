functions {
 #include correlation_angles.stan
 #include triangular.stan
}
data {
  int N;
  matrix[N, N] R;
}
transformed data {
  int K = N * (N - 1) / 2;
}
parameters {
  vector[K] theta;
  real<lower=0> sigma;
}
transformed parameters {
  matrix[N, N] R_hat = multiply_lower_tri_self_transpose(angle2chol(angle_vec2angle_mat(theta, N)));
}
model {
  theta ~ std_normal();
  to_vector(R) ~ normal(to_vector(R_hat), sigma);
}

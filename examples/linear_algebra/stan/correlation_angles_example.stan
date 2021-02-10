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
  vector<lower = -pi(), upper = pi()>[K] theta;
}
transformed parameters {
  matrix[N, N] R_hat = multiply_lower_tri_self_transpose(angle2chol(angle_vec2angle_mat(theta, N)));
}
model {
  vector[N] R_vec = to_vector(R);
  vector[N] R_hat_vec = to_vector(R_hat);
  squared_distance(R_vec, R_hat_vec) ~ chi_square(1);
  target += sum(log(fabs(2 * (R_vec - R_hat_vec))));
}

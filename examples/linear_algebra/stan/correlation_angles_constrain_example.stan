functions {
 #include correlation_angles.stan
 #include triangular.stan
 matrix build_angle_mat (vector where_zero, vector angle_raw, int K) {
  int N = num_elements(angle);
  matrix[K, K] mat;
  int count = 1;
  int raw_count = 1;

  for (k in 2:K){
    mat[k, k] = 0;
    for (j in 1:k - 1) {
      raw_count += 1;
      if (where_zero[raw_count] != 1) {
         mat[k, j] = angle_raw[count];
         count += 1;
      }
    }
  }
    
  return mat;
}
 
 vector zero_constrain (matrix angle_raw, int N){
  matrix[N, N] inv_mat;
  matrix[N, N] angle = angle_raw;
  
  inv_mat[1, 1] = 1;
  
  for (i in 2:N) {
    // constrain first column
    // if C = BB^T then for the first column
    // c_{i, 1} = b_{i, 1} 
    // if c_{i, 1} = 0 then cos(\theta_{i, 1}) = 0
    // acos(0) = \theta_{i, 1} = pi / 2
     if (is_nan(angle_raw[i, 1]) == 1) {
       inv_mat[i, 1] = 0;
       angle[i, 1] = pi() / 2;
     } else inv_mat[i, 1] = cos(angle[i, 1]);
    
    inv_mat[i, i] = prod(sin(angle[i, 1:(i - 1)]));
   
    if (i > 2) {
      for (j in 2:(i - 1)) {
        real prod_sines = prod(sin(angle_mat[i, 1:j - 1]));
        if (is_nan(angle_raw[i, j]) == 1) {
        angle[i, j] = acos( -( inv_mat[j, 1:j - 1]' * inv_mat[i, 1:j - 1] ) /
                             ( inv_mat[j, j] * prod_sines ) );
        
        }
        inv_mat[i, j] = cos(angle[i, j]) * prod_sines;
      }
    }
  }
  return inv_mat;
}
data {
  int N;
  matrix[N, N] R;
  int<lower=0, upper=1> is_symmetric;
  int Z; // number of zeros
  int K;
  vector[K] where_zero;
}
transformed data {
 // int K = N * (N - 1) / 2;
  row_vector[N] m;
  
  for (n in 1:N)
    m[n] = N - n + 1;
}
parameters {
  vector<lower = 0, upper = pi()>[K - Z] theta_free;
} 
model {
   matrix[N, N] theta_mat = build_angle_mat(where_zero, theta_free, N);
   matrix[N, N] L = zero_constrain(theta_mat, N);
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

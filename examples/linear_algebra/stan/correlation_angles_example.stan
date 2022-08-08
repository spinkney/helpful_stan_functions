functions {
vector lower_elements(matrix M, int tri_size){
    int n = rows(M);
    int counter = 1;
    vector[tri_size] out;

    for (i in 2:n){
      for (j in 1:(i - 1)) {
        out[counter] = M[i, j];
        counter += 1;
      }
    }
    return out;
  }

 matrix angle2chol_lp(matrix angle_mat) {
  int N = rows(angle_mat);
  matrix[N, N] inv_mat = identity_matrix(N);
  
  inv_mat[ : , 1] = cos(angle_mat[ : , 1]);
  
  for (i in 2 : N) {
    inv_mat[i, i] = prod(sin(angle_mat[i, 1 : (i - 1)]));
    target += (N - i + 1) * log(sin(angle_mat[i:N, i - 1]));
    if (i > 2) {
      for (j in 2 : (i - 1)) {
        inv_mat[i, j] = cos(angle_mat[i, j]) * prod(sin(angle_mat[i, 1 : (j - 1)]));
      }
    }
  }
  return inv_mat;
}
matrix angle_vec2angle_mat(vector angle, int K) {
  int N = num_elements(angle);
  matrix[K, K] mat = add_diag(identity_matrix(K), rep_vector(-1, K));
  int count = K - 1;
  int pos = 1;

  for (k in 1 : K - 1) {
    mat[k + 1:K, k] = segment(angle, pos, count);
    pos += count;
    count -= 1;
  }
  
  return mat;
}

matrix tfp_cholesky_lp(vector raw, int N) {
    matrix[N, N] L = identity_matrix(N);
    int counter = 1;
    real s;
    
    for (i in 2 : N) {
      for (j in 1 : (i - 1)) {
        L[i, j] = raw[counter];
        counter += 1;
      }
      s = norm2(L[i,  : i]);
      L[i,  : i] = L[i,  : i] / s;
      target += -(i + 1) * log(s);
    }
    return L;
  }
}
data {
  int N;
  matrix[N, N] R;
}
transformed data {
  int K = (N * (N - 1)) %/% 2;
  row_vector[N] m;

  for (n in 1:N)
    m[n] = N - n + 1;
}
parameters {
  vector<lower = 0, upper = pi()>[K] theta;
  cholesky_factor_corr[N] L_stan;
  vector[K] raw;
}
transformed parameters {
  matrix[N, N] L = angle2chol_lp(angle_vec2angle_mat(theta, N));
  matrix[N, N] R_hat = multiply_lower_tri_self_transpose(L);
  matrix[N, N] L_tfp = tfp_cholesky_lp(raw, N);
  matrix[N, N] R_tfp_hat = multiply_lower_tri_self_transpose(L_tfp);
  matrix[N, N] R_stan = multiply_lower_tri_self_transpose(L_stan);
}
model {
  lower_elements(R, K) ~ normal(lower_elements(R_hat, K), 0.05);
  lower_elements(R, K) ~ normal(lower_elements(R_stan, K) , 0.05);
  lower_elements(R, K) ~ normal(lower_elements(R_tfp_hat, K), 0.05);
}
generated quantities {
    matrix[N, N] angle_mat = angle_vec2angle_mat(theta, N);
}
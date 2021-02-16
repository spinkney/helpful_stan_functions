functions {
 vector lower_elements(matrix M){
    int n = rows(M);
    int k = n * (n - 1) / 2;
    int counter = 1;
    vector[k] out;
 
    for (i in 2:n){
      for (j in 1:(i - 1)) {
        out[counter] = M[i, j];
        counter += 1;
      }
    }
    return out;
  }
  
   vector lower_elements_constrain(matrix M, int A){
    int n = rows(M);
    int counter = 1;
    vector[A] out;
 
    for (i in 2:n){
      for (j in 1:(i - 1)) {
        if(M[i, j] > 0){
         out[counter] = M[i, j];
         counter += 1;
        }
      }
    }
    return out;
  }
  
 matrix build_angle_mat (vector where_zero, vector angle_raw, int K) {
  int N = num_elements(where_zero);
  matrix[K, K] mat;
  int count = 1;
  int raw_count = 0;
  
  mat[1, 1] = 0.0;
  for (k in 2:K){
    mat[k, k] = 0.0;
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
 
 matrix zero_constrain (vector angle_raw, int[] csr_rows, int[] csr_cols, int Z, int N){
  matrix[N, N] inv_mat;
  int count = 1;
//  matrix[N, N] angle = angle_raw;
  
  inv_mat[1, 1] = 1;
  
  for (i in 2:N) {
    // constrain first column
    // if C = BB^T then for the first column
    // c_{i, 1} = b_{i, 1} = 0 so 
    // if c_{i, 1} = 0 then cos(\theta_{i, 1}) = 0
    // acos(0) = \theta_{i, 1} = pi / 2
     inv_mat[i, 1] = cos(angle[i, 1]);
   
    if (i > 2) {
      for (j in 2:(i - 1)) {
        real prod_sines = prod(sin(angle[i, 1:j - 1]));
        real cos_theta;
        if (is_nan(angle_raw[i, j]) == 1) {
          cos_theta = -(dot_product(inv_mat[j, 1:j - 1], inv_mat[i, 1:j - 1]) ) / ( inv_mat[j, j] * prod_sines );
            if ( cos_theta < -1 || cos_theta > 1 ) reject("cos_theta is ", cos_theta, " and must be in [-1, 1]"); // cos_theta = 0;
           // else if( cos_theta > 1) cos_theta = 1;
            angle[i, j] = acos( cos_theta );
          }
        inv_mat[i, j] = cos(angle[i, j]) * prod_sines;
      }
    }
    inv_mat[i, i] = prod(sin(angle[i, 1:(i - 1)]));
  }
  return inv_mat;
}
}
data {
  int N;
  matrix[N, N] R;
  matrix[N, N] R_lower_tri; // strictly lower tri elemenets /n
                            // all other values set to 0 including diagonal
  int<lower=0, upper=1> is_symmetric;
  int Z; // number of non-zero values in lower tri
  
  int K; // number of lower tri elements
  vector[K] where_zero; // vector of lower tri elements by row indicating /n
                        // with a 1 which elements are zero and 0 otherwise
}
transformed data {
 // int K = N * (N - 1) / 2;
  vector[Z] R_values = csr_extract_w(R_lower_tri);
  int col_index[Z] = csr_extract_v(R_lower_tri);
  int row_index[Z] = csr_extract_u(R_lower_tri);
  vector[K - Z] pi_over_two = rep_vector(pi() / 2, K - Z);
  row_vector[N] m;
  
  for (n in 1:N)
    m[n] = N - n + 1;
}
parameters {
  vector<lower = 0, upper = pi()>[Z] theta_free;
} 
transformed parameters {
  vector<lower = 0, upper = pi()>[K] theta;
  matrix[N, N] theta_mat = build_angle_mat(where_zero, theta_free, N);
  matrix[N, N] L = zero_constrain(theta_mat, N);
  matrix[N, N] R_hat = multiply_lower_tri_self_transpose(L);
  
  theta[where_zero] = pi_over_two;
  theta[1 - where_zero] = theta_free;
}
model {
  // log_absd_jacs 
  // sin(theta) is always > 0 since theta in (0, pi)
  // because of this the diagonal of L is always > 0 also
   target += 2 * sum(log(sin(csr_extract_w(theta_mat))));            // angle2chol
   target += N * log(2) + m * diagonal(log(L));   // cholesky tcrossprod
   if (is_symmetric == 1)
    lower_elements(R) ~ normal(lower_elements(R_hat), 0.001);
    else to_vector(R) ~ normal(to_vector(R_hat), 0.001);
}
generated quantities {
   matrix[N, N] theta_out = build_angle_mat(where_zero, theta_free, N);
   matrix[N, N] R_out = multiply_lower_tri_self_transpose( zero_constrain(theta_out, N) );
}

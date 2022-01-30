functions {
#include correlation_angles.stan
#include triangular.stan
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

 // set the upper triangle to 0
 // only looking at strictly lower tri part
 vector sparse_cholesky_lp (vector angle_raw, int[] csr_rows, int[] csr_cols, int Z, int N){
  vector[Z + N] sparse_chol; // Z + N is the number of non-zero values in lower tri plus
                             // the diagonal since the angles do not include the diagonal
  int R = size(csr_rows);
  int C = size(csr_cols);
  int S[R, C] = append_array(csr_rows, csr_cols);
  int count = 1;

  sparse_chol[count] = 1;

  // traversing in col-major order
  // S[2, 1]           skips S[2, 2]
  // S[3, 1], S[3, 2]  skips S[3, 3]
  // etc
  for (r in S) {
    int this_rows_column_num = 0; 
    for (c in r) {
      count += 1;
      this_rows_column_num += 1;

      if(this_rows_column_num == 1)
        sparse_chol[count] = cos(angle_raw[count]); // constrain first column


    }

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
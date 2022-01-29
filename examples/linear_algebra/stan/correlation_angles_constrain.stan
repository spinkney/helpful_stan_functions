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

 matrix zero_constrain (matrix angle_raw, int N){
  matrix[N, N] inv_mat;
  matrix[N, N] angle = angle_raw;

  inv_mat[1, 1] = 1;

  for (i in 2:N) {
    // constrain first column
    // if C = BB^T then for the first column
    // c_{i, 1} = b_{i, 1} = 0 so 
    // if c_{i, 1} = 0 then cos(\theta_{i, 1}) = 0
    // acos(0) = \theta_{i, 1} = pi / 2
     if (is_nan(angle_raw[i, 1]) == 1) {
       inv_mat[i, 1] = 0;
       angle[i, 1] = pi() / 2;
     } else inv_mat[i, 1] = cos(angle[i, 1]);

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
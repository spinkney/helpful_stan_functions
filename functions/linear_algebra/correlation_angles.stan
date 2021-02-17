  /* Angle matrix to Cholesky correlation
   *
   * Copyright Sean Pinkney, Feb. 10, 2021
   *
   * Peter J. Rousseeuw & Geert Molenberghs (1993) 
   * Transformation of non positive semidefinite correlation matrices, 
   * Communications in Statistics - Theory and Methods, 
   * 22:4, 965-984, DOI: 10.1080/03610928308831068
   *
   * @param angle_mat Matrix of angles 
   * @param Cholesky factor of correlation matrix
   */
matrix angle2chol (matrix angle_mat){
  int N = rows(angle_mat);
  matrix[N, N] inv_mat = identity_matrix(N);
    
  inv_mat[, 1] = cos(angle_mat[, 1]);
  
  for (i in 2:N) {
    inv_mat[i, i] = prod(sin(angle_mat[i, 1:(i - 1)]));
    if (i > 2) {
      for (j in 2:(i - 1)) {
        inv_mat[i, j] =
          cos(angle_mat[i, j]) * prod(sin(angle_mat[i, 1:(j - 1)]));
      }
    }
  }
  return inv_mat;
}

  /* Angle vector to angle matrix
   *
   * Copyright Sean Pinkney, Feb. 10, 2021
   *
   * @param angle Vector of angles length of the number of elements of the lower triangle elemements n/
     of the output angle matrix. If the size of the output matrix is K x K then the n/
     length is K * (K - 1) / 2
   * @param angle matrix
   */
matrix angle_vec2angle_mat (vector angle, int K) {
  int N = num_elements(angle);
  matrix[K, K] mat = add_diag(identity_matrix(K), rep_vector(-1, K));
  int count = K;
  int pos = 1;
    
  for (k in 1:K - 1){
    count -= 1;
    mat[k + 1:K, k] = segment(angle, pos, count);
    pos += count;
  }
    
  return mat;
}


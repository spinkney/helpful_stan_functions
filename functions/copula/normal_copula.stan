  /* Normal copula log density
   *
   * Copywrite Andre Pfeuffer, Sean Pinkney
   * https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/XdX3u1vDAAAJ
   * Accessed and modified Feb. 5, 2021 
   *
   * Meyer, Christian. "The Bivariate Normal Copula." 
   * arXiv preprint arXiv:0912.2816 (2009). Eqn 3.3.
   * accessed Feb. 6, 2021.
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param rho Real number [-1, 1]
   * @param log density
   */
real normal_copula(real u, real v, real rho) {
  real rho_sq = square(rho);
  
  return (0.5 * rho * (-2. * u * v + square(u) * rho + square(v) * rho)) /
      (-1. + rho_sq) - 0.5 * log1m(rho_sq);
}

  /* Normal copula log density vectorized
   *
   * Copywrite Sean Pinkney, Feb. 6, 2021
   *
   * Meyer, Christian. "The Bivariate Normal Copula." 
   * arXiv preprint arXiv:0912.2816 (2009). Eqn 3.3.
   * accessed Feb. 6, 2021.
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param rho Real number [-1, 1]
   * @param log density
   */
real normal_copula_vector(vector u, vector v, real rho){
   int N = num_elements(u);
   real rho_sq = square(rho);
   
   real a1 = 0.5 * rho;
   real a2 = rho_sq - 1;
   real a3 = 0.5 * log1m(rho_sq);
   real x = -2 * u' * v + rho * (dot_self(u) + dot_self(v));
  
  return a1 * x / a2 - N * a3;
}

  /* Multi-Normal Cholesky copula log density
   *
   * Copywrite Sean Pinkney, Feb. 8, 2021
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param rho Real number [-1, 1]
   * @param log density
   */

real multi_normal_copula_lpdf(matrix u, matrix L){
   int K = rows(u);
   int N = cols(u);
   real inv_sqrt_det_log = sum(log(diagonal(L)));
   matrix[K, N] x = mdivide_left_tri_low(L, u);

   return -N * inv_sqrt_det_log - 0.5 * sum(columns_dot_self(x) - columns_dot_self(u));
}

// slow version
// real multi_normal_copula_lpdf(matrix u, matrix L){
//    int K = rows(L);
//    int N = cols(u);
//    real inv_sqrt_det_log = sum(log(diagonal(L)));
//    matrix[K, K] L_inv1m = add_diag(chol2inv(L), rep_vector(-1.0, K));
//    real acc = 0;
//    
//    for (n in 1:N)
//       acc += quad_form(L_inv1m, u[, n]);
// 
//    return -N * inv_sqrt_det_log - 0.5 * acc;
// }
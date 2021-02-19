  /** @addtogroup normal Normal Copula Functions
   * \ingroup copula
   *  @{ */
   
 /**
   * Normal copula log density
   *
   * Copyright Andre Pfeuffer, Sean Pinkney 2017, 2021 \n
   * https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/XdX3u1vDAAAJ \n
   * Accessed and modified Feb. 5, 2021 
   *
   * \f$x = 100\f$
   *
   * Meyer, Christian. "The Bivariate Normal Copula." 
   * arXiv preprint arXiv:0912.2816 (2009). Eqn 3.3.
   * accessed Feb. 6, 2021.
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param rho Real number [-1, 1]
   * @return log density
   */
real normal_copula_lpdf(real u, real v, real rho) {
  real rho_sq = square(rho);
  
  return (0.5 * rho * (-2. * u * v + square(u) * rho + square(v) * rho)) /
      (-1. + rho_sq) - 0.5 * log1m(rho_sq);
}

 /**
   * Normal copula log density vectorized
   *
   * Copyright 2021, Sean Pinkney
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
real normal_copula_vector_lpdf(vector u, vector v, real rho){
   int N = num_elements(u);
   real rho_sq = square(rho);
   
   real a1 = 0.5 * rho;
   real a2 = rho_sq - 1;
   real a3 = 0.5 * log1m(rho_sq);
   real x = -2 * u' * v + rho * (dot_self(u) + dot_self(v));
  
  return a1 * x / a2 - N * a3;
}

 /**
   * Multi-Normal Cholesky copula log density
   *
   * Copyright 2021, Sean Pinkney
   *
   * @param u Matrix
   * @param L Cholesky factor matrix
   * @param log density
   */

real multi_normal_copula_lpdf(matrix u, matrix L){
   int K = rows(u);
   int N = cols(u);
   real inv_sqrt_det_log = sum(log(diagonal(L)));
   matrix[K, N] x = mdivide_left_tri_low(L, u);

   return -N * inv_sqrt_det_log - 0.5 * sum(columns_dot_self(x) - columns_dot_self(u));
}

 /**
   * Bivariate Normal Copula cdf
   *
   * Copyright 2021, Sean Pinkney
   *
   * @param u Vector size 2
   * @param rho Real on [-1, 1]
   * @param cumulative density
   */
real bivariate_normal_copula_cdf(vector u, real rho){
   real a = 1 / sqrt(1 - square(rho));
   real avg_uv = mean(u);
   real pu = inv_Phi(u[1]);
   real pv = inv_Phi(u[2]);
   real alpha_u = a * (pv / pu - rho);
   real alpha_v = a * (pu / pv - rho);
   real d = 0;
   
   if ( u[1] < 0.5 && u[2] >= 0.5 || u[1] >= 0.5 && u[2] < 0.5 )
      d = 0.5;
   
   return avg_uv - owens_t(pu, alpha_u) - owens_t(pv, alpha_v) - d;
}

 /** @} */
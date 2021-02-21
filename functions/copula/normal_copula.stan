  /** @addtogroup normal Normal Copula Functions
   * \ingroup copula
   *  @{ */
   
 /**
   * Normal copula log density
   *
   * \f[
   *   c(u,\,v;\, \rho) = \frac{1}{\sqrt{1 - \rho^2}} \exp \bigg(
   *  \frac{2\rho \Phi^{-1}(u) \Phi^{-1}(v) - \rho^2 (\Phi^{-1}(u)^2 + \Phi^{-1}(v)^2)}{2 (1 - \rho^2)}\bigg)
   * \f]
   *
   *  Meyer, Christian. "The Bivariate Normal Copula." 
   * arXiv preprint arXiv:0912.2816 (2009). Eqn 3.3.
   * accessed Feb. 6, 2021.
   *
   * @copyright Andre Pfeuffer, Sean Pinkney 2017, 2021
   * https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/XdX3u1vDAAAJ \n
   * Accessed and modified Feb. 5, 2021 
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param rho Real number (-1, 1)
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
   * @param rho Real number (-1, 1)
   * @return log density
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
   * \f{aligned}{
   *    c(\mathbf{u}) &= \frac{\partial^d C}{\partial \mathbf{\Phi_1}\cdots \partial \mathbf{\Phi_d}} \\
   *                   &= \frac{1}{\sqrt{\det L}} \exp \Bigg(-\frac{1}{2} 
   *                              \begin{pmatrix} \mathbf{\Phi}^{-1}(u_1) \\ vdots  \\ \mathbf{\Phi}^{-1}(u_d)  \end{pmatrix}^T
   *                              (L^{-1} - I)
   *                              \begin{pmatrix} \mathbf{\Phi}^{-1}(u_1) \\ vdots  \\ \mathbf{\Phi}^{-1}(u_d)  \end{pmatrix}
   *                              \Bigg)
   *\f}
   *
   * @copyright Sean Pinkney, 2021
   *
   * @param u Matrix
   * @param L Cholesky factor matrix
   * @return log density
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
   * Meyer, Christian. "The Bivariate Normal Copula." \n
   * arXiv preprint arXiv:0912.2816 (2009). Eqn 3.11. \n
   * accessed Feb. 6, 2021.
   *
   * @copyright Sean Pinkney, 2021
   *
   * @param u Vector size 2
   * @param rho Real on (-1, 1)
   * @return cumulative density
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
   /** @addtogroup multi_skewed_elliptical Multivariate Skewed Elliptical Distribution Functions
   *
   * The multivariate skew elliptical distributions take the form
   *  \f[
   *    f(\mathbf{y} \mid \mathbf{\mu}, \mathbf{\Sigma}, \mathbf{\lambda}) =
   *     2 f_g^k(\mathbf{y}) F_g_{q(\mathbf{y})}(\mathbf{\lambda})^T (\mathbf{y - \mu})
   *  \f]
   * where \f$ f_g^k(\mathbf{y}) \f$ is the pdf of the \f$k\f$-th marginal from the
   * elliptical class of dsitributions. These distributions are characterized by
   * the density generator function \f$ g^k(u), \, u /ge 0 \f$ such that
   * \f[
   *   \int_0^\infty u^{k/2 - 1} g^k(u) du = \frac{\Gamma(k/2)}{\pi^{k/2}} 
   * \f] and the pdf is of the following form
   * \f[
   *   f_{\mu, \Sigma}(x) = |\Sigma |^{-1/2} g^k((x - \mu)\Sigma^{-1} (x - mu)).
   * \f]
   *
   * \ingroup multivariate
   *  @{ */

   /**
   * Multivariate Skewed T Distribution
   *
   *  The density is defined as
   * \f[
   *  f_z(z) = 2 * t_d(z \mid \Omega, \nu) T \bigg( \alpha^T z \sqrt{\frac{\nu + d}{\nu + z^T \Omega^{-1} z}} \mid \nu + d \bigg)
   * \f]
   * where \f$ T(\cdot \mid \rho) \f$ denotes the univariate Student's \f$ t \f$ distribution on \f$ \rho \f$ d.f. 
   * @copyright Sean Pinkney 2021 
   *
   * @param z Matrix
   * @param L Cholesky factor of correlation matrix 
   * @param sigma Vector real \f$ > 0 \f$
   * @param lambda Vector
   * @param nu real \f$ > 0 \f$
   * @return integrand
   */
  real multi_skewed_t_lpdf(matrix z, matrix L, vector sigma, vector lambda, real nu) {
    int p = rows(z);
    int N = cols(z);
    row_vector[N] Q = columns_dot_self(mdivide_left_tri_low(L, z));
    row_vector[N] alpha = (lambda ./ sigma)' * z;
    real lp = N * log2() + 
                // multi_student_t lpdf
                N * ( lgamma( (p + nu) / 2) -
                    lgamma(nu / 2) - 0.5 * p * log(nu) - 0.5 * p * log(pi()) -
                    sum(log(diagonal(L)))) - 
                0.5 * sum( (nu + p) * log1p(Q / nu) );
    
    
    for (n in 1:N)
     lp += student_t_lcdf(alpha[n] * sqrt( (p + nu) / (nu + Q[n])) | p + nu, 0, 1);
     
    return lp;
  }
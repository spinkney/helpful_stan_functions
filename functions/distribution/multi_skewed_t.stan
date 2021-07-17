   /** @addtogroup multi_skewed_elliptical Multivariate Skewed Elliptical Distribution Functions
   *
   * The probability density function of the multivariate skewed t distribution is given by
   *  \f[
   *    f(\mathbf{x} \mid \mathbf{\mu}, \mathbf{\omega}) \left(\prod_{i=1}^c \binom{m_i}{x_i} \right) \int_0^1 \prod_{i=1}^c (1-t^{\omega_i/D})^{x_i} \operatorname{d}t
   *  \f]
   * where \f$ \mathbf{m} = (m_1, \ldots, m_c) \in \mathbb{N}^c \f$ is the size of each $c$ group
   * of the population. The population is given by \f$ N = \sum_{i=1}^c m_i \f$ where the realizations
   * from each group is given by \f$ \mathbf{x} \f$. The weights of each group are contained in the \f$c\f$-sized
   * simplex \f$ \mathbf{\omega} \f$. Where \f$ D = \mathbf{\omega} \cdot (\mathbf{m} - \mathbf{x})) = \sum_{i=1}^c 
   * \omega_i(m_i - x_i) \f$. 
   * 
   * \ingroup multivariate
   *  @{ */

   /**
   * Multivariate Skewed T Distribution
   *
   * @copyright Sean Pinkney 2021 
   *
   * @param t Real number on [0,1]
   * @param xc 
   * @param theta Array of real parameters
   * @param x_r Array of data (real) 
   * @param x_i Array of data (integer)
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
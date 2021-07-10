  /** @addtogroup cholesky Cholesky Correlation Matrix Constraint functions
   *
   * Original code is in Stan-math 
   * https://github.com/stan-dev/math/blob/develop/stan/math/prim/fun/cholesky_corr_constrain.hpp#L48
   *
   * \ingroup correlation
   *  @{ */
  /**
   * Correlation Constraint
   *
   * See https://mc-stan.org/docs/2_27/reference-manual/cholesky-factors-of-correlation-matrices-1.html
   *
   * @copyright Sean Pinkney, 2021 
   * @author Sean Pinkney
   * @param x Vector 
   * @return vector \f$\in (-1, 1)\f$
   */
  vector corr_constrain_lp(vector x) {
    int N = num_elements(x);
    vector[N] tanh_x = tanh(x);
    
    target += sum(log1m(square(tanh_x)));
    
    return tanh_x;
  }
  /**
   * Cholesky Correlation Matrix Constraint
   *
   * See https://mc-stan.org/docs/2_27/reference-manual/cholesky-factors-of-correlation-matrices-1.html
   *
   * @copyright Sean Pinkney, 2021 
   * @author Sean Pinkney
   * @param y Vector 
   * @param K Int 
   * @return Cholesky factor of correlation matrix
   */
 matrix cholesky_corr_constrain_lp (vector y, int K) {
    int k_choose_2 = (K * (K - 1)) / 2;
    vector[k_choose_2] z = corr_constrain_lp(y);
    matrix[K, K] x = rep_matrix(0, K, K);
    int iter = 0;
    
    x[1, 1] = 1;

    for (i in 2:K) {
      real sum_sqs;
      iter += 1;
      x[i, 1] = z[iter];
      sum_sqs = square(x[i, 1]);
        if (i > 2) {
            for (j in 2:i - 1) {
                iter += 1;
                target += 0.5 * log1m(sum_sqs);
                x[i, j] = z[iter] * sqrt(1.0 - sum_sqs);
                sum_sqs += square(x[i, j]);
            }
        }
        x[i, i] = sqrt(1.0 - sum_sqs);
    }
  return x;
  }
   /** @} */
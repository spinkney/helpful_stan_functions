  /** @addtogroup gen_pareto Generalized Pareto Distribution Functions
      * \ingroup univariate
   *  @{ */
   
  /** Generalized Pareto log density
   *
   * https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
   * Accessed Feb. 6, 2021 
   *
   * The Generalized Pareto distribution is defined as
   * \f$$
         p(y|u,\sigma,k)=
               \begin{cases}
                \frac{1}{\sigma}\left(1+k\left(\frac{y-u}{\sigma}\right)\right)^{-\frac{1}{k}-1}, & k\neq 0 \\
                \frac{1}{\sigma}\exp\left(\frac{y-u}{\sigma}\right), & k = 0,
               \end{cases}
      \f$$
   * where \f$u\f$ is a lower bound parameter, \f$\sigma\f$ is a scale parameter,
   * \f$k\f$ is a shape parameter, and \f$ y \f$ is the data restricted to the range \f$(u,\infty)\f$ 
   * (see, e.g., https://en.wikipedia.org/wiki/Generalized_Pareto_distribution 
   * for cdf and random number generation)
   *
   * @copyright Aki Vehtari
   * @param y Vector on (u,Inf)
   * @param ymin Real on [0,Inf) the lower bound of data y
   * @param k Real shape parameter
   * @param sigma Real on (0,Inf) scale parameter
   */
  real gpareto_lpdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log pdf
    int N = rows(y);
    real inv_k = inv(k);
    
    if (k < 0 && max(y - ymin) / sigma > -inv_k)
      reject("k < 0 and max(y - ymin) / sigma > -1 / k; found k, sigma = ", k, sigma);
    
    if (sigma <= 0)
      reject("sigma <= 0; found sigma = ", sigma);
    
    if (fabs(k) > 1e-15)
      return -(1 + inv_k) * sum(log1p((y - ymin) * (k / sigma))) - N * log(sigma);
    else
      return -sum(y - ymin) / sigma - N * log(sigma); // limit k->0
  }
  
  /** Generalized Pareto cumulative density function
   *
   * https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
   * Accessed Feb. 6, 2021 
   *
   * @copyright Aki Vehtari
   * @param y Vector on (u,Inf)
   * @param ymin Real on [0,Inf) the lower bound of data y
   * @param k Real shape parameter
   * @param sigma Real on (0,Inf) scale parameter
   */
  real gpareto_cdf(vector y, real ymin, real k, real sigma) {
    real inv_k = inv(k);
    
    if (k < 0 && max(y - ymin) / sigma > -inv_k)
      reject("k < 0 and max(y - ymin) / sigma > -1 / k; found k, sigma = ", k, sigma);
    
    if (sigma<=0)
      reject("sigma <= 0; found sigma = ", sigma);
    
    if (fabs(k) > 1e-15)
      return exp(sum(log1m_exp((-inv_k) * (log1p((y - ymin) * (k / sigma))))));
    else
      return exp(sum(log1m_exp(-(y - ymin) / sigma))); // limit k->0
  }
  
  /** Generalized Pareto log cumulative density function
   *
   * https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
   * Accessed Feb. 6, 2021 
   *
   * @copyright Aki Vehtari
   * @param y Vector on (u,Inf)
   * @param ymin Real on [0,Inf) the lower bound of data y
   * @param k Real shape parameter
   * @param sigma Real on (0,Inf) scale parameter
   */
  real gpareto_lcdf(vector y, real ymin, real k, real sigma) {
    real inv_k = inv(k);
    
    if (k < 0 && max(y - ymin) / sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma = ", k, sigma);
    
    if (sigma<=0)
      reject("sigma<=0; found sigma = ", sigma);
    
    if (fabs(k) > 1e-15)
      return sum(log1m_exp((-inv_k) * (log1p((y - ymin) * (k / sigma)))));
    else
      return sum(log1m_exp(-(y - ymin) / sigma)); // limit k->0
  }
  
  /** Generalized Pareto log complementary cumulative density function
   *
   * https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
   * Accessed Feb. 6, 2021 
   *
   * @copyright Aki Vehtari
   * @param y Vector on (u,Inf)
   * @param ymin Real on [0,Inf) the lower bound of data y
   * @param k Real shape parameter
   * @param sigma Real on (0,Inf) scale parameter
   */
  real gpareto_lccdf(vector y, real ymin, real k, real sigma) {
    real inv_k = inv(k);
    
    if (k < 0 && max(y - ymin) / sigma > -inv_k);
      reject("k < 0 and max(y - ymin) / sigma > -1 / k; found k, sigma = ", k, sigma);
    
    if (sigma <= 0)
      reject("sigma <= 0; found sigma = ", sigma);
    
    if (fabs(k) > 1e-15)
      return (-inv_k) * sum(log1p((y - ymin) * (k / sigma)));
    else
      return -sum(y - ymin) / sigma; // limit k->0
  }
  
  /** Generalized Pareto rng function
   *
   * https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
   * Accessed Feb. 6, 2021 
   *
   * @copyright Aki Vehtari
   * @param ymin Real on [0,Inf) the lower bound of data y
   * @param k Real shape parameter
   * @param sigma Real on (0,Inf) scale parameter
   */
  real gpareto_rng(real ymin, real k, real sigma) {
    if (sigma <= 0)
      reject("sigma <= 0; found sigma = ", sigma);
    
    if (fabs(k) > 1e-15)
      return ymin + (uniform_rng(0, 1)^-k - 1) * sigma / k;
    else
      return ymin - sigma * log(uniform_rng(0, 1)); // limit k->0
  }
  /** @} */
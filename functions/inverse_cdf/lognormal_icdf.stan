  /** @ingroup icdf
   * **Lognormal Inverse CDF**
   * 
   * \f{aligned}{
   *  F^{-1}(x) &= \exp(\mu + \sqrt{2} \sigma \text{inv_erf}(2p - 1)) \\
   *         &= \exp(\mu + \sigma \Phi^{−1}(x))
   * \f}
   *  because
   *  \f{aligned}{
        \Phi(x) &= 1 / (2 * \pi) * \int e^{t^2} \; dt = 0.5 + 0.5 \, \text{erf}(x / \sqrt{2})  \\
   *  \implies \Phi^{−1}(x) &= \sqrt{2} \, \text{inv_erf}(2x − 1)
   * \f}
   
   * @author Sean Pinkney
   * @param p Real on \f$[0,\,1]\f$
   * @param mu Real \f$(-\infty, +\infty)\f$
   * @param sigma Real \f$(0, +\infty)\f$
   * @return inverse CDF value
   * @throws reject if \f$ p \notin [0, 1] \f$ 
   */
real lognormal_icdf (real p, real mu, real sigma){
   if (is_nan(p) || p < 0 || p > 1)
     reject("lognormal_icdf: p must be between 0 and 1; ",
           "found p = ", p);
   if (is_nan(mu) || is_inf(mu) )
      reject("lognormal_icdf: mu must be finite; ",
           "found mu = ", mu);
   if (is_nan(sigma) || is_inf(sigma) || sigma <= 0)
      reject("lognormal_icdf: sigma must be finite and > 0; ",
           "found sigma = ", sigma);
           
    return exp( mu + sigma * inv_Phi(p) );
  }
  
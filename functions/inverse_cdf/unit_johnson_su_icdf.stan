  /** @ingroup icdf
   * Unit Johnson SU Inverse CDF
   * 
   * \f{aligned}{
   *  F^{-1}(u, \, \mu, \, \sigma) &= \text{inv_logit}\bigg[ \sinh\bigg(\frac{\Phi^{-1}(u) - \mu}{\sigma}\bigg) \bigg]\\
   * \f}
   *
   * @author Sean Pinkney
   * @param u Real on \f$[0,\,1]\f$
   * @param mu Real \f$(-\infty, +\infty)\f$
   * @param sigma Real \f$(0, +\infty)\f$
   * @return inverse CDF value
   * @throws reject if \f$ p \notin [0, 1] \f$ 
   */
  real unit_johnson_icdf (real u, real mu, real sigma){
   if (is_nan(u) || u < 0 || u > 1)
     reject("unit_johnson_icdf: p must be between 0 and 1; ",
           "found u = ", u);
   if (is_nan(mu) || is_inf(mu) )
      reject("unit_johnson_icdf: mu must be finite; ",
           "found mu = ", mu);
   if (is_nan(sigma) || is_inf(sigma) || sigma <= 0)
      reject("unit_johnsonl_icdf: sigma must be finite and > 0; ",
           "found sigma = ", sigma);
    real x = sinh( (inv_Phi(u) - mu) / sigma);
    return inv_logit(x);
  }
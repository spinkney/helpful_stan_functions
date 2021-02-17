   /* Lognormal Inverse CDF
   *
   * Copyright Sean Pinkney, Feb. 2021
   *
   * icdf(x) = exp(μ + √2 σ inv_erf(2p - 1))
   *         = exp(μ + σ Φ^{−1}(x))
   *
   *  because
   *  Φ(x)= 1 / (2 * pi) * ∫ e^{t^2} dt = 0.5 + 0.5 * erf(x / √2)
   *  => Φ^{−1}(x) = √2 inv_erf(2x − 1)
   * 
   * @param p Real on [0,1]
   * @param mu Real (-∞, +∞)
   * @param sigma Real (0, +∞)
   * @param return inverse error result
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
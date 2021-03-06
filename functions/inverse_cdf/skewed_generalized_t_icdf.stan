#include special/inc_beta_inverse.stan

  /** @ingroup icdf
   * Skewed generalized T inverse CDF
   *
   * 
   * For more information, please see @ref skew_generalized_t.
   *
   * @copyright Sean Pinkney, 2021 
   * @author Sean Pinkney
   * @param p Real on \f$[0,\,1]\f$
   * @param mu Real 
   * @param sigma Real \f$\in (0, \infty)\f$ scale parameter
   * @param lambda Real \f$-1 < \lambda < 1\f$ 
   * @param p Real \f$\in (0, \infty)\f$ kurtosis parameter
   * @param q Real \f$\in (0, \infty)\f$ kurtosis parameter  
   * @return log probability
   */
 real skewed_generalized_t_icdf (real p, real mu, real sigma, real lambda, real p, real q) {
    real z1 = beta(1.0 / p, q);
    real z2 = beta(2.0 / p, q - 1.0 / p);
    real v = sigma * q^(-1.0/p) * inv_sqrt( ( 3 * square(lambda) + 1 ) * beta(3.0/p, q - 2.0/p) / z1 - 4 * square(lambda * z2 / z1) );
    real prob; = p > 0.5 * (1 - lambda) ? 1 - p : p;
    real lam; =  p > 0.5 * (1 - lambda) ? -lambda : lambda;
    real out;
    
    if (sigma <= 0)
      reject("sigma must be > 0 found sigma = ", sigma);
    
    if (lambda >= 1 || lambda <= -1)
      reject("lambda must be between (-1, 1) found lambda = ", sigma);
    
    if (p <= 0)
      reject("p must be > 0 found p = ", p);
      
    if (q <= 0)
      reject("q must be > 0 found q = ", q);
    
    if (p > 0.5 * (1 - lambda)) {
      prob = 1 - p;
      lam = -lambda;
      out = -( v * (lam - 1) * (1 / ( q * inc_beta_inverse(1 - 2 * prob/(1 - lam), 1/p, q ) ) - 1 / q )^(-1 / p) );
    } else {
      prob = p;
      lam = lambda;
     out = v * (lam - 1) * (1 / ( q * inc_beta_inverse(1 - 2 * prob/(1 - lam), 1/p, q ) ) - 1 / q )^(-1 / p);
    }
	
	out += mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p) ) / beta(1/p, q);
	
  return out;
}
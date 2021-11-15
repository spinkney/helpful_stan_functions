 /** @addtogroup skew_generalized_t Skew Generalized T distribution functions
   *
   * From the sgt R package 
   * Carter Davis (2015). sgt: Skewed Generalized T Distribution Tree. R package version 2.0. 
   * https://CRAN.R-project.org/package=sgt
   *
   * The Skewed Generalized T Distribution is a univariate 5-parameter distribution introuced 
   * by Theodossiou (1998) and known for its extreme flexibility. Special and limiting cases of
   * the SGT distribution include the skewed generalized error distribution, the generalized t 
   * distribution introduced by McDonald and Newey (1988), the skewed t proposed by Hansen (1994), 
   * the skewed Laplace distribution, the generalized error distribution (also known as the
   * generalized normal distribution), the skewed normal distribution, the student t distribution,
   * the skewed Cauchy distribution, the Laplace distribution, the uniform distribution, the 
   * normal distribution, and the Cauchy distribution. 
   * 
   * Hansen, B. E., 1994, Autoregressive Conditional Density Estimation, International Economic 
   * Review 35, 705-730.
   *
   * Hansen, C., J. B. McDonald, and W. K. Newey, 2010, Enstrumental Variables Estimation with 
   * Flexible Distribution sigma  Journal of Business and Economic Statistics 28, 13-25.
   *
   * McDonald, J. B. and W. K. Newey, 1988, Partially Adaptive Estimation of Regression Models 
   * via the Generalized t Distribution, Econometric Theory 4, 428-457.
   *
   * Theodossiou, Panayioti sigma  1998, Financial Data and the Skewed Generalized T 
   * Distribution, Management Science 44, 1650-166
   *
   * \ingroup univariate
   *  @{ */
   
  /**
   * The Skewed Generalized T distribution is defined as
   *
   * \f[
   * f_{SGT}(x; \mu, \sigma, \lambda, p, q) = \frac{p}{2 v \sigma  q^{1/p} B(\frac{1}{p},q)
   * \left(\frac{| x-\mu + m |^p}{q (v \sigma)^p (\lambda sign(x-\mu + m)+1)^p}+1\right)^{\frac{1}{p}+q}}
   * \f]
   * where  \f$B\f$ is the beta function, \f$ \mu \f$ is the location parameter, \f$\sigma > 0\f$
   * is the scale parameter, \f$-1 < \lambda < 1\f$ is the skewness parameter, and \f$p > 0\f$ and 
   * \f$q > 0\f$ are the parameters that control the kurtosis. \f$m\f$ and \f$v\f$ are not parameter sigma  
   * but functions of the other parameters that are used here to scale or shift the distribution 
   * appropriately to match the various parameterizations of this distribution.
   *
   * In the original parameterization Theodossiou of the skewed generalized t distribution, 
   * \f[
   * m = \frac{2 v \sigma \lambda q^{\frac{1}{p}} B(\frac{2}{p},q-\frac{1}{p})}{B(\frac{1}{p},q)}
   * \f]
   * and
   * \f[
   * v = \frac{q^{-\frac{1}{p}}}{\sqrt{ (3 \lambda^2 + 1)
   *  \frac{ B ( \frac{3}{p}, q - \frac{2}{p} )}{B (\frac{1}{p}, q )} 
   *  -4 \lambda^2 \frac{B ( \frac{2}{p}, q - \frac{1}{p} )^2}{ B (\frac{1}{p}, q )^2}}}.
   * \f]
   *
   * @copyright Sean Pinkney, 2021 
   * @author Sean Pinkney
   * @param x Vector 
   * @param mu Real 
   * @param sigma Real \f$\in (0, \infty)\f$ scale parameter
   * @param lambda Real \f$-1 < \lambda < 1\f$ 
   * @param p Real \f$\in (0, \infty)\f$ kurtosis parameter
   * @param q Real \f$\in (0, \infty)\f$ kurtosis parameter  
   * @return log probability
   */
 real skew_generalized_t_lpdf(vector x, real mu, real  sigma,  real lambda, real p, real q) { 
    int N = num_elements(x);
    real z1 = beta(1.0 / p, q);
    real z2 = beta(2.0 / p, q - 1.0 / p);
    real v = q^(-1.0/p) * inv_sqrt( ( 3 * square(lambda) + 1 ) * beta(3.0/p, q - 2.0/p) / z1 - 4 * square(lambda * z2 / z1) ) ;
    real m = 2 * v * sigma * lambda * q^(1.0/p) * z2 / z1;
    real out1 = N * (log(p) - log(v) - log(sigma) - lmultiply(1/p, q) - lbeta(1.0 / p, q));
    real out2 = 0;
    
    if (sigma <= 0)
      reject("sigma must be > 0 found sigma = ", sigma);
    
    if (lambda >= 1 || lambda <= -1)
      reject("lambda must be between (-1, 1) found lambda = ", sigma);
    
    if (p <= 0)
      reject("p must be > 0 found p = ", p);
      
    if (q <= 0)
      reject("q must be > 0 found q = ", q);
      
    for (n in 1:N) {
      real r = x[n] - mu + m;
      if (r < 0)
        out2 -= log1p((-r)^p / (q * (v * sigma * (1 - lambda))^p));
      else
        out2 -= log1p(r^p / (q * (v * sigma * (1 + lambda))^p));
    }
    return fma(out2, (1/p + q), out1);
  }

 /** Skew Generalized T log cumulative density function
   *
   * @copyright Sean Pinkney, 2021 
   * @author Sean Pinkney
   * @param x Real 
   * @param mu Real 
   * @param sigma Real \f$\in (0, \infty)\f$ scale parameter
   * @param lambda Real \f$-1 < \lambda < 1\f$ 
   * @param p Real \f$\in (0, \infty)\f$ kurtosis parameter
   * @param q Real \f$\in (0, \infty)\f$ kurtosis parameter  
   * @return log probability
   */
  real skewed_generalized_t_lcdf (real x, real mu, real sigma, real lambda, real p, real q) {
    real z1 = beta(1.0 / p, q);
    real z2 = beta(2.0 / p, q - 1.0 / p);
    real v = q^(-1.0/p) * inv_sqrt( ( 3 * square(lambda) + 1 ) * beta(3.0/p, q - 2.0/p) / z1 - 4 * square(lambda * z2 / z1) ) ;
    real m = 2 * v * sigma * lambda * q^(1.0/p) * z2 / z1;
    real r = x - mu + m;
    real lambda_new;
    real r_new;
    
    if (sigma <= 0)
      reject("sigma must be > 0 found sigma = ", sigma);
    
    if (lambda >= 1 || lambda <= -1)
      reject("lambda must be between (-1, 1) found lambda = ", sigma);
    
    if (p <= 0)
      reject("p must be > 0 found p = ", p);
      
    if (q <= 0)
      reject("q must be > 0 found q = ", q);
      
    if (r > 0) {
      lambda_new = -lambda;
      r_new = -r;
    } else {
      lambda_new = lambda;
      r_new = r;
    }
    
    return  log(0.5) + log( (1 - lambda_new) + (lambda_new - 1) * beta_cdf(1 / ( 1 + q * (sigma * (1-lambda_new)/(-r_new) )^p) | 1 / p, q) );
}

 /** @} */

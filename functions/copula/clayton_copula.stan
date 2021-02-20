  /** @addtogroup clayton Clayton Copula Functions
   *
   * The Clayton copula is defined as 
   * \f[
   *      C(u,v) = (u^{-\theta} + v^{-\theta} - 1)^{-1/\theta}
   * \f]
   * for \f$0 < \theta < \infty\f$. 
   * 
   * \ingroup copula
   *  @{ */
   
  /** 
   * Clayton Bivariate Copula Log Density
   *
   * The probability density function is defined as
   * \f[
   *     c(u, v) = \frac{ \partial^2 C(u, v) }{ \partial u \partial v} = (\theta + 1)(uv)^{-(\theta + 1)}(u^{-\theta} + v^{-\theta} - 1)^{-\frac{2 \theta + 1}{\theta}}
   * \f]
   *
   * Copyright Andre Pfeuffer, Sean Pinkney \n
   * https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/XdX3u1vDAAAJ \n
   * Accessed and modified Feb. 5, 2021 
   *
   * The copula is defined \f$0 < \theta < \infty\f$
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param theta Real number on (0, Inf)
   * @return log density
   * @throws reject if theta <= 0
   */
real clayton_copula_lpdf(real u, real v, real theta) {
    if(theta <= 0.0) 
      reject("clayton_copula: theta must > 0");
    
     return log1p(theta) 
            - (theta + 1) * (log(u) + log(v))
            - (1 + 2 * theta) / theta 
              * log(pow(u, -theta) + pow(v, -theta) - 1);
  }

 /** @} */
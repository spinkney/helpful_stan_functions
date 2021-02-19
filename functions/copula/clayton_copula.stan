/** \defgroup copula Copula Functions
  *
  * A copula is a multivariate cumulative distribution function for which the 
  * marginal probability distribution of each variable is uniform on the interval \f$[0, 1]\f$.
  * Copulas are used to describe the dependence between random variables.
  *
  **/
  
  /** @addtogroup clayton  Clayton Copula Functions
   * \ingroup copula
   *  @{ */
   
  /** 
   * Clayton Bivariate Copula Log Density
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
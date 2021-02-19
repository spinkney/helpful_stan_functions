  /** @addtogroup frank Frank Copula Functions
   * \ingroup copula
   *  @{ */
   
  real frank_copula_lpdf(real u, real v, real theta) {
    return log(theta) + log1m(exp(-theta))-theta*(u+v)
    - 2* log(1 - exp(-theta) - (1-exp(-theta*u))*(1-exp(-theta*v)));
  }
  
   /** @} */
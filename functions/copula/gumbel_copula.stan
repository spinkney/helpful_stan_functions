  /* Gumbel copula log density
   *
   * Copyright Ben Goodrich
   * https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/UZURBv2_AAAJ
   * Accessed Feb. 5, 2021 
   *
   * Archimedean family of Gumbel with parametric generator
   *        exp(−t^1/θ), t ∈ [0, ∞]
   * with θ ∈ [1, ∞). The range of admissible Kendall’s tau, 
   * as well as that of the upper taildependence coefficient, is [0,1). 
   * Note that this copula does not allow for lower tail dependence.
   * 
   * https://cran.r-project.org/web/packages/copula/copula.pdf
   * "copula" R package in the copFamilies documentation accessed Feb. 5, 2021
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param theta Real number >= 1, will throw otherwise
   * @param log density
   */
  real gumbel_copula(real u, real v, real theta) {
    
    real neg_log_u = -log(u); 
    real log_neg_log_u = log(neg_log_u);
    real neg_log_v = -log(v); 
    real log_neg_log_v = log(neg_log_v);
    real log_temp = log_sum_exp(theta * log_neg_log_u,
                                theta * log_neg_log_v);
    real theta_m1 = theta - 1;
    
    if (theta < 1) reject("theta must be >= 1");
    
    if (is_inf(theta)) {
      if (u == v) return 0;
      else return negative_infinity();
    }
    
    return theta_m1 * log_neg_log_u + theta_m1 * log_neg_log_v + 
            neg_log_u + neg_log_v - exp(log_temp / theta) + 
            log_sum_exp( 2 * theta_m1  / -theta * log_temp, log(theta_m1) +
                  (1 - 2 * theta) / theta * log_temp);
  }
  
  
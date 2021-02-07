  /* Clayton copula log density
   *
   * Copywrite Andre Pfeuffer, Sean Pinkney
   * https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/XdX3u1vDAAAJ
   * Accessed and modified Feb. 5, 2021 
   *
   *
   * @param u Real number on (0,1], not checked but function will return NaN
   * @param v Real number on (0,1], not checked but function will return NaN
   * @param theta Real number on (0, Inf)
   * @param log density
   */
real clayton_copula(real u, real v, real theta) {
    if(theta == 0.0) 
      reject("clayton_copula: theta must != 0");
    
     return log1p(theta) 
            - (theta + 1) * (log(u) + log(v))
            - (1 + 2 * theta) / theta 
              * log(pow(u, -theta) + pow(v, -theta) - 1);
  }


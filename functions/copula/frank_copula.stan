  real frank_copula(real u, real v, real theta) {
    return log(theta) + log1m(exp(-theta))-theta*(u+v)
    - 2* log(1 - exp(-theta) - (1-exp(-theta*u))*(1-exp(-theta*v)));
  }
real normal_copula(real u, real v, real rho) {
  real rho_sq = square(rho);
  
  return (0.5 * rho * (-2. * u * v + square(u) * rho + square(v) * rho)) /
      (-1. + rho_sq) - 0.5 * log1m(rho_sq);
  

  // return (0.5 * rho * (-2. * u * v + 1. * square(u) * rho + 1. * square(v) * rho)) /
  //     (-1. + rho_sq) - 0.5 * log1m(rho_sq);
  // 
  }

real normal_copula_vector(vector u, vector v, real rho){
   real rho_sq = square(rho);
   
   real a1 = -0.5 * log1m(rho_sq);
   
   real a2 = 2 * rho * u' * v;
   real a3 = -rho_sq * (dot_self(u) + dot_self(v));
   real a4 = 1 / (2 * (1 - rho_sq));
  
  // real a = 0.5 * rho * (-2 * u' * v + rho * (dot_self(u) + dot_self(v)));
  // real b = (-1. + square(rho)) - 0.5 * log1m(square(rho));
  // 
  return a1 + a4 * (a2 - a3);
}

functions {
real bivariate_normal_copula(real u, real v, real rho){
   real a = 1 / sqrt(1 - square(rho));
   real avg_uv = 0.5 * (u + v);
   real pu = inv_Phi(u);
   real pv = inv_Phi(v);
   real alpha_u = a * (pv / pu - rho);
   real alpha_v = a * (pu / pv - rho);
   real d = 0;
   
   if ( u < 0.5 && v >= 0.5 || u >= 0.5 && v < 0.5 )
      d = 0.5;
   
   return avg_uv - owens_t(pu, alpha_u) - owens_t(pv, alpha_v) - d;
}
}

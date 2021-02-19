
  /* 1d interpolation using cubic Hermite splines
  *
  * Useful for interpolating solutions of autonomous ODEs. Assumes that
  * derivative_fun() is defined in the functions block before this function.
  * It must have signature
  *
  * vector derivative_fun(real t, vector y, int[] a0, vector a1);
  *  
  * i.e. same as what can be used with ODE solvers (see
  * https://mc-stan.org/docs/2_26/stan-users-guide/coding-the-ode-system-function.html
  * ), but it shouldn't actually use the t argument. It will be always called
  * with t = 0.0.
  *
  * Info: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
  *
  * Author: Juho Timonen
  *
  * @param y array of D-vectors, length N_in
  * @param x incresing array of reals, length N_in
  * @param x_out increasing array of reals, length N_out, values 
  * must be in (min(x), max(x)]
  * @param a0 array of integer inputs given to derivative_fun()
  * @param theta parameter vector given to derivative_fun()
  * @return array of D-vectors, length N_out, corresponding to
  * interpolated values y_out
  */
vector[] interp_1d_cubic(vector[] y, data real[] x, data real[] x_out,
    int[] a0, vector theta){
  int left = 1;
  int right = 1;
  real h = 0.0;
  real w = 0.0;
  int N_in = size(x);
  int N_out = size(x_out);
  int D = size(y[1]);
  vector[D] f_left;
  vector[D] f_right;
  real h00; real h10; real h01; real h11;
  vector[D] y_out[N_out];
  for (j in 1:N_out) {
    
    // Find left and right point indices
    while(x[right] < x_out[j]) {
      right = right + 1;
    }
    while(x[left+1] < x_out[j]) {
      left = left + 1;
    }
    
    // Evaluate derivatives
    f_left = derivative_fun(0.0, y[left], a0, theta);
    f_right = derivative_fun(0.0, y[right], a0, theta);
    
    // Hermite basis functions
    h = x[right] - x[left];
    w = (x_out[j] - x[left]) / h;
    h00 = 2.0 * w^3 - 3.0 * w^2 + 1.0;
    h10 = w^3 - 2.0 * w^2 + w;
    h01 = -2.0 * w^3 + 3.0 * w^2;
    h11 = w^3 - w^2;
    
    // Compute interpolation
    y_out[j] = h00 * y[left] + h10 * h * f_left + h01 * y[right] + h11 * h * f_right;
  }
  return(y_out);
}

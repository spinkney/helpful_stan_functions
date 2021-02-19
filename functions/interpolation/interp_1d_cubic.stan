  /** \defgroup interpolation Interpolation Functions
  *
  * Interpolation is a type of estimation, a method of constructing new data points within the range of a discrete set of known data points.
  *
  **/
  
  /** @addtogroup 1d_interpolation 1-dimensional Interpolation Functions
   * \ingroup interpolation
   *  @{ */

 /** 
  * 1d interpolation using cubic Hermite splines
  *
  * Useful for interpolating solutions of autonomous ODEs. Assumes that
  * derivative_fun() is defined in the functions block before this function.
  * It must have signature
  *
  *     vector derivative_fun(real t, vector y, vector theta, data real[] x_r, 
  *      data int[] x_i),
  *  
  * i.e. same as what can be used with ODE solvers, but it shouldn't actually
  * use the t argument. 
  * Info: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
  *
  * @author Juho Timonen
  *
  * @param y array of D-vectors, length N_in
  * @param x incresing array of reals, length N_in
  * @param x_out increasing array of reals, length N_out, values 
  * must be in (min(x), max(x)]
  * @param theta parameter vector given to derivative_fun()
  * @param x_r array of real inputs given to derivative_fun()
  * @param x_i array of integer inputs given to derivative_fun()
  * @return array of D-vectors, length N_out, corresponding to
  * interpolated values y_out
  */
vector[] interp_1d_cubic(vector[] y, data real[] x, data real[] x_out,
    vector theta, data real[] x_r, data int[] x_i){
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
    f_left = derivative_fun(0.0, y[left], theta, x_r, x_i);
    f_right = derivative_fun(0.0, y[right], theta, x_r, x_i);
    
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

 /** @} */

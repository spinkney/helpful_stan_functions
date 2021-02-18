
  /* Explicit midpoint method
  *
  * Assumes that derivative_fun() is defined in the functions block before
  * this function. It must have signature
  *
  * vector derivative_fun(real t, vector y, vector theta, data real[] x_r, 
  *     data int[] x_i),
  *  
  * i.e. same as what can be used with the built-in ODE solvers.
  *
  * Author: Juho Timonen
  *
  * @param t0 initial time
  * @param y0 initial state, D-vector
  * @param h step size (positive)
  * @param num_steps number of steps to take
  * @param theta parameter vector given to derivative_fun()
  * @param x_r array of real inputs given to derivative_fun()
  * @param x_i array of integer inputs given to derivative_fun()
  * @return array of D-vectors, length equal to num_steps + 1
  */
vector[] odeint_midpoint(real t0, vector y0, data real h, data int num_steps,
    vector theta, data real[] x_r, data int[] x_i){

  int d = num_elements(y0);
  vector[d] y[num_steps+1];
  vector[d] y_mid;
  real t = t0;
  y[1] = y0;
  
  // Integrate at grid of time points
  for(i in 1:num_steps){
    
    // Half-Euler step
    y_mid = y[i] + 0.5 * h * derivative_fun(t, y[i], theta, x_r, x_i);
    
    // Full step using derivative at midpoint
    y[i+1] = y[i] + h * derivative_fun(t + 0.5 * h, y_mid, theta, x_r, x_i);
    t = t + h;
  }
  
  return(y);
}

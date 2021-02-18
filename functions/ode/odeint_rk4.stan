
  /* 4th-order Runge-Kutta method
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
vector[] odeint_rk4(real t0, vector y0, data real h, data int num_steps,
    vector theta, data real[] x_r, data int[] x_i){

  int d = num_elements(y0);
  vector[d] y[num_steps+1];
  vector[d] k1; vector[d] k2; vector[d] k3; vector[d] k4;
  real t = t0;
  y[1] = y0;
  
  // Integrate at grid of time points
  for(i in 1:num_steps){
    k1 = h * derivative_fun(t          , y[i]           , theta, x_r, x_i);
    k2 = h * derivative_fun(t + 0.5 * h, y[i] + 0.5 * k1, theta, x_r, x_i);
    k3 = h * derivative_fun(t + 0.5 * h, y[i] + 0.5 * k2, theta, x_r, x_i);
    k4 = h * derivative_fun(t + h      , y[i] + k3      , theta, x_r, x_i);
    y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    t = t + h;
  }
  
  return(y);
}

  /* Forward Euler method
  *
  * Assumes that derivative_fun() is defined in the functions block before
  * this function. It must have signature
  *
  *   vector derivative_fun(real t, vector y, int[] a0, vector a1);
  *  
  * i.e. same as what can be used with ODE solvers (see
  * https://mc-stan.org/docs/2_26/stan-users-guide/coding-the-ode-system-function.html
  * ).
  *
  * Author: Juho Timonen
  *
  * @param t0 initial time
  * @param y0 initial state, D-vector
  * @param h step size (positive)
  * @param num_steps number of steps to take
  * @param a0 array of integer inputs given to derivative_fun()
  * @param theta parameter vector given to derivative_fun()
  * @return array of D-vectors, length equal to num_steps + 1
  */
vector[] odeint_euler(real t0, vector y0, data real h, data int num_steps,
    int[] a0, vector theta){

  int d = num_elements(y0);
  vector[d] y[num_steps+1];
  real t = t0;
  y[1] = y0;
  
  // Integrate at grid of time points
  for(i in 1:num_steps){
    y[i+1] = y[i] + h * derivative_fun(t, y[i], a0, theta);
    t = t + h;
  }
  return(y);
}

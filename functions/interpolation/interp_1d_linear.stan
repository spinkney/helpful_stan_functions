  /* Linear 1d interpolation
  *
  * Author: Juho Timonen
  *
  * @param y array of D-vectors, length N_in
  * @param x incresing array of reals, length N_in
  * @param x_out increasing array of reals, length N_out, values 
  * must be in (min(x), max(x)]
  * @return array of D-vectors, length N_out, corresponding to
  * interpolated values y_out
  */
vector[] interp_1d_linear(vector[] y, data real[] x, data real[] x_out){
  int left = 1;
  int right = 1;
  real w = 1.0;
  int N_in = size(x);
  int N_out = size(x_out);
  int D = size(y[1]);
  vector[D] y_out[N_out];
  for (j in 1:N_out) {
    while(x[right] < x_out[j]) {
      right = right + 1;
    }
    while(x[left+1] < x_out[j]) {
      left = left + 1;
    }
    w = (x[right] - x_out[j]) / (x[right] - x[left]);
    y_out[j] = w*y[left] + (1 - w)*y[right];
  }
  return(y_out);
}

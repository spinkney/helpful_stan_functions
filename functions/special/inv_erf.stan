   /* Inverse error function
   *
   * Copywrite Sean Pinkney, Feb. 2021
   *
   * Giles, Mike. (2012). "Approximating the Erfinv Function."
   * GPU Computing Gems Jade Edition. 10.1016/B978-0-12-385963-1.00010-1. 
   * accessed Feb. 15, 2021. at
   * https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
   *
   * @param x Real number on [0,1]
   * @param return inverse error result
   */
 real inv_erf (real x) {
    if (is_nan(x) || is_infinite(x) || x < 0 || x > 1)
      reject("inv_erf: x must be finite and between 0 and 1; ",
           "found x = ", x);
           
    real w = -log1m(square(x));
    real s = sqrt(w);
    real p;
  
   if (w < 5.0) {
    w = w - 2.5;
    p =   2.81022636e-08;
    p =   3.43273939e-07 + p * w;
    p =   -3.5233877e-06 + p * w;
    p =  -4.39150654e-06 + p * w;
    p =    0.00021858087 + p * w;
    p =   -0.00125372503 + p * w;
    p =   -0.00417768164 + p * w;
    p =      0.246640727 + p * w;
    p =       1.50140941 + p * w;
  } else {
    w = sqrt(w) - 3.000000;
    p =  -0.000200214257;
    p =   0.000100950558 + p * w;
    p =    0.00134934322 + p * w;
    p =   -0.00367342844 + p * w;
    p =    0.00573950773 + p * w;
    p =    -0.0076224613 + p * w;
    p =    0.00943887047 + p * w;
    p =       1.00167406 + p * w;
    p =       2.83297682 + p * w;
  }
   return p * x;
 }
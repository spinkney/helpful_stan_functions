functions {
  #include expect_equal.stan
  
  real my_fun(real x) {
    return 2 * x;
  }
  vector my_fun2(vector x) {
    return 2 * x;
  }
}
generated quantities {
  expect_equal_real(6.0, my_fun(3.0), 1E-8);
  expect_equal_vector([5.99999999999, 4.0, 10.0]', my_fun2([3.0, 2.0, 5.0]'), 1E-8);
}

functions {
  #include expect_equal.stan
  #include normal_copula.stan
}
data {
}
model {
}
generated quantities {
  real y = 1 - bivariate_normal_copula_cdf([0.5, 0.2]', 0.3);
  expect_equal_real(0.8663527, y, 1e-6);
}


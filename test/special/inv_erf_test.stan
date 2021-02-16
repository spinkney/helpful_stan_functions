functions {
  #inv_erf.stan
  #expect_equal.stan
}
generated quantities {
  expect_equal_real(0.3708072, inv_erf(0.4), 1E-7);
}

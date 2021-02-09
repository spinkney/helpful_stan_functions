// Copywrite Milad Kharratzadeh with modification suggested by Martin Modrak
// https://github.com/milkha/Splines_in_Stan
// Accessed Feb. 9, 2021 
// modification suggested in https://github.com/milkha/Splines_in_Stan/issues/3#issuecomment-319671070

function {
  #include build_b_spline.stan
}
data {
  int num_data;             // number of data points
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  real Y[num_data];
  real X[num_data];
  int<lower=0, upper=1> penalized_ind;
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_data] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp = append_row(rep_vector(knots[1], 
                                                                spline_degree), knots);
  vector[2 * spline_degree + num_knots] ext_knots = append_row(ext_knots_temp, 
                                                              rep_vector(knots[num_knots], 
                                                              spline_degree)); // set of extended knots

  for (ind in 1:num_basis)
    B[ind] = to_row_vector(build_b_spline(X, 
                                          to_array_1d(ext_knots), 
                                          ind, 
                                          spline_degree + 1));
  
  for(i in 1:num_data) 
    if(X[i] == ext_knots[num_knots + 2 * spline_degree - 1]) 
      B[num_basis, i] = 1;
    
}

parameters {
  row_vector[num_basis] a_raw; 
  real a0;  // intercept
  real<lower=0> sigma; 
  real<lower=0> tau;   
}

transformed parameters {
  row_vector[num_basis] a; // spline coefficients
  vector[num_data] Y_hat;
  
  if (penalized_ind == 1) {
    a[1] = a_raw[1];
    
    for (i in 2:num_basis)
      a[i] = a[i - 1] + a_raw[i] * tau; 
      
  } else a = a_raw * tau; // spline coefficients
  
  Y_hat = a0 * to_vector(X) + to_vector(a * B);
}

model {
  // Priors
  a_raw ~ normal(0, 1);
  a0 ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  //Likelihood
  Y ~ normal(Y_hat, sigma);
}

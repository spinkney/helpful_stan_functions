  real multi_student_t_Q_lpdf(matrix z, row_vector Q, matrix L, real nu) {
    int p = rows(z);
    int N = cols(z);
    real lp = N * ( lgamma( (p + nu) / 2) -
                    lgamma(nu / 2) - 0.5 * p * log(nu) - 0.5 * p * log(pi()) -
                    sum(log(diagonal(L))));

    return lp - 0.5 * sum( (nu + p) * log1p(Q / nu) );
  }
  
  real multi_skew_generalized_t_lpdf(matrix z, matrix L, vector sigma, vector lambda, real nu) {
    int p = rows(z);
    int N = cols(z);
    row_vector[N] Q = columns_dot_self(mdivide_left_tri_low(L, z));
    row_vector[N] alpha = (lambda ./ sigma)' * z;
    real lp = N * log2() + multi_student_t_Q_lpdf(z | Q, L, nu);
    
    
    for (n in 1:N)
     lp += student_t_lcdf(alpha[n] * sqrt( (p + nu) / (nu + Q[n])) | p + nu, 0, 1);
     
    return lp;
  }
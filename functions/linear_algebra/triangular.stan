  /* Lower triangular elements
   *
   * Copyright 2021, Sean Pinkney
   *
   * @param mat Matrix
   * @param K Int if mat is size N x N then \n
   K is N * N / 2
   * @param vector of lower triangular elements
   */
vector lower_ele (matrix mat, int K){
  int N = rows(mat);
  vector[K] out;
  
  int count = N;
  int pos = 1;
    
  for (k in 1:K - 1){
    count -= 1;
    out[pos:pos + count - 1] = sub_col(mat, k + 1, k, count);
    pos += count;
  }
  
  return out;
  }

  /* Lower triangular elements
   *
   * Copywrite Sean Pinkney, Feb. 10, 2021
   *
   * @param mat Matrix
   * @param K Int if mat is size N x N then \n
   K is N * N / 2
   * @param vector of lower triangular elements
   */
  vector lower_elements(matrix M, int tri_size){
    int n = rows(M);
    int counter = 1;
    vector[tri_size] out;
 
    for (i in 2:n){
      for (j in 1:(i - 1)) {
        out[counter] = M[i, j];
        counter += 1;
      }
    }
    return out;
  }

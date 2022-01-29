/** @addtogroup splines B-Spline Basis Interpolation Functions
*
* @include \interpolation\build_b_spline.stanfunctions
* 
* \ingroup interpolation
*  @{ */

/**
* B-Spline Basis interpolation
* 
* "Splines In Stan", Milad Kharratzadeh, Oct. 23, 2017,
*  https://github.com/milkha/Splines_in_Stan/blob/master/splines_in_stan.pdf
* accessed Feb. 9, 2021.
* 
* @copyright Milad Kharratzadeh
* https://github.com/milkha/Splines_in_Stan/blob/master/b_spline.stan
* Accessed Feb. 9, 2021 
*
* @param t Array of real numbers, the points at which the b_spline is calculated
* @param ext_knots Array of real numbers, the set of extended knots
* @param ind Int the index of the b_spline
* @param order Int the order of the b-spline
* @return vector of B-spline bases
*/
vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));

    if (order==1)

      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind + 1]); 

    else {

      if (ext_knots[ind] != ext_knots[ind + order - 1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind + order - 1] - ext_knots[ind]);

      if (ext_knots[ind + 1] != ext_knots[ind + order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind + 1], size(t))) / 
                 (ext_knots[ind + order] - ext_knots[ind + 1]);

      // Calculating the B-spline recursively as linear interpolation of two lower-order splines 
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order - 1) + 
                 w2 .* build_b_spline(t, ext_knots, ind + 1, order - 1);
    }

    return b_spline;
  }
/** @} */
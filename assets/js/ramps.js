"use strict";

function makeRamp(n) {
  if (!Number.isInteger(n) || n<0)
    throw new Error("Argument to makeRamp must be a non-negative integer");

  var binom = binomcoef;
  var cs = new Array(n+1);
  var sign = 1;

  for (var k=0; k<=n; k++, sign *= -1) {
    cs[k] = sign * binom( n+k , k ) * binom( 2*n+1, n-k );
  };

  return function(x) {
    if (x<0) return 0.;
    if (x>1) return 1.;

    var xk = 1.;
    var result = 0;

    for (var k=0; k<=n; k++, xk *= x) {
      result += cs[k] * xk;
    };

    // For prefactor of x^{n+1}
    result *= xk;

    return result;

  };

};

var dpq = require('./dpq');

/**
 * The distribution function of the uniform distribution.
 */

module.exports = function punif(x, a, b, lower_tail, log_p){

  if (isNaN(x) || isNaN(a) || isNaN(b)){
    return x + a + b;
  }

  if (b < a) return NaN;
  if (!isFinite(a) || !isFinite(b)) return NaN;

  if (x >= b)
    return dpq.R_DT_1(lower_tail);
 
  if (x <= a)
    return dpq.R_DT_0(lower_tail);

  if (lower_tail) {
    return dpq.R_D_val((x - a) / (b - a), log_p);
  } else {
    return dpq.R_D_val((b - x) / (b - a), log_p);
  }
};

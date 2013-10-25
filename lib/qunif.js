var dpq = require('./dpq');

/**
 *    The quantile function of the uniform distribution.
 */

module.exports = function qunif(p, a, b, lower_tail, log_p){

  if (isNaN(p) || isNaN(a) || isNaN(b)){
    return p + a + b;
  }

  dpq.R_Q_P01_check(p, log_p);
  if (!isFinite(a) || !isFinite(b)) return NaN;
  if (b < a) return NaN;
  if (b == a) return a;

  return a + dpq.R_DT_qIv(p, lower_tail, log_p) * (b - a);
};





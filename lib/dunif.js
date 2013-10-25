var dpq = require('./dpq');

/**
 * The density of the uniform distribution.
 */

module.exports = function dunif(x, a, b, give_log){
  if (isNaN(x) || isNaN(a) || isNaN(b)){
    return x + a + b;
  }
  
  if (b <= a) return NaN;

  if (a <= x && x <= b){
    return give_log ? -Math.log(b - a) : 1. / (b - a);
  }

  return dpq.R_D__0(give_log);
};

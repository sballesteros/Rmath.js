var cst = require('./cst')
  , dpq = require('./dpq');

/**
 * Compute the density of the normal distribution.
 */
module.exports = function dnorm(x, mu, sigma, give_log){

  if(!isFinite(sigma)) return dpq.R_D__0(give_log);
  if(!isFinite(x) && mu == x) return NaN;/* x-mu is NaN */
  if (sigma <= 0) {
    if (sigma < 0) return NaN;
    /* sigma == 0 */
    return (x == mu) ? Number.POSITIVE_INFINITY : dpq.R_D__0(give_log);
  }
  x = (x - mu) / sigma;

  if(!isFinite(x)) return dpq.R_D__0(give_log);
  return (give_log ?
	  -(cst.M_LN_SQRT_2PI  +	0.5 * x * x + Math.log(sigma)) :
	  cst.M_1_SQRT_2PI * Math.exp(-0.5 * x * x)  /	  sigma);
  /* M_1_SQRT_2PI = 1 / sqrt(2 * pi) */
};

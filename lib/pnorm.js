var cst = require('./cst')
  , dpq = require('./dpq');

/** 
 * The main computation evaluates near-minimax approximations derived
 * from those in "Rational Chebyshev approximations for the error
 * function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
 * transportable program uses rational functions that theoretically
 * approximate the normal distribution function to at least 18
 * significant decimal digits.  The accuracy achieved depends on the
 * arithmetic system, the compiler, the intrinsic functions, and
 * proper selection of the machine-dependent constants.
 *
 *  REFERENCE
 *
 *	Cody, W. D. (1993).
 *	ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
 *	Special Function Routines and Test Drivers".
 *	ACM Transactions on Mathematical Software. 19, 22-32.
 *
 *  EXTENSIONS
 *
 *  The "_both" , lower, upper, and log_p  variants were added by
 *  Martin Maechler, Jan.2000;
 *  as well as log1p() and similar improvements later on.
 *
 *  James M. Rath contributed bug report PR#699 and patches correcting SIXTEN
 *  and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line
 *  with the original Cody code.
 */


module.exports = function pnorm(x, mu, sigma, lower_tail, log_p){

  var res ={
    p: undefined,
    cp: undefined
  };

  /* Note: The structure of these checks has been carefully thought through.
   * For example, if x == mu and sigma == 0, we get the correct answer 1.
   */
  if(isNaN(x) || isNaN(mu) || isNaN(sigma)){
    return x + mu + sigma;
  }

  if(!isFinite(x) && mu == x) return NaN;/* x-mu is NaN */
  if (sigma <= 0) {
    if(sigma < 0) return NaN;
    /* sigma = 0 : */
    return (x < mu) ? dpq.R_DT_0(lower_tail) : dpq.R_DT_1(lower_tail);
  }
  res.p = (x - mu) / sigma;
  if(!isFinite(res.p))
    return (x < mu) ? dpq.R_DT_0(lower_tail) : dpq.R_DT_1(lower_tail);
  x = res.p;

  res = pnorm_both(x, res, (lower_tail ? 0 : 1), log_p);

  return(lower_tail ? res.p : res.cp);
};




var SIXTEN = 16; /* Cutoff allowing exact "*" and "/" */

var a = [
  2.2352520354606839287,
  161.02823106855587881,
  1067.6894854603709582,
  18154.981253343561249,
  0.065682337918207449113
];

var b = [
  47.20258190468824187,
  976.09855173777669322,
  10260.932208618978205,
  45507.789335026729956
];

var c = [
  0.39894151208813466764,
  8.8831497943883759412,
  93.506656132177855979,
  597.27027639480026226,
  2494.5375852903726711,
  6848.1904505362823326,
  11602.651437647350124,
  9842.7148383839780218,
  1.0765576773720192317e-8
];

var d = [
  22.266688044328115691,
  235.38790178262499861,
  1519.377599407554805,
  6485.558298266760755,
  18615.571640885098091,
  34900.952721145977266,
  38912.003286093271411,
  19685.429676859990727
];

var p = [
  0.21589853405795699,
  0.1274011611602473639,
  0.022235277870649807,
  0.001421619193227893466,
  2.9112874951168792e-5,
  0.02307344176494017303
];

var q = [
  1.28426009614491121,
  0.468238212480865118,
  0.0659881378689285515,
  0.00378239633202758244,
  7.29751555083966205e-5
];

/* Difference between 1.0 and the minimum double greater than 1.0 */
var DBL_EPSILON = 2.2204460492503131e-16; //TODO FIX;


function pnorm_both(x, res, i_tail, log_p)
{
  /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
     if(lower) return  res.p := P[X <= x]
     if(upper) return res.cp := P[X >  x] = 1 - P[X <= x]
  */

  var xden, xnum, temp, del, eps, xsq, y;
  var i, lower, upper;


  if(isNaN(x)) { res.p = res.cp = x; return res; }


  /* Consider changing these : */
  eps = DBL_EPSILON * 0.5;

  /* i_tail in {0,1,2} =^= {lower, upper, both} */
  lower = i_tail != 1;
  upper = i_tail != 0;

  y = Math.abs(x);
  if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
    if (y > eps) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i) {
	xnum = (xnum + a[i]) * xsq;
	xden = (xden + b[i]) * xsq;
      }
    } else xnum = xden = 0.0;

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if(lower)  res.p = 0.5 + temp;
    if(upper) res.cp = 0.5 - temp;
    if(log_p) {
      if(lower)  res.p = Math.log(res.p);
      if(upper) res.cp = Math.log(res.cp);
    }
  }
  else if (y <= cst.M_SQRT_32) {

    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

    xsq = trunc(y * SIXTEN) / SIXTEN;				
    del = (y - xsq) * (y + xsq);
    if(log_p) {	
      res.p = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
      if((lower && x > 0.) || (upper && x <= 0.))
	res.cp = Math.log1p(-Math.exp(-xsq * xsq * 0.5) *
		       Math.exp(-del * 0.5) * temp);
    }								
    else {	
      res.p = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
      res.cp = 1.0 - res.p;
    }

    if (x > 0.) {/* swap  res.cp <--> res.p */
      temp = res.p; 
      if(lower) {res.p = res.cp;}
      res.cp = temp;
    }

  }

  /* else	  |x| > sqrt(32) = 5.657 :
   * the next two case differentiations were really for lower=T, log=F
   * Particularly	 *not*	for  log_p !

   * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
   *
   * Note that we do want symmetry(0), lower/upper -> hence use y
   */
  else if((log_p && y < 1e170) /* avoid underflow below */
	  /*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	   * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

	   xsq = x*x;

	   if(xsq * DBL_EPSILON < 1.)
	   del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	   else
	   del = 0.;
	   res.p = -.5*xsq - cst.M_LN_SQRT_2PI - log(x) + log1p(-del);
	   res.cp = log1p(-exp(res.p)); /.* ~ log(1) = 0 *./

 	   swap_tail;

	   [Yes, but xsq might be infinite.]

	  */
	  || (lower && -37.5193 < x  &&  x < 8.2924)
	  || (upper && -8.2924  < x  &&  x < 37.5193)
	 ) {

    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
    xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (cst.M_1_SQRT_2PI - temp) / y;


    xsq = trunc(x * SIXTEN) / SIXTEN;				
    del = (x - xsq) * (x + xsq);
    if(log_p) {	
      res.p = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
      if((lower && x > 0.) || (upper && x <= 0.))
	res.cp = Math.log1p(-Math.exp(-xsq * xsq * 0.5) *
		       Math.exp(-del * 0.5) * temp);
    }								
    else {	
      res.p = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
      res.cp = 1.0 - res.p;
    }

    if (x > 0.) {/* swap  res.cp <--> res.p */
      temp = res.p; 
      if(lower) {res.p = res.cp;}
      res.cp = temp;
    }


  } else { /* large x such that probs are 0 or 1 */
    if(x > 0) {	res.p = dpq.R_D__1(log_p); res.cp = dpq.R_D__0(log_p);	}
    else {	        res.p = dpq.R_D__0(log_p); res.cp = dpq.R_D__1(log_p);	}
  }


  //???
  //#ifdef NO_DENORMS
  //    /* do not return "denormalized" -- we do in R */
  //    if(log_p) {
  //	if(res.p > -Number.MIN_VALUE)	 res.p = -0.;
  //	if(res.cp > -Number.MIN_VALUE)res.cp = -0.;
  //    }
  //    else {
  //	if(res.p < Number.MIN_VALUE)	 res.p = 0.;
  //	if(res.cp < Number.MIN_VALUE)	res.cp = 0.;
  //    }
  //#endif
  return res;
}

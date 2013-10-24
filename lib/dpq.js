/* Utilities for "dpq" handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */

function R_D__0 (log_p){ /* 0 */
  return log_p ? Number.NEGATIVE_INFINITY : 0.0; 
};
exports.R_D__0 = R_D__0;

function R_D__1(log_p){ /* 1 */
  return log_p ? 0. : 1.;
};
exports.R_D__1 = R_D__1;

function R_DT_0(lower_tail){/* 0 */
  if(lower_tail){
    return R_D__0;
  } else {
    return R_D__1;
  }
};
exports.R_DT_0 = R_DT_0;

function R_DT_1(lower_tail){ /* 1 */
  if(lower_tail){
    return R_D__1;
  } else {
    return R_D__0;
  }
};

exports.R_DT_1 = R_DT_1;


/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy */
function R_D_Lval(p, lower_tail){ return (lower_tail ? (p) : (0.5 - (p) + 0.5));}; /*  p  */
exports.R_D_Lval = R_D_Lval;

function R_D_Cval(p, lower_tail){ return (lower_tail ? (0.5 - (p) + 0.5) : (p));}; /*  1 - p */
exports.R_D_Cval = R_D_Cval;

function R_D_val(x, log_p) {return (log_p ? Math.log(x) : (x));}; /*  x  in pF(x,..) */
exports.R_D_val = R_D_val;

function R_D_qIv(p, log_p) {return (log_p ? Math.exp(p) : (p));}; /*  p  in qF(p,..) */
exports.R_D_qIv = R_D_qIv;

function R_D_exp(x, log_p) {return (log_p ? (x): Math.exp(x));};  /* exp(x) */
exports.R_D_exp = R_D_exp;

function R_D_log(p, log_p) {return (log_p ? (p): Math.log(p));};  /* log(p) */
exports.R_D_log = R_D_log;

function R_D_Clog(p, log_p){return (log_p ? Math.log1p(-(p)) : (0.5 - (p) + 0.5));}; /* [log](1-p) */
exports.R_D_Clog = R_D_Clog;



/* log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x))) : */
function R_Log1_Exp(x){ return ((x) > -Math.LN2 ? Math.log(-Math.expm1(x)) : Math.log1p(-Math.exp(x)));};
exports.R_Log1_Exp = R_Log1_Exp;

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
function R_D_LExp(x, log_p){ return (log_p ? R_Log1_Exp(x) : Math.log1p(-x))};
exports.R_D_LExp = R_D_LExp;

function R_DT_val(x, lower_tail, log_p) {(lower_tail ? R_D_val(x, log_p)  : R_D_Clog(x, log_p))};
exports.R_DT_val = R_DT_val;

function R_DT_Cval(x, lower_tail, log_p) {return (lower_tail ? R_D_Clog(x, log_p) : R_D_val(x, log_p));};
exports.R_DT_Cval = R_DT_Cval;

/*function R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
function R_DT_qIv(p, lower_tail, log_p){return (log_p ? (lower_tail ? Math.exp(p) : - Math.expm1(p)) 
			                        : R_D_Lval(p, lower_tail)); };
exports.R_DT_qIv = R_DT_qIv;

/*function R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
function R_DT_CIv(p, lower_tail, log_p)	{return (log_p ? (lower_tail ? -Math.expm1(p) : Math.exp(p)) 
			                         : R_D_Cval(p, lower_tail));};
exports.R_DT_CIv = R_DT_CIv;

function R_DT_exp(x, lower_tail, log_p) {return R_D_exp(R_D_Lval(x, lower_tail), log_p);}; /* exp(x) */
exports.R_DT_exp = R_DT_exp;

function R_DT_Cexp(x, lower_tail, log_p) {return R_D_exp(R_D_Cval(x, lower_tail), log_p);}; /* exp(1 - x) */
exports.R_DT_Cexp = R_DT_Cexp;

function R_DT_log(p, lower_tail, log_p) {return	(lower_tail? R_D_log(p, log_p) : R_D_LExp(p, log_p));}; /* log(p) in qF */
exports.R_DT_log = R_DT_log;

function R_DT_Clog(p, lower_tail, log_p) {return (lower_tail? R_D_LExp(p, log_p): R_D_log(p, log_p));} /* log(1-p) in qF*/
exports.R_DT_Clog = R_DT_Clog;

function R_DT_Log(p, lower_tail) {return (lower_tail? (p) : R_Log1_Exp(p));}; /* ==   R_DT_log when we already "know" log_p == TRUE :*/
exports.R_DT_Log = R_DT_Log;


function R_Q_P01_check(p, log_p){  
  if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)) ){
    return NaN;
  }
};
exports.R_Q_P01_check = R_Q_P01_check;

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
function R_Q_P01_boundaries(p, _LEFT_, _RIGHT_, lower_tail, log_p){
  if (log_p) {					
    if(p > 0)					
      return NaN;
    if(p == 0) /* upper bound*/
      return lower_tail ? _RIGHT_ : _LEFT_;
    if(p == Number.NEGATIVE_INFINITY)	
      return lower_tail ? _LEFT_ : _RIGHT_;
  }		
  else { /* !log_p */
    if(p < 0 || p > 1)
      return NaN;
    if(p == 0)
      return lower_tail ? _LEFT_ : _RIGHT_;
    if(p == 1)
      return lower_tail ? _RIGHT_ : _LEFT_;
  }
};
exports.R_Q_P01_boundaries = R_Q_P01_boundaries;


function R_P_bounds_01(x, x_min, x_max, lower_tail){
  if(x <= x_min) return R_DT_0(lower_tail);
  if(x >= x_max) return R_DT_1(lower_tail);
};
exports.R_P_bounds_01 = R_P_bounds_01;

/* is typically not quite optimal for (-Inf,Inf) where
 * you'd rather have */

function R_P_bounds_Inf_01(x, lower_tail){
  if(!isFinite(x)) {
    if (x > 0) return R_DT_1(lower_tail);
    /* x < 0 */return R_DT_0(lower_tail);
  }
};
exports.R_P_bounds_Inf_01 = R_P_bounds_Inf_01;

/* additions for density functions (C.Loader) */
function R_D_fexp(f, x, give_log) {return (give_log ? -0.5*Math.log(f)+(x) : Math.exp(x)/Math.sqrt(f));};
exports.R_D_fexp = R_D_fexp;

function R_D_forceint(x) {return Math.floor((x) + 0.5);};
exports.R_D_forceint = R_D_forceint;

function R_D_nonint(x) {return (Math.abs((x) - Math.floor((x)+0.5)) > 1e-7);};
exports.R_D_nonint = R_D_nonint;

/* [neg]ative or [non int]eger : */
function R_D_negInonint(x) {return (x < 0. || R_D_nonint(x));};
exports.R_D_negInonint = R_D_negInonint;

function R_D_nonint_check(x, log_p) {
  if(R_D_nonint(x)) {
    return R_D__0(log_p);
  } else {
    return x;
  }
}
exports.R_D_nonint_check = R_D_nonint_check;

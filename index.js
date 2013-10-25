require('es6-shim'); //For the Math object (Math.log1p and Math.expm1)

module.exports = {
  dnorm: require('./lib/dnorm'),
  qnorm: require('./lib/qnorm'),
  pnorm: require('./lib/pnorm'),

  dunif: require('./lib/dunif'),
  qunif: require('./lib/qunif'),
  punif: require('./lib/punif')
};

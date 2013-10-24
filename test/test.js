var rstat = require('..')
  , assert = require('assert');

describe('rmath', function(){

  describe('dnorm', function(){
    it('should compare with R', function(){
      assert.equal(rstat.dnorm(2, 1, 2, false).toFixed(7), 0.1760327);
    });
    
    it('should compare with R in log scale', function(){
      assert.equal(rstat.dnorm(2, 1, 2, true).toFixed(6), -1.737086);
    });
  });

  describe('pnorm', function(){
    it('should compare with R', function(){
      assert.equal(rstat.pnorm(2, 1, 2, true).toFixed(7), 0.6914625);
    });

    it('should compare with R with lower_tail == false', function(){
      assert.equal(rstat.pnorm(2, 1, 2, false).toFixed(7), 0.3085375);
    });

    it('should compare with R with log_p == true', function(){
      assert.equal(rstat.pnorm(2, 1, 2, false, true).toFixed(6), -1.175912);
    });
  });

  describe('qnorm', function(){
    it('should compare with R', function(){
      assert.equal(rstat.qnorm(0.5, 1, 2, true), 1);
    });

    it('should compare with R', function(){
      assert.equal(rstat.qnorm(0.1, 1, 2, true).toFixed(6), -1.563103);
    });

    it('should compare with R with lower_tail == false', function(){
      assert.equal(rstat.qnorm(0.5, 1, 2, false), 1);
    });

    it('should compare with R with log_p == true', function(){
      assert(isNaN(rstat.qnorm(0.5, 1, 2, false, true)));
    });
  });

});

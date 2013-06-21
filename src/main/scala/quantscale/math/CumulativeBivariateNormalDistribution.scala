package quantscale.math

import quantscale.analytic.CodyCND

trait CumulativeBivariateNormalDistribution {
  def value(x: Double, y: Double, corr: Double, cnd: Function1D = CodyCND): Double
}

/**
 * based on West adaptation of Genz code
 */
object GenzWestLeFlochBVD extends CumulativeBivariateNormalDistribution {
  val XX = Array(
    Array(-0.932469514203152, -0.981560634246719, -0.993128599185095),
    Array(-0.661209386466265, -0.904117256370475, -0.963971927277914),
    Array(-0.238619186083197, -0.769902674194305, -0.912234428251326),
    Array(0., -0.587317954286617, -0.839116971822219),
    Array(0., -0.36783149899818, -0.746331906460151),
    Array(0., -0.125233408511469, -0.636053680726515),
    Array(0., 0., -0.510867001950827),
    Array(0., 0., -0.37370608871542),
    Array(0., 0., -0.227785851141645),
    Array(0., 0., -0.0765265211334973))
  val W = Array(
    Array(0.17132449237917, 0.0471753363865118, 0.0176140071391521),
    Array(0.360761573048138, 0.106939325995318, 0.0406014298003869),
    Array(0.46791393457269, 0.160078328543346, 0.0626720483341091),
    Array(0., 0.203167426723066, 0.0832767415767048),
    Array(0., 0.233492536538355, 0.10193011981724),
    Array(0., 0.249147045813403, 0.118194531961518),
    Array(0., 0., 0.131688638449177),
    Array(0., 0., 0.142096109318382),
    Array(0., 0., 0.149172986472604),
    Array(0., 0., 0.152753387130726))

  def value(xi: Double, yi: Double, corr: Double, cndi: Function1D): Double = {
    var res = 0.0
    var NG = 0
    var LG = 0
    var cnd = if (cndi == null) CodyCND else cndi

    var x = if (xi == Double.PositiveInfinity || xi == Double.NegativeInfinity) math.signum(xi) * 1.0e100 else xi
    var y = if (yi == Double.NegativeInfinity || yi == Double.PositiveInfinity) math.signum(yi) * 1.0e100 else yi

    if (math.abs(corr) < 0.3) {
      NG = 1; LG = 3
    } else if (math.abs(corr) < 0.75) {
      NG = 2; LG = 6
    } else {
      NG = 3; LG = 10
    }

    var h = -x
    var k = -y
    var hk = h * k
    if (math.abs(corr) < 0.925) {
      if (math.abs(corr) > 0.0) {
        var hs = (h * h + k * k) / 2.0;
        var asr = math.asin(corr)
        var i = 1
        while (i <= LG) {
          var iss = -1
          while (iss <= 1) {
            var sn = math.sin(asr * (iss * XX(i - 1)(NG - 1) + 1) * 0.5)
            res = res + W(i - 1)(NG - 1) * math.exp((sn * hk - hs) / (1.0 - sn * sn))
            iss += 2
          }
          i += 1
        }
        res = res * asr * 0.795774715459476678e-1;
      }
      res = res + cnd.value(-h) * cnd.value(-k);
    } else {
      if (corr < 0.0) {
        k *= -1.0; hk *= -1.0;
      }
      if (Math.abs(corr) < 1.0) {
        val Ass = (1.0 - corr) * (1.0 + corr);
        var a = Math.sqrt(Ass);
        val bs = (h - k) * (h - k);
        val c = (4.0 - hk) / 8.0;
        val d = (12.0 - hk) / 16.0;
        var asr = -(bs / Ass + hk) / 2.0;
        if (asr > -100.0) {
          res = a * Math.exp(asr) *
            (1.0 - c * (bs - Ass) * (1.0 - d * bs / 5.0) / 3.0 +
              c * d * Ass * Ass / 5.0);
        }
        if (-hk < 100.0) {
          val B = Math.sqrt(bs);
          res = res -
            Math.exp(-hk / 2.0) * 2.506628274631 *
            cnd.value(-B / a) * B *
            (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
        }
        a /= 2.0;
        var i = 1
        while (i <= LG) {
          var iss = -1

          while (iss <= 1) {
            var xs = a * (iss * XX(i - 1)(NG - 1) + 1.0)
            xs = math.abs(xs * xs)
            val rs = math.sqrt(1 - xs)
            asr = -(bs / xs + hk) / 2.0
            if (asr > -100.0) {
              res = res + a * W(i - 1)(NG - 1) * math.exp(asr) *
                (math.exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) /
                  rs - (1.0 + c * xs * (1.0 + d * xs)));
            }
            iss += 2
          }
          i += 1
        }
        res *= -0.159154943091895336;
      }
      if (corr > 0.0) {
        res += cnd.value(-math.max(h, k));
      } else {
        res *= -1.0;
        if (k > h) {
          if (h >= 0) {
            res += (cnd.value(-h) - cnd.value(-k));
          } else {
            res += (cnd.value(k) - cnd.value(h));
          }
        }
      }
    }
    return res
  }
}
object GenzWestBVD extends CumulativeBivariateNormalDistribution {
  val XX = Array(
    Array(-0.932469514203152, -0.981560634246719, -0.993128599185095),
    Array(-0.661209386466265, -0.904117256370475, -0.963971927277914),
    Array(-0.238619186083197, -0.769902674194305, -0.912234428251326),
    Array(0., -0.587317954286617, -0.839116971822219),
    Array(0., -0.36783149899818, -0.746331906460151),
    Array(0., -0.125233408511469, -0.636053680726515),
    Array(0., 0., -0.510867001950827),
    Array(0., 0., -0.37370608871542),
    Array(0., 0., -0.227785851141645),
    Array(0., 0., -0.0765265211334973))
  val W = Array(
    Array(0.17132449237917, 0.0471753363865118, 0.0176140071391521),
    Array(0.360761573048138, 0.106939325995318, 0.0406014298003869),
    Array(0.46791393457269, 0.160078328543346, 0.0626720483341091),
    Array(0., 0.203167426723066, 0.0832767415767048),
    Array(0., 0.233492536538355, 0.10193011981724),
    Array(0., 0.249147045813403, 0.118194531961518),
    Array(0., 0., 0.131688638449177),
    Array(0., 0., 0.142096109318382),
    Array(0., 0., 0.149172986472604),
    Array(0., 0., 0.152753387130726))

  def value(xi: Double, yi: Double, corr: Double, cndi: Function1D): Double = {
//    println("bivariate "+xi+" "+yi+" "+corr)
    var res = 0.0
    var NG = 0
    var LG = 0
    var cnd = if (cndi == null) CodyCND else cndi

    var x = if (xi == Double.PositiveInfinity || xi == Double.NegativeInfinity) math.signum(xi) * 1.0e100 else xi
    var y = if (yi == Double.NegativeInfinity || yi == Double.PositiveInfinity) math.signum(yi) * 1.0e100 else yi

    if (math.abs(corr) < 0.3) {
      NG = 1; LG = 3
    } else if (math.abs(corr) < 0.75) {
      NG = 2; LG = 6
    } else {
      NG = 3; LG = 10
    }

    var h = -x
    var k = -y
    var hk = h * k
    if (math.abs(corr) < 0.925) {
      if (math.abs(corr) > 0.0) {
        var hs = (h * h + k * k) / 2.0;
        var asr = math.asin(corr)
        var i = 1
        while (i <= LG) {
          var iss = -1
          while (iss <= 1) {
            var sn = math.sin(asr * (iss * XX(i - 1)(NG - 1) + 1) * 0.5)
            res = res + W(i - 1)(NG - 1) * math.exp((sn * hk - hs) / (1.0 - sn * sn))
            iss += 2
          }
          i += 1
        }
        res = res * asr * 0.795774715459476678e-1;
      }
      res = res + cnd.value(-h) * cnd.value(-k);
    } else {
      if (corr < 0.0) {
        k *= -1.0; hk *= -1.0;
      }
      if (Math.abs(corr) < 1.0) {
        val Ass = (1.0 - corr) * (1.0 + corr);
        var a = Math.sqrt(Ass);
        val bs = (h - k) * (h - k);
        val c = (4.0 - hk) / 8.0;
        val d = (12.0 - hk) / 16.0;
        var asr = -(bs / Ass + hk) / 2.0;
        if (asr > -100.0) {
          res = a * Math.exp(asr) *
            (1.0 - c * (bs - Ass) * (1.0 - d * bs / 5.0) / 3.0 +
              c * d * Ass * Ass / 5.0);
        }
        if (-hk < 100.0) {
          val B = Math.sqrt(bs);
          res = res -
            Math.exp(-hk / 2.0) * 2.506628274631 *
            cnd.value(-B / a) * B *
            (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
        }
        a /= 2.0;
        var i = 1
        while (i <= LG) {
          var iss = -1

          while (iss <= 1) {
            var xs = a * (iss * XX(i - 1)(NG - 1) + 1.0)
            xs = math.abs(xs * xs)
            val rs = math.sqrt(1 - xs)
            asr = -(bs / xs + hk) / 2.0
            if (asr > -100.0) {
              res = res + a * W(i - 1)(NG - 1) * math.exp(asr) *
                (math.exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) /
                  rs - (1.0 + c * xs * (1.0 + d * xs)));
            }
            iss += 2
          }
          i += 1
        }
        res *= -0.159154943091895336;
      }
      if (corr > 0.0) {
        res += cnd.value(-math.max(h, k));
      } else {
        res *= -1.0;
        if (k > h) {
//          println("bivariate KH "+k+" "+h+" "+corr)
          res += (cnd.value(k) - cnd.value(h));
        }
      }
    }
    return res
  }
}


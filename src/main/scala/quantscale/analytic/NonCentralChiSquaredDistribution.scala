package quantscale.analytic

import org.apache.commons.math3.special.Gamma

import quantscale.fdm.Epsilon

//TODO compare perf with opengamma version
object NonCentralChiSquaredDistribution {

  def cdfFraser(x: Double, degrees: Double, nonCentrality: Double): Double = {
    val _dofOverTwo = degrees / 2.0;
    val _lambdaOverTwo = nonCentrality / 2.0;
    val s = math.sqrt(_lambdaOverTwo * 2.0);
    val mu = math.sqrt(x);
    var z: Double = 0.0;
    if (java.lang.Double.doubleToLongBits(mu) == java.lang.Double.doubleToLongBits(s)) {
      z = (1 - _dofOverTwo * 2.0) / 2 / s;
    } else {
      z = mu - s - (_dofOverTwo * 2.0 - 1) / 2 * (math.log(mu) - math.log(s)) / (mu - s);
    }
    return CumulativeNormalDistribution.value(z)
  }

  def cdfOpenGamma(x: Double, degrees: Double, nonCentrality: Double, errtol: Double = 1e-8): Double = {
    if (x < 0) {
      return 0.0;
    }
    val _dofOverTwo = degrees / 2.0;
    val _lambdaOverTwo = nonCentrality / 2.0;
    if ((_dofOverTwo + _lambdaOverTwo) > 1000) {
      return cdfFraser(x, degrees, nonCentrality);
    }
    val _k: Int = math.round(_lambdaOverTwo).toInt
    var _pStart: Double = 0.0
    if (_lambdaOverTwo == 0) {
      _pStart = 0.0;
    } else {
      val logP = -_lambdaOverTwo + _k * math.log(_lambdaOverTwo) - Gamma.logGamma(_k + 1);
      _pStart = math.exp(logP);
    }
    var regGammaStart = 0.0;
    val halfX = x / 2.0;
    val logX = math.log(halfX);
    regGammaStart = Gamma.regularizedGammaP(_dofOverTwo + _k, halfX);

    var sum = _pStart * regGammaStart;
    var oldSum = Double.NegativeInfinity;
    var p = _pStart;
    var regGamma = regGammaStart;
    var temp = 0.0;
    var i = _k;

    // first add terms below _k
    while (i > 0 && math.abs(sum - oldSum) / sum > errtol) {
      i -= 1;
      p *= (i + 1) / _lambdaOverTwo;
      temp = (_dofOverTwo + i) * logX - halfX - Gamma.logGamma(_dofOverTwo + i + 1);
      regGamma += math.exp(temp);
      oldSum = sum;
      sum += p * regGamma;
    }

    p = _pStart;
    regGamma = regGammaStart;
    oldSum = Double.NegativeInfinity;
    i = _k;
    while (math.abs(sum - oldSum) / sum > errtol) {
      i+=1;
      p *= _lambdaOverTwo / i;
      temp = (_dofOverTwo + i - 1) * logX - halfX - Gamma.logGamma(_dofOverTwo + i);
      regGamma -= math.exp(temp);
      oldSum = sum;
      sum += p * regGamma;
    }

    return sum;
  }

  def cdfDing(x: Double, f: Double, theta: Double, max_iter: Int = 5000, errtol: Double = 1e-8): Double = {
    //
    // This is an implementation of:
    //
    // Algorithm AS 275:
    // Computing the Non-Central #2 Distribution Function
    // Cherng G. Ding
    // Applied Statistics, Vol. 41, No. 2. (1992), pp. 478-482.
    //
    // This uses a stable forward iteration to sum the
    // CDF, unfortunately this can not be used for large
    // values of the non-centrality parameter because:
    // * The first term may underfow to zero.
    // * We may need an extra-ordinary number of terms
    //   before we reach the first *significant* term.
    //
    // Special case:
    if (x < Epsilon.MACHINE_EPSILON)
      return 0.0;
    var tk = GammaFunction.gammaPDerivative(f / 2 + 1, x / 2);
    val lambda = theta / 2;
    var vk = math.exp(-lambda);
    var uk = vk;
    var sum = tk * vk;
    if (sum == 0)
      return sum;
    var lterm, term = 0.0
    var i = 1
    var isFinished = false
    while (i <= max_iter && !isFinished) {
      tk = tk * x / (f + 2 * i);
      uk = uk * lambda / i;
      vk = vk + uk;
      lterm = term;
      term = vk * tk;
      sum += term;
      if ((math.abs(term / sum) < errtol) && (term <= lterm)) {
        isFinished = true
      }
      i += 1
    }
    //Error check:
    if (i >= max_iter)
      throw new RuntimeException(
        "Series did not converge, closest value was " + sum);
    return sum;
  }

  def cdfBentonKrishnamoorthy(y: Double, n: Double, lambda: Double, max_iter: Int = 5000, errtol: Double = 1e-8): Double = {
    //
    // This is taken more or less directly from:
    //
    // Computing discrete mixtures of continuous
    // distributions: noncentral chisquare, noncentral t
    // and the distribution of the square of the sample
    // multiple correlation coeficient.
    // D. Benton, K. Krishnamoorthy.
    // Computational Statistics & Data Analysis 43 (2003) 249 - 267
    //
    // We're summing a Poisson weighting term multiplied by
    // a central chi squared distribution.
    //
    // Special case:
    if (y < Epsilon.MACHINE_EPSILON)
      return 0;
    var errorf, errorb = 0.0;

    val x = y / 2;
    val del = lambda / 2;
    //
    // Starting location for the iteration, we'll iterate
    // both forwards and backwards from this point.  The
    // location chosen is the maximum of the Poisson weight
    // function, which ocurrs *after* the largest term in the
    // sum.
    //
    val k: Int = math.round(del).toInt
    val a = n / 2 + k;
    // Central chi squared term for forward iteration:
    var gamkf = Gamma.regularizedGammaP(a, x); //gammaP?

    if (lambda == 0)
      return gamkf;
    // Central chi squared term for backward iteration:
    var gamkb = gamkf;
    // Forwards Poisson weight:
    var poiskf = GammaFunction.gammaPDerivative(k + 1, del);
    // Backwards Poisson weight:
    var poiskb = poiskf;
    // Forwards gamma function recursion term:
    var xtermf = GammaFunction.gammaPDerivative(a, x);
    // Backwards gamma function recursion term:
    var xtermb = xtermf * x / a;
    var sum = poiskf * gamkf;
    if (sum == 0)
      return sum;
    var i = 1;
    //
    // Backwards recursion first, this is the stable
    // direction for gamma function recurrences:
    //
    var isBreak = false
    while (i <= k && !isBreak) {
      xtermb *= (a - i + 1) / x;
      gamkb += xtermb;
      poiskb = poiskb * (k - i + 1) / del;
      errorf = errorb;
      errorb = gamkb * poiskb;
      sum += errorb;
      if ((math.abs(errorb / sum) < errtol) && (errorb <= errorf))
        isBreak = true
      i += 1;
    }
    i = 1;
    //
    // Now forwards recursion, the gamma function
    // recurrence relation is unstable in this direction,
    // so we rely on the magnitude of successive terms
    // decreasing faster than we introduce cancellation error.
    // For this reason it's vital that k is chosen to be *after*
    // the largest term, so that successive forward iterations
    // are strictly (and rapidly) converging.
    //
    do {
      xtermf = xtermf * x / (a + i - 1);
      gamkf = gamkf - xtermf;
      poiskf = poiskf * del / (k + i);
      errorf = poiskf * gamkf;
      sum += errorf;
      i += 1;
    } while ((math.abs(errorf / sum) > errtol) && (i < max_iter));

    //Error check:
    if (i >= max_iter)
      throw new RuntimeException(
        "Series did not converge, closest value was " + sum);

    return sum;
  }

  def cdfQBentonKrishnamoorthy(x: Double, f: Double, theta: Double, max_iter: Int = 5000, errtol: Double = 1e-8): Double = {
    //
    // Computes the complement of the Non-Central Chi-Square
    // Distribution CDF by summing a weighted sum of complements
    // of the central-distributions.  The weighting factor is
    // a Poisson Distribution.
    //
    // This is an application of the technique described in:
    //
    // Computing discrete mixtures of continuous
    // distributions: noncentral chisquare, noncentral t
    // and the distribution of the square of the sample
    // multiple correlation coeficient.
    // D. Benton, K. Krishnamoorthy.
    // Computational Statistics & Data Analysis 43 (2003) 249 - 267
    //

    // Special case:
    if (x == 0)
      return 1;

    //
    // Initialize the variables we'll be using:
    //
    val lambda = theta / 2;
    val del = f / 2;
    val y = x / 2;
    var sum = 0.0;
    //
    // k is the starting location for iteration, we'll
    // move both forwards and backwards from this point.
    // k is chosen as the peek of the Poisson weights, which
    // will occur *before* the largest term.
    //
    val k: Int = math.round(lambda).toInt
    // Forwards and backwards Poisson weights:
    var poisf = GammaFunction.gammaPDerivative(1 + k, lambda);
    var poisb = poisf * k / lambda;
    // Initial forwards central chi squared term:
    var gamf = Gamma.regularizedGammaQ(del + k, y);
    // Forwards and backwards recursion terms on the central chi squared:
    var xtermf = GammaFunction.gammaPDerivative(del + 1 + k, y);
    var xtermb = xtermf * (del + k) / y;
    // Initial backwards central chi squared term:
    var gamb = gamf - xtermb;

    //
    // Forwards iteration first, this is the
    // stable direction for the gamma function
    // recurrences:
    //
    var i = k;
    var isBreak = false
    while ((i - k) < max_iter && !isBreak) {
      val term = poisf * gamf;
      sum += term;
      poisf *= lambda / (i + 1);
      gamf += xtermf;
      xtermf *= y / (del + i + 1);
      if (((sum == 0) || (math.abs(term / sum) < errtol)) && (term >= poisf * gamf))
        isBreak = true;
      i += 1
    }
    //Error check:
    if ((i - k) >= max_iter)
      throw new RuntimeException(
        "Series did not converge, closest value was " + sum);
    //
    // Now backwards iteration: the gamma
    // function recurrences are unstable in this
    // direction, we rely on the terms deminishing in size
    // faster than we introduce cancellation errors.
    // For this reason it's very important that we start
    // *before* the largest term so that backwards iteration
    // is strictly converging.
    //
    i = k - 1
    isBreak = false
    while (i >= 0 && !isBreak) {
      val term = poisb * gamb;
      sum += term;
      poisb *= i / lambda;
      xtermb *= (del + i) / y;
      gamb -= xtermb;
      if ((sum == 0) || (math.abs(term / sum) < errtol))
        isBreak = true
      i -= 1
    }

    return sum;
  }

  def cdfBoost(x: Double, k: Double, l: Double, max_iter: Int = 5000, errtol: Double = 1e-8): Double = {
    var result = 0.0
    if (l == 0) {
      result = GammaFunction.cdfChiSquared(k, x); //central FIXME
    } else if (x > k + l) {
      // Complement is the smaller of the two:
      result = 1 - cdfQBentonKrishnamoorthy(x, k, l, max_iter, errtol)
    } else if (l < 200) {
      result = cdfDing(x, k, l, max_iter, errtol)
    } else {
      // For largers values of the non-centrality
      // parameter Ding's method will consume an
      // extra-ordinary number of terms, and worse
      // may return zero when the result is in fact
      // finite, use Krishnamoorthy's method instead:
      result = cdfBentonKrishnamoorthy(x, k, l, max_iter, errtol)
    }
    return result
  }

  def cdf(x: Double, df: Double, ncp: Double, max_iter: Int = 5000, errtol: Double = 1e-8): Double = {
//    return cdfOpenGamma(x, df, ncp)
    return cdfBoost(x, df, ncp, max_iter, errtol)
//    return cdfBentonKrishnamoorthy(x,df,ncp)
  }

}
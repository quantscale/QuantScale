package quantscale.analytic

import quantscale.fdm.Epsilon
import org.apache.commons.math3.special.Gamma

object GammaFunction {
  private val c1_ = 76.18009172947146;
  private val c2_ = -86.50532032941677;
  private val c3_ = 24.01409824083091;
  private val c4_ = -1.231739572450155;
  private val c5_ = 0.1208650973866179e-2;
  private val c6_ = -0.5395239384953e-5;

  def gammaPDerivative(a: Double, x: Double): Double = {
    if (a <= 0) throw new RuntimeException("argument a must be > 0 but was " + a)
    if (x < 0) throw new RuntimeException("argument x must be >= 0 but was " + x)
    if (x == 0) {
      if (a > 1) return 0 else {
        if (a == 1) return 1 else throw new RuntimeException("Overflow")
      }
    }
    var f1 = 0.0
    if (x <= Epsilon.MACHINE_EPSILON_SQRT) {
      // Oh dear, have to use logs, should be free of cancellation errors though:
      f1 = math.exp(a * math.log(x) - x - Gamma.logGamma(a));
    } else {
      // direct calculation, no danger of overflow as gamma(a) < 1/a
      // for small a.
      f1 = math.pow(x, a) * math.exp(-x) / Gamma.gamma(a);
    }
    return f1 / x
  }

  def logValue(x: Double): Double = {
    if (x <= 0.0) throw new RuntimeException("positive argument required");
    var temp = x + 5.5;
    temp -= (x + 0.5) * math.log(temp);
    var ser = 1.000000000190015;
    ser += c1_ / (x + 1.0);
    ser += c2_ / (x + 2.0);
    ser += c3_ / (x + 3.0);
    ser += c4_ / (x + 4.0);
    ser += c5_ / (x + 5.0);
    ser += c6_ / (x + 6.0);

    return -temp + math.log(2.5066282746310005 * ser / x);
  }

  def cdfChiSquared(degrees_of_freedom: Double, x: Double): Double = {
    return Gamma.regularizedGammaP(degrees_of_freedom / 2, x / 2); //is it same as gammap?
  }
}
package quantscale.vol

import quantscale.analytic.CumulativeNormalDistribution
import quantscale.math.AS241InvCND

trait BlackVolatilitySolver {

}

object Li2011ImpliedVolatilityGuess {
  private val m00 = -0.00006103098165; private val n00 = 1;
  private val m01 = 5.33967643357688; private val n01 = 22.96302109010794;
  private val m10 = -0.40661990365427; private val n10 = -0.48466536361620;
  private val m02 = 3.25023425332360; private val n02 = -0.77268824532468;
  private val m11 = -36.19405221599028; private val n11 = -1.34102279982050;
  private val m20 = 0.08975394404851; private val n20 = 0.43027619553168;
  private val m03 = 83.84593224417796; private val n03 = -5.70531500645109;
  private val m12 = 41.21772632732834; private val n12 = 2.45782574294244;
  private val m21 = 3.83815885394565; private val n21 = -0.04763802358853;
  private val m30 = -0.21619763215668; private val n30 = -0.03326944290044;

  //TODO restrict x to the well defined region where this approx is valid
  def impliedVolatilitySqrtTime(x: Double, c: Double): Double = {
    val num = m00 + m10 * x + m20 * x * x + m30 * x * x * x + c * (m01 + m11 * x + m21 * x * x) + c * c * (m02 + m12 * x) + c * c * c * m03
    val den = n00 + n10 * x + n20 * x * x + n30 * x * x * x + c * (n01 + n11 * x + n21 * x * x) + c * c * (n02 + n12 * x) + c * c * c * n03
    return num / den
  }
}

abstract class Li2011BlackVolatilitySolver(val tolerance: Double = 1e-6, val maxIterations: Int = 100) {
  var iterations = 0
  var c = 0.0
  var ex = 0.0
  var x = 0.0

  var reportException = false

  private[vol] def init(isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double) {
    c = price / forward * df
    ex = forward / strike
    if (!isCall) c = c + 1 - 1.0 / ex //put call parity
    val x0 = math.log(ex)
    x = x0
    if (x0 > 0) {
      //use c(-x0, v)
      c = ex * c + 1 - ex //in out duality
      x = -x0
      ex = 1.0 / ex
    }
  }
}

class Li2011SORBlackVolatilitySolver(tolerance: Double = 1e-6, val omega: Double = 1.0, maxIterations: Int = 100) extends Li2011BlackVolatilitySolver(tolerance, maxIterations) {

  private def compute(v0: Double, sqrttte: Double): Double = {
    var v = v0
    var Nm = CumulativeNormalDistribution.value(x / v - v / 2) / ex
    var Np = CumulativeNormalDistribution.value(x / v + v / 2)
    var cEstimate = Np - Nm
    iterations = 0
    while (math.abs(c - cEstimate) > tolerance && !v.isNaN()) {
      val F = c + Nm + omega * Np
      val Nm1 = AS241InvCND.value(F / (1 + omega))
      val G = Nm1 + math.sqrt(Nm1 * Nm1 - 2 * x)
      v = G //vk+1 = G(vk)
      Nm = CumulativeNormalDistribution.value(x / v - v / 2) / ex
      Np = CumulativeNormalDistribution.value(x / v + v / 2)
      cEstimate = Np - Nm
      iterations += 1
      if (iterations > maxIterations) {

        if (reportException) {
          throw new RuntimeException("maximum number of iterations exceeded")
        } else {
          v = Double.NaN
        }
      }
    }
    val sigma = v / sqrttte
    return sigma
  }

  def impliedVolatility(impliedVolGuess: Double, isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double, tte: Double): Double = {
    init(isCall, strike, price, forward, df)
    val sqrttte = math.sqrt(tte)
    val v0 = impliedVolGuess * sqrttte
    return compute(v0, sqrttte)
  }

  def impliedVolatility(isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double, tte: Double): Double = {
    init(isCall, strike, price, forward, df)
    val sqrttte = math.sqrt(tte)
    val v0 = Li2011ImpliedVolatilityGuess.impliedVolatilitySqrtTime(x, c)
    return compute(v0, sqrttte)
  }
}

class Li2011SORDRBlackVolatilitySolver(tolerance: Double = 1e-6, maxIterations: Int = 100) extends Li2011BlackVolatilitySolver(tolerance, maxIterations) {

  def impliedVolatility(impliedVolGuess: Double, isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double, tte: Double): Double = {
    init(isCall, strike, price, forward, df)
    val sqrttte = math.sqrt(tte)
    val v0 = impliedVolGuess * sqrttte
    return compute(v0, sqrttte)
  }

  def impliedVolatility(isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double, tte: Double): Double = {
    init(isCall, strike, price, forward, df)
    val sqrttte = math.sqrt(tte)
    val v0 = Li2011ImpliedVolatilityGuess.impliedVolatilitySqrtTime(x, c)
    return compute(v0, sqrttte)
  }

  private def compute(v0: Double, sqrttte: Double): Double = {
    var v = v0
    var Nm = CumulativeNormalDistribution.value(x / v - v / 2) / ex
    var Np = CumulativeNormalDistribution.value(x / v + v / 2)

    var cEstimate = Np - Nm
    var phi = (v * v + 2 * x) / (v * v - 2 * x)
    var omega = phi
    iterations = 0
    while (math.abs(c - cEstimate) > tolerance && !v.isNaN()) {
      val F = c + Nm + omega * Np
      val Nm1 = AS241InvCND.value(F / (1 + omega))
      val G = Nm1 + math.sqrt(Nm1 * Nm1 - 2 * x)
      v = G //vk+1 = G(vk)
      Nm = CumulativeNormalDistribution.value(x / v - v / 2) / ex
      Np = CumulativeNormalDistribution.value(x / v + v / 2)
      cEstimate = Np - Nm
      phi = (v * v + 2 * x) / (v * v - 2 * x)
      omega = phi
      iterations += 1
      if (iterations > maxIterations) {
        if (reportException) {
          throw new RuntimeException("maximum number of iterations exceeded")
        } else {
          v = Double.NaN
        }
      }
    }
    val sigma = v / sqrttte
    return sigma
  }
}

class Li2011SORTSBlackVolatilitySolver(tolerance: Double = 1e-6, maxIterations: Int = 100) extends Li2011BlackVolatilitySolver(tolerance, maxIterations) {
  def impliedVolatility(impliedVolGuess: Double, isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double, tte: Double): Double = {
    init(isCall, strike, price, forward, df)
    val sqrttte = math.sqrt(tte)
    val v0 = impliedVolGuess * sqrttte
    return compute(v0, sqrttte)
  }

  def impliedVolatility(isCall: Boolean, strike: Double, price: Double, forward: Double, df: Double, tte: Double): Double = {
    init(isCall, strike, price, forward, df)
    val sqrttte = math.sqrt(tte)
    val v0 = Li2011ImpliedVolatilityGuess.impliedVolatilitySqrtTime(x, c)
    return compute(v0, sqrttte)
  }

  private def compute(v0: Double, sqrttte: Double): Double = {
    var v = v0
    var Nm = CumulativeNormalDistribution.value(x / v - v / 2) / ex
    var Np = CumulativeNormalDistribution.value(x / v + v / 2)

    var cEstimate = Np - Nm

    var omega = 1.0
    iterations = 0
    while (math.abs(c - cEstimate) > tolerance && !v.isNaN) {
      var phi = (v * v + 2 * x) / (v * v - 2 * x)
      var alpha = (1 + omega) / (1 + phi)
      val F = c + Nm + omega * Np
      val Nm1 = AS241InvCND.value(F / (1 + omega))
      val G = Nm1 + math.sqrt(Nm1 * Nm1 - 2 * x)
      v = alpha * G + (1 - alpha) * v
      Nm = CumulativeNormalDistribution.value(x / v - v / 2) / ex
      Np = CumulativeNormalDistribution.value(x / v + v / 2)
      cEstimate = Np - Nm
      iterations += 1
      if (iterations > maxIterations) {
        if (reportException) {
          throw new RuntimeException("maximum number of iterations exceeded")
        } else {
          v = Double.NaN
        }
      }
    }
    val sigma = v / sqrttte
    return sigma
  }
}
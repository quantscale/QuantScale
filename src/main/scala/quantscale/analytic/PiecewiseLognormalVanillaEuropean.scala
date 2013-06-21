package quantscale.analytic

import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.solvers.SecantSolver
import org.apache.commons.math3.analysis.DifferentiableUnivariateFunction
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver
import org.apache.commons.math3.analysis.solvers.NewtonSolver

trait DiscountFactorProvider {
  def getDf(t: Double): Double
}

private class F1(strike: Double,
                 spot: Double,
                 variance: Double,
                 driftDf: Double,
                 dividends: Array[Dividend],
                 discount: DiscountFactorProvider,
                 maturity: Double) extends DifferentiableUnivariateFunction {
  private val sigma = math.sqrt(variance / maturity)
  private val discountDf = discount.getDf(maturity)

  def value(d1: Double): Double = {
    var sum = 0.0
    for (dividend <- dividends) {
      sum += dividend.value * discount.getDf(dividend.paymentTime) * math.exp(sigma * dividend.exTime / math.sqrt(maturity) * d1 - sigma * sigma * 0.5 * dividend.exTime * dividend.exTime / maturity)
    }
    return spot / driftDf * discountDf - sum - strike * discountDf * math.exp(math.sqrt(variance) * d1 - 0.5 * variance)
  }
  
  def derivative() : UnivariateFunction = {
    val f = this
    val eps = 1e-5
    return new UnivariateFunction() {
      def value(x: Double) : Double = {
        return (f.value(x+eps)-f.value(x-eps))/(2*eps)
      }
    }
  }
}
/**
 * Zhang "Fast Valuation for European Options with Cash dividends" published in Wilmott May 2011.
 */
object PiecewiseLognormalVanillaEuropean {
  //maturity, dividend times and discountCurve assume same daycount convention
  def priceEuropeanVanilla(
    isCall: Boolean,
    strike: Double,
    spot: Double,
    variance: Double,
    driftDf: Double,
    dividends: Array[Dividend],
    discountCurve: DiscountFactorProvider,
    maturity: Double): Double = {

    val sign = if (isCall) 1 else -1;
    val sqrtVar = math.sqrt(variance);
    val forward = spot / driftDf;
    val d1 = computeD1(strike, spot, variance, driftDf, dividends, discountCurve, maturity)
    val d2 = d1 - sqrtVar;
    val nd1 = CumulativeNormalDistribution.value(sign * d1);
    val nd2 = CumulativeNormalDistribution.value(sign * d2);
    var sum = 0.0
    for (dividend <- dividends) {
      val d3 = d1 - sqrtVar * dividend.exTime / maturity
      sum += dividend.value * discountCurve.getDf(dividend.paymentTime) * CumulativeNormalDistribution.value(d3)
    }
    val discountDf = discountCurve.getDf(maturity)
    val price = sign * (discountDf * (forward * nd1 - strike * nd2) - sum);
    return price;
  }

  private def computeD1(strike: Double,
                        spot: Double,
                        variance: Double,
                        driftDf: Double,
                        dividends: Array[Dividend],
                        discountCurve: DiscountFactorProvider,
                        maturity: Double): Double = {
    val solver = new NewtonSolver()
    val f = new F1(strike, spot, variance, driftDf, dividends, discountCurve, maturity)
    val forward = spot / driftDf
    var minValue = Double.MaxValue
    for (dividend <- dividends) {
      minValue = math.min(minValue, (math.log(spot / (dividend.value * discountCurve.getDf(dividend.exTime))) + 0.5 * variance * dividend.exTime / maturity * dividend.exTime / maturity) / math.sqrt(variance * dividend.exTime / maturity))
    }
    val startValue = math.min((math.log(forward / strike) + 0.5 * variance) / math.sqrt(variance), minValue)
    return solver.solve(10, f, startValue)
  }
  
}
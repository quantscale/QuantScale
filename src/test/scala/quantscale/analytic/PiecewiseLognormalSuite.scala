package test.quantscale.analytic

import org.scalatest.FunSuite
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import quantscale.analytic.Dividend
import quantscale.analytic.DiscountFactorProvider
import quantscale.analytic.PiecewiseLognormalVanillaEuropean
import quantscale.analytic.BlackScholesVanillaEuropean
@RunWith(classOf[JUnitRunner])
class PiecewiseLognormalSuite extends FunSuite {
  test("vanilla-against-zhang") {
    val spot = 100.0
    val strike = 100.0
    val isCall = true
    val vol = 0.30
    val maturity = 0.25
    val dividends = Array(new Dividend(1.0, 0.1, 0.1))
    val variance = vol * vol * maturity
    val r = 0.05
    val driftDf = math.exp(-r * maturity)
    val discountCurve = new DiscountFactorProvider() {
      def getDf(t: Double) = math.exp(-r * t)
    }
    val price = PiecewiseLognormalVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance, driftDf, dividends, discountCurve, maturity)

    val q = dividends(0).value / spot
    val stdPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot - dividends(0).value * math.exp(-r * dividends(0).paymentTime), variance, driftDf, driftDf)
    println(stdPrice)
    println(price)
    assert(math.abs(stdPrice - 6.03608) < 1e-5, "equal forward price does not match, was " + stdPrice)
    assert(math.abs(price - 6.05984) < 1e-5, "first order price does not match, was " + price)

  }
}
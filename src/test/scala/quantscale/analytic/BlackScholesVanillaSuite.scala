package quantscale.analytic

import org.scalatest.{FunSuite}

class BlackScholesVanillaSuite extends FunSuite {

  test("adjoint-price") {
    val isCall = true
    val strike = 99.0
    val spot = 100.0
    val tte = 0.9
    val variance = new Variance(20.0, tte)
    val drift = new DiscountFactor(0.01, tte)
    val discount = new DiscountFactor(0.02, tte)
    val refPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance.value, drift.value, discount.value)

    val measures: Array[BSMMeasure] = Array(new Price)
    BlackScholesVanillaEuropean.priceEuropeanVanillaAdjoint(isCall, strike, spot, variance, drift, discount, measures)
    println(refPrice, " ", measures(0).value)
    assert(math.abs(refPrice - measures(0).value) < 1e-8)
  }

  test("adjoint-greeks") {
    val isCall = true
    val strike = 99.0

    val tte = 0.9
    val variance = new Variance(0.20, tte)
    val drift = new DiscountFactor(0.01, tte)
    val discount = new DiscountFactor(0.02, tte)

    val price = new Price()
    val delta = new DeltaGreek()
    val vega = new VegaGreek
    val rho = new RhoGreek
    val rho2 = new Rho2Greek
    val theta = new ThetaGreek

    var i = 1000000
    var refElapsed = 0.0
    var elapsed = 0.0
    while (i > 0) {
      val spot = 100.0 + i*50/1000000.0
      val eps = 1e-7
      val refStart = System.nanoTime()
      val refPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance.value, drift.value, discount.value)
      val refDelta = (BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot + eps, variance.value, drift.value, discount.value) - refPrice) / eps
      val refVega = (BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, new Variance(variance.vol + eps, variance.time).value, drift.value, discount.value) - refPrice) / eps
      val refRho = (BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance.value, drift.value, new DiscountFactor(discount.rate + eps, discount.time).value) - refPrice) / eps
      val refRho2 = (BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance.value, new DiscountFactor(drift.rate + eps, drift.time).value, discount.value) - refPrice) / eps
      val refTheta = (BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, new Variance(variance.vol, variance.time + eps).value, new DiscountFactor(drift.rate, drift.time + eps).value, new DiscountFactor(discount.rate, discount.time + eps).value) - refPrice) / eps
      if (i>10) refElapsed += System.nanoTime() - refStart
      val start = System.nanoTime()
      val measures: Array[BSMMeasure] = Array(price, delta, vega, rho, rho2, theta)
      BlackScholesVanillaEuropean.priceEuropeanVanillaAdjoint(isCall, strike, spot, variance, drift, discount, measures)
      //      println(refPrice, " ", price.value)
      //      println(refDelta, " ", delta.value)
      //      println(refVega, " ", vega.value)
      //      println(refRho, " ",rho.value)
      //      println(refRho2, " ", rho2.value)
      //      println(refTheta, " ", theta.value)
      if (i >10) elapsed += System.nanoTime() - start
      assert(math.abs(refPrice - price.value) < 1e-6, "price")
      //expectResult(refDelta, delta.value, 1e-6, "delta")
      assert(math.abs(refDelta - delta.value) < 1e-5, "delta "+refDelta + " "+ delta.value)
      assert(math.abs(refVega - vega.value) < 1e-5, "vega")
      assert(math.abs(refRho - rho.value) < 1e-5, "rho")
      assert(math.abs(refRho2 - rho2.value) < 1e-4, "rho2")
      assert(math.abs(refTheta - theta.value) < 1e-5, "theta")
      i -= 1
    }
    println("FD time=" + refElapsed / 1e9 + " Adjoint time=" + elapsed / 1e9)

  }

}

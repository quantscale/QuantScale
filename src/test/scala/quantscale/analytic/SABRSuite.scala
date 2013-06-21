package test.quantscale.analytic

import org.junit.runner.RunWith
import org.scalatest.FunSuite
import quantscale.analytic.BachelierVanillaEuropean
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.analytic.SABRModelSpec
import quantscale.analytic.SABRVanilla
import quantscale.fdm.mesh.Mesh1DBoundaries
import quantscale.fdm.mesh.UniformMesh1D
import quantscale.fdm.sabr._
import quantscale.math.CubicPP
import quantscale.math.CubicSpline
import quantscale.vol.Li2011SORBlackVolatilitySolver
import quantscale.vol.Li2011SORDRBlackVolatilitySolver
import quantscale.vol.Li2011SORTSBlackVolatilitySolver
import org.scalatest.junit.JUnitRunner
import quantscale.fdm.Epsilon

@RunWith(classOf[JUnitRunner])
class SABRSuite extends FunSuite {

  test("li-sordr-iv") {
    val forward = 100.0
    val strikes = Array(0.2 * forward, 0.5 * forward, 0.9 * forward, 1.0 * forward, 1.1 * forward, 1.5 * forward)
    val sigma = 0.4
    val isCall = true
    val tte = 0.5
    val df = 1.0
    for (strike <- strikes) {
      val price = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, sigma * sigma * tte, 1.0, df)
      val solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
      val vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte);
      println(strike + " " + vol + " in " + solver.iterations)
      assert(math.abs(vol - sigma) < 1e-4)
    }
  }

  test("li-sorts-iv") {
    val forward = 100.0
    val strikes = Array(0.2 * forward, 0.5 * forward, 0.9 * forward, 1.0 * forward, 1.1 * forward, 1.5 * forward)
    val sigma = 0.4
    val isCall = true
    val tte = 0.5
    val df = 1.0
    for (strike <- strikes) {
      val price = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, sigma * sigma * tte, 1.0, df)
      val solver = new Li2011SORTSBlackVolatilitySolver(1e-12)
      val vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte);
      println(strike + " " + vol + " in " + solver.iterations)
      assert(math.abs(vol - sigma) < 1e-4)
    }
  }
  test("li-sor-iv") {
    val forward = 1.0
    val strikes = Array(0.5, 0.9, 1.0, 1.1, 1.5)
    val sigma = 0.4
    val isCall = true
    val tte = 1.0
    val df = 1.0
    for (strike <- strikes) {
      val price = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, sigma * sigma * tte, 1.0, df)
      val solver = new Li2011SORBlackVolatilitySolver(1e-6, 0.5)
      val vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte);
      println(strike + " " + vol + " in " + solver.iterations)
      assert(math.abs(vol - sigma) < 1e-4)
    }
  }

  test("li-hagan") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 2.0;
    val strike = forward / 10;
    val df = 1.0

    val isCall = false
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val blackVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    val price = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, blackVol * blackVol * tte, 1.0, df)
    println(blackVol + " " + price)
    val solver = new Li2011SORTSBlackVolatilitySolver(1e-12)
    val solvedVol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
    println(solvedVol)
    assert(math.abs(solvedVol - blackVol) < 1e-4)
  }

  test("hagan") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 2.0;
    val strike = 0.01;
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val blackVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    println(blackVol)
    assert(math.abs(0.6444840588326359 - blackVol) < 1e-15)
  }

  test("hagan-fd-price") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 2.0;
    val strike = 0.01;
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val blackVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    val isCall = true
    val blackPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, blackVol * blackVol * tte, 1.0, 1.0)
    println(blackPrice)

    val pdePrice = new HaganSABRDensitySolver(spec, forward, tte).price(isCall, strike)
    println(pdePrice)
  }

  test("hagan-benhamou") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 2.0;
    val strike = 0.03;
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val blackVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    println(blackVol)
    val isCall = true
    val blackPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, blackVol * blackVol * tte, 1.0, 1.0)
    val benPrice = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
    println(blackPrice + " " + benPrice + " " + math.abs(blackPrice - benPrice))
    val normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
    val normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
    println(blackPrice + " " + normalPrice + " " + (blackPrice - normalPrice))
    val ahPrice = SABRVanilla.priceAndreasenHuge(isCall, strike, spec, forward, tte)
    println(blackPrice + " " + ahPrice + " " + (blackPrice - ahPrice))

  }

  test("hagan-normal") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 2.0;
    val strike = 0.02;
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val blackVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    val isCall = false
    val blackPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, blackVol * blackVol * tte, 1.0, 1.0)
    val normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
    val normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
    println(blackPrice + " " + normalPrice + " " + (blackPrice - normalPrice))
    assert(math.abs(normalPrice - blackPrice) < 2e-6)
    val solver = new Li2011SORTSBlackVolatilitySolver(1e-12)
    val volB = solver.impliedVolatility(isCall, strike, blackPrice, forward, 1.0, tte)
    val volN = solver.impliedVolatility(isCall, strike, normalPrice, forward, 1.0, tte)
    println(volB + " " + volN + " " + math.abs(volB - volN))

  }

  test("andreasenhuge") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 0.1
    val strike = 0.01;
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)

    val ahVol = SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, strike, tte)
    println(haganVol + " " + ahVol + " " + (haganVol - ahVol))
    assert(math.abs(haganVol - ahVol) < 3e-4)
  }

  test("andreasenhuge-price") {
    val alpha = 0.0758194;
    val nu = 0.1;
    val beta = 0.5;
    val rho = -0.1;
    val forward = 0.02;
    val tte = 2.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val splineCall = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte)
    val splinePut = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte)
    println("IsCall Strike HaganPrice AndreasenHugePrice PriceError HaganVol NormalHaganVol AndreasenHugeVol VolError")
    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 50.0)
      val isCall = strike >= forward
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      val haganPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, haganVol * haganVol * tte, 1.0, df)
      val normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      val normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)

      val ahPrice = if (isCall) splineCall.value(strike) else splinePut.value(strike)
      print(isCall + " " + strike + " " + haganPrice + " " + ahPrice + " " + (haganPrice - ahPrice))
      var solver = new Li2011SORDRBlackVolatilitySolver(1e-8)
      val haganVol2 = solver.impliedVolatility(isCall, strike, haganPrice, forward, df, tte)
      assert(math.abs(haganVol - haganVol2) < 1e-5)
      val normalVol2 = solver.impliedVolatility(isCall, strike, normalPrice, forward, df, tte)
      val ahVol = solver.impliedVolatility(isCall, strike, ahPrice, forward, df, tte)
      println(" " + haganVol2 + " " + normalVol2 + " " + ahVol + " " + (haganVol - ahVol))
    }
  }

  test("andreasenhuge-density") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val splineCall = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte)
    val splinePut = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte)
    val splineCallRaw = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte, false)
    val splinePutRaw = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte, false)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    //vols
    val pde = new HaganTruncatedSABRDensitySolver(spec, forward, tte, 50, 10)
    pde.useRannacher = true
    println("vols")
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeApprox AndreasenHugeRaw AndreasenHuge HaganPDE")
    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 10.0)
      val isCall = strike >= forward
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      val normalVolScaled = normalVol * math.log(forward / strike) / (forward - strike)
      var normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      normalVol = solver.impliedVolatility(isCall, strike, normalPrice, forward, df, tte)
      val ahVolApprox = SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, strike, tte)
      val benPrice = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
      var benVol = Double.NaN
      try {
        benVol = solver.impliedVolatility(isCall, strike, benPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      val pdePrice = pde.price(isCall, strike)
      var pdeVol = Double.NaN
      try {
        pdeVol = solver.impliedVolatility(isCall, strike, pdePrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      val spline = if (isCall) splineCall else splinePut
      val ahPrice = spline.value(strike)
      var ahVol = Double.NaN
      try {
        ahVol = solver.impliedVolatility(isCall, strike, ahPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }

      val splineRaw = if (isCall) splineCallRaw else splinePutRaw
      val ahPriceRaw = splineRaw.value(strike)
      var ahVolRaw = Double.NaN
      try {
        ahVolRaw = solver.impliedVolatility(isCall, strike, ahPriceRaw, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      println(f"$strike%2.4f $haganVol%2.3f $normalVol%2.3f $normalVolScaled%2.3f $ahVolApprox%2.3f $ahVolRaw%2.3f $ahVol%2.3f $pdeVol%2.3f")
    }
    println("densities")
    //density
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeRaw AndreasenHuge Benhamou")

    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 20.0)
      val eps = 1e-3 * strike
      val isCall = strike >= forward
      val haganVolUp = SABRVanilla.impliedVolatilityHagan(spec, forward, strike + eps, tte)
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      val haganVolDown = SABRVanilla.impliedVolatilityHagan(spec, forward, strike - eps, tte)
      var variance = haganVolUp * haganVolUp * tte
      val callUp = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike + eps, forward, variance, 1.0, df)
      variance = haganVol * haganVol * tte
      val call = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, variance, 1.0, df)
      variance = haganVolDown * haganVolDown * tte
      val callDown = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike - eps, forward, variance, 1.0, df)
      val densityHagan = (callUp - 2 * call + callDown) / (eps * eps)

      var normalVolUp = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike + eps, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      var normalVolDown = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike - eps, tte)
      val callUpNormal = BachelierVanillaEuropean.price(isCall, strike + eps, forward, normalVolUp, tte)
      val callNormal = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      val callDownNormal = BachelierVanillaEuropean.price(isCall, strike - eps, forward, normalVolDown, tte)
      val densityNormal = (callUpNormal - 2 * callNormal + callDownNormal) / (eps * eps)

      normalVolUp *= math.log(forward / (strike + eps)) / (forward - (strike + eps))
      normalVol *= math.log(forward / strike) / (forward - strike)
      normalVolDown *= math.log(forward / (strike - eps)) / (forward - (strike - eps))
      variance = normalVolUp * normalVolUp * tte
      val callUpNormalS = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike + eps, forward, variance, 1.0, df)
      variance = normalVol * normalVol * tte
      val callNormalS = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, variance, 1.0, df)
      variance = normalVolDown * normalVolDown * tte
      val callDownNormalS = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike - eps, forward, variance, 1.0, df)
      val densityNormalS = (callUpNormalS - 2 * callNormalS + callDownNormalS) / (eps * eps)

      val spline = if (isCall) splineCall else splinePut
      val densityAh = spline.secondDerivative(strike)

      val splineRaw = if (isCall) splineCallRaw else splinePutRaw
      val densityAhRaw = splineRaw.secondDerivative(strike)

      val callUpBen = SABRVanilla.priceBenhamou(spec, isCall, strike + eps, forward, tte)
      val callBen = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
      val callDownBen = SABRVanilla.priceBenhamou(spec, isCall, strike - eps, forward, tte)
      val densityBen = (callUpBen - 2 * callBen + callDownBen) / (eps * eps)

      val callUpPDE = pde.price(isCall, strike + eps)
      val callPDE = pde.price(isCall, strike)
      val callDownPDE = pde.price(isCall, strike - eps)
      val densityPDE = (callUpPDE - 2 * callPDE + callDownPDE) / (eps * eps)
      //       println(f"$strike%2.4f $call%2.2e $callNormal%2.2e $callBen%2.2e")
      println(f"$strike%2.4f $densityHagan%2.2e $densityNormal%2.2e $densityNormalS%2.2e $densityAhRaw%2.2e $densityAh%2.2e $densityPDE%2.2e")
    }
  }

  test("hagan-density") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val splineCall = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte)
    val splinePut = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte)
    val splineCallRaw = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte, false)
    val splinePutRaw = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte, false)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    //vols
    //    for (i <- 0 to 10) {
    //      val start = System.nanoTime()
    //    val pde = new HaganSABRDensitySolver(spec, forward, tte,360,40)
    //    pde.solve()
    //    val mid = System.nanoTime()
    //    val pde2 = new HaganTRBDF2SABRDensitySolver(spec, forward, tte,360,40)
    //    pde2.solve()
    //    val end = System.nanoTime()
    //    println((mid-start)*1e-9+" "+(end-mid)*1e-9)
    //    }
    //    System.exit(1)
    val pde = new HaganSABRDensitySolver(spec, forward, tte, 500, 5, 5.0)
    pde.useSmoothing = false
    pde.useRannacher = false
    pde.solve()
    println(pde.computeCourantNumber + " " + pde.h)
    val pde2 = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, 500, 5, 5.0)

    println("vols")
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeApprox AndreasenHugeRaw AndreasenHuge HaganPDE HaganEulerPDE")
    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 10.0)
      val isCall = strike >= forward
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      val normalVolScaled = normalVol * math.log(forward / strike) / (forward - strike)
      var normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      normalVol = solver.impliedVolatility(isCall, strike, normalPrice, forward, df, tte)
      val ahVolApprox = SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, strike, tte)
      val benPrice = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
      var benVol = Double.NaN
      try {
        benVol = solver.impliedVolatility(isCall, strike, benPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      val pdePrice = pde.price(isCall, strike)
      var pdeVol = Double.NaN
      try {
        pdeVol = solver.impliedVolatility(isCall, strike, pdePrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }

      val pde2Price = pde2.price(isCall, strike)
      var pde2Vol = Double.NaN
      try {
        pde2Vol = solver.impliedVolatility(isCall, strike, pde2Price, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      val spline = if (isCall) splineCall else splinePut
      val ahPrice = spline.value(strike)
      var ahVol = Double.NaN
      try {
        ahVol = solver.impliedVolatility(isCall, strike, ahPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }

      val splineRaw = if (isCall) splineCallRaw else splinePutRaw
      val ahPriceRaw = splineRaw.value(strike)
      var ahVolRaw = Double.NaN
      try {
        ahVolRaw = solver.impliedVolatility(isCall, strike, ahPriceRaw, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      println(f"$strike%2.4f $haganVol%2.3f $normalVol%2.3f $normalVolScaled%2.3f $ahVolApprox%2.3f $ahVolRaw%2.3f $ahVol%2.3f $pdeVol%2.3f $pde2Vol%2.3f ")
    }
    println("densities")
    //density
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeRaw AndreasenHuge HaganPDE HaganEulerPDE")

    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 20.0)
      val eps = 1e-3 * strike
      val isCall = strike >= forward
      val haganVolUp = SABRVanilla.impliedVolatilityHagan(spec, forward, strike + eps, tte)
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      val haganVolDown = SABRVanilla.impliedVolatilityHagan(spec, forward, strike - eps, tte)
      var variance = haganVolUp * haganVolUp * tte
      val callUp = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike + eps, forward, variance, 1.0, df)
      variance = haganVol * haganVol * tte
      val call = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, variance, 1.0, df)
      variance = haganVolDown * haganVolDown * tte
      val callDown = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike - eps, forward, variance, 1.0, df)
      val densityHagan = (callUp - 2 * call + callDown) / (eps * eps)

      var normalVolUp = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike + eps, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      var normalVolDown = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike - eps, tte)
      val callUpNormal = BachelierVanillaEuropean.price(isCall, strike + eps, forward, normalVolUp, tte)
      val callNormal = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      val callDownNormal = BachelierVanillaEuropean.price(isCall, strike - eps, forward, normalVolDown, tte)
      val densityNormal = (callUpNormal - 2 * callNormal + callDownNormal) / (eps * eps)

      normalVolUp *= math.log(forward / (strike + eps)) / (forward - (strike + eps))
      normalVol *= math.log(forward / strike) / (forward - strike)
      normalVolDown *= math.log(forward / (strike - eps)) / (forward - (strike - eps))
      variance = normalVolUp * normalVolUp * tte
      val callUpNormalS = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike + eps, forward, variance, 1.0, df)
      variance = normalVol * normalVol * tte
      val callNormalS = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, variance, 1.0, df)
      variance = normalVolDown * normalVolDown * tte
      val callDownNormalS = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike - eps, forward, variance, 1.0, df)
      val densityNormalS = (callUpNormalS - 2 * callNormalS + callDownNormalS) / (eps * eps)

      val spline = if (isCall) splineCall else splinePut
      val densityAh = spline.secondDerivative(strike)

      val splineRaw = if (isCall) splineCallRaw else splinePutRaw
      val densityAhRaw = splineRaw.secondDerivative(strike)

      val callUpBen = SABRVanilla.priceBenhamou(spec, isCall, strike + eps, forward, tte)
      val callBen = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
      val callDownBen = SABRVanilla.priceBenhamou(spec, isCall, strike - eps, forward, tte)
      val densityBen = (callUpBen - 2 * callBen + callDownBen) / (eps * eps)

      val callUpPDE = pde.price(isCall, strike + eps)
      val callPDE = pde.price(isCall, strike)
      val callDownPDE = pde.price(isCall, strike - eps)
      val densityPDE = (callUpPDE - 2 * callPDE + callDownPDE) / (eps * eps)

      val callUpPDE2 = pde2.price(isCall, strike + eps)
      val callPDE2 = pde2.price(isCall, strike)
      val callDownPDE2 = pde2.price(isCall, strike - eps)
      val densityPDE2 = (callUpPDE2 - 2 * callPDE2 + callDownPDE2) / (eps * eps)
      //       println(f"$strike%2.4f $call%2.2e $callNormal%2.2e $callBen%2.2e")
      println(f"$strike%2.4f $densityHagan%2.2e $densityNormal%2.2e $densityNormalS%2.2e $densityAhRaw%2.2e $densityAh%2.2e $densityPDE%2.2e $densityPDE2%2.2e")
    }
  }

  test("alan") {
    val forward = 1.0
    val K = 0.8
    val tte = 10.0
    val nu = 1.00
    val beta = 0.30
    val alpha = 1.00
    val rho = 0.90
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val splineCall = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte)
    val splinePut = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte)
    val splineCallRaw = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte, false)
    val splinePutRaw = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte, false)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val pde =
    //      new HaganTruncatedSABRDensitySolver(spec, forward, tte, 10000, 1000)
      new HaganLMG3SABRDensitySolver(spec, forward, tte, 10000, 10, 1000.0)

    val volHagan = 100 * SABRVanilla.impliedVolatilityHagan(spec, forward, K, tte)
    val volAHShort = 100 * SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, K, tte)
    val priceAH = splinePut.value(K)
    val priceAHRaw = splinePutRaw.value(K)
    val pricePDE = pde.price(false, K)
    val normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, K, tte)
    val normalPrice = BachelierVanillaEuropean.price(true, K, forward, normalVol, tte)
    //    val volNormal = 100*solver.impliedVolatility(true, K, normalPrice, forward, df, tte)
    val volAH = 100 * solver.impliedVolatility(false, K, priceAH, forward, df, tte)
    val volAHRaw = 100 * solver.impliedVolatility(false, K, priceAHRaw, forward, df, tte)
    val volPDE = 100 * solver.impliedVolatility(false, K, pricePDE, forward, df, tte)

    val priceBen = SABRVanilla.priceBenhamou(spec, true, K, forward, tte) + K - forward
    println(BlackScholesVanillaEuropean.priceEuropeanVanilla(true, K, forward, volHagan * volHagan * tte, 1.0, df))
    println(pricePDE)
    println(priceBen)
    println(normalPrice)
    val volBen = Double.NaN //100* solver.impliedVolatility(false, K, priceBen, forward, df, tte)

    println(f"$volHagan%2.3f $volBen%2.3f $volAHShort%2.3f $volAH%2.3f $volPDE%2.3f")
  }

  def computeMaxError(solver: HaganSABRDensitySolver, strikes: Array[Double], forward: Double, tte: Double, refSolver: HaganSABRDensitySolver): Double = {
    var ivsolver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    var maxError = 0.0
    for (strike <- strikes) {
      val isCall = strike >= forward
      val price = solver.price(isCall, strike)
      val vol = ivsolver.impliedVolatility(isCall, strike, price, forward, 1.0, tte)
      val refprice = refSolver.price(isCall, strike)
      val refVol = ivsolver.impliedVolatility(isCall, strike, refprice, forward, 1.0, tte)
      val error = math.abs(vol - refVol)
      maxError = math.max(maxError, error)
    }
    return maxError
  }

  def computeMaxError(solver: HaganSABRDensitySolver, Qspline: CubicPP): (Double, Double) = {
    var maxError = 0.0
    var Fmax = 0.0
    val F = solver.F
    val Q0 = solver.Q0

    for (i <- 1 until F.size - 2) {
      val Qref = Qspline.value(F(i))
      var error = Q0(i) - Qref
      if (error < 0) error = -error
      if (error > maxError) {
        maxError = error
        Fmax = F(i)
      }
    }
    return (maxError, Fmax)
  }

  test("pde-density-accuracy-ah") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0

    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0 * forward

    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280*4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    println("courant=" + refPDE.computeCourantNumber)
    val yPrime = Array.ofDim[Double](refPDE.F.size)
    //    CubicSpline.computeHarmonicFirstDerivativePCHIM(refPDE.F, refPDE.Q0, yPrime)
    CubicSpline.computeC2FirstDerivative(refPDE.F, refPDE.Q0, yPrime)
    val Qspline = CubicSpline.makeHermiteSpline(refPDE.F, refPDE.Q0, yPrime)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
    for (l <- 0 to 5) {
      println("SpaceSteps TimeSteps Scheme MaxError FError Time")

      for (spaceStep <- spaceSteps) {
        var buffer = ""
        for (timeStep <- timeSteps) {
          var startTime = System.nanoTime()
          val cn = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          cn.useSmoothing = false
          cn.useRannacher = false
          cn.solve()
          var endTime = System.nanoTime()
          var elapsed = (endTime - startTime) * 1e-9
          var err = computeMaxError(cn, Qspline)
          var error = err._1
          var f = err._2
          buffer += f"$spaceStep $timeStep CN $error%2.1e $f%2.4f $elapsed%2.1e\n"
          startTime = System.nanoTime()
          var ran = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          ran.useRannacher = true
          ran.solve()
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          err = computeMaxError(ran, Qspline)
          error = err._1
          f = err._2
          buffer += f"$spaceStep $timeStep RAN $error%2.1e $f%2.4f $elapsed%2.1e\n"
          startTime = System.nanoTime()
          val lmg2 = new HaganLMG2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          lmg2.solve()
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          err = computeMaxError(lmg2, Qspline)
          error = err._1
          f = err._2
          buffer += f"$spaceStep $timeStep LMG2 $error%2.1e $f%2.4f $elapsed%2.1e\n"
          startTime = System.nanoTime()
          val lmg3 = new HaganLMG3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          lmg3.solve()
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          err = computeMaxError(lmg3, Qspline)
          error = err._1
          f = err._2
          buffer += f"$spaceStep $timeStep LMG3 $error%2.1e $f%2.4f $elapsed%2.1e\n"
          startTime = System.nanoTime()
          val trbdf2 = new HaganTRBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          trbdf2.solve()
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          err = computeMaxError(trbdf2, Qspline)
          error = err._1
          f = err._2
          buffer += f"$spaceStep $timeStep TRBDF2 $error%2.1e $f%2.4f $elapsed%2.1e\n"
          startTime = System.nanoTime()
          val cnbdf2 = new HaganCNBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          cnbdf2.solve()
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          err = computeMaxError(cnbdf2, Qspline)
          error = err._1
          f = err._2
          buffer += f"$spaceStep $timeStep CNBDF2 $error%2.1e $f%2.4f $elapsed%2.1e\n"
          startTime = System.nanoTime()
          val re = new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          re.solve()
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          err = computeMaxError(re, Qspline)
          error = err._1
          f = err._2
          buffer += f"$spaceStep $timeStep RE $error%2.1e $f%2.4f $elapsed%2.1e\n"
        }
        print(buffer)
      }
    }
  }

  test("pde-vol-accuracy-ah") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0 * forward
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.05, forward * 3)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 500, 1280 * 4, FmaxTruncation)
    //  val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)

    val schemes = Array("CN","RAN","LMG2","LMG3","LS","TRBDF2","TRBDF3")
    for (l <- 0 to 7) {
      println("SpaceSteps TimeSteps Scheme MaxError Time")
      for (scheme <- schemes) {
      for (spaceStep <- spaceSteps) {
        var buffer = ""
        for (timeStep <- timeSteps) {
          var startTime = System.nanoTime()
          val cn = createHaganSABRSolver(scheme,spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          cn.solve()
          var endTime = System.nanoTime()
          var elapsed = (endTime - startTime) * 1e-9
          var error = computeMaxError(cn, strikes, forward, tte, refPDE)
          buffer += f"$spaceStep $timeStep $scheme $error%2.1e $elapsed%2.1e\n"
        }
        print(buffer)
      }
      }
    }
  }

  private def createErrorString(isLatex: Boolean, scheme: String, spaceStep: Int, timeStep: Int, error: Double, errorLocation: Double, elapsed: Double): String = {
    return if (isLatex) f"$spaceStep & $timeStep & $scheme & $error%2.1e & $elapsed%2.1e\\\\\n"
    else f"$spaceStep $timeStep $scheme $error%2.1e $errorLocation%2.4f $elapsed%2.1e\n"
  }

  test("pde-density-accuracy-hagan") {
    val isLatex = false
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 500, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    val Qspline = CubicSpline.makeCubicSpline(refPDE.F, refPDE.Q0)

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "CNBDF3", "RE")
    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
    for (l <- 0 to 9) {
      //            println("SpaceSteps TimeSteps Scheme MaxError FError Time")
      println("SpaceSteps & TimeSteps & Scheme & MaxError & Time\\\\")
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {
          var buffer = ""
          for (timeStep <- timeSteps) {
            var startTime = System.nanoTime()
            val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
            cn.solve()
            var endTime = System.nanoTime()
            var elapsed = (endTime - startTime) * 1e-9
            var err = computeMaxError(cn, Qspline)
            var error = err._1
            var f = err._2
            buffer += createErrorString(isLatex, scheme, spaceStep, timeStep, error, f, elapsed)
          }
          print(buffer)
        }
      }
    }
  }

  test("pde-vol-accuracy-hagan") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganSABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 64, FmaxTruncation)

    refPDE.solve()
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "RE")
    for (l <- 0 until 5) {
      println("SpaceSteps TimeSteps Scheme MaxError Time")
      for (spaceStep <- spaceSteps) {
        var buffer = ""
        for (scheme <- schemes) {
          for (timeStep <- timeSteps) {
            var startTime = System.nanoTime()
            val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
            cn.solve()
            var endTime = System.nanoTime()
            var elapsed = (endTime - startTime) * 1e-9
            var error = computeMaxError(cn, strikes, forward, tte, refPDE)
            buffer += f"$spaceStep $timeStep $scheme $error%2.1e $elapsed%2.1e\n"
          }
          print(buffer)
        }
      }
    }
  }

  test("pde-vs-formula-performance-hagan") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val strike = forward

    val FmaxTruncation = 5
    //warmup
    var vol = 0.0
    for (i <- 0 until 100) {
      val pde =
      //      new HaganSABRDensitySolver(spec, forward, tte, 125, 5, FmaxTruncation)
      //    pde.useRannacher=true
        new HaganTRBDF2SABRDensitySolver(spec, forward, tte, 62, 5, FmaxTruncation)
      val price = pde.price(true, strike)
      vol = solver.impliedVolatility(true, strike, price, forward, df, tte)
    }
    println(vol)
    var startTime = System.nanoTime()
    for (i <- 0 until 1000) {
      val pde =
      //           new HaganSABRDensitySolver(spec, forward, tte, 125, 5, FmaxTruncation)
      //    pde.useRannacher=true

        new HaganTRBDF2SABRDensitySolver(spec, forward, tte, 62, 5, FmaxTruncation)
      val price = pde.price(true, strike)
      vol = solver.impliedVolatility(true, strike, price, forward, df, tte)
    }
    var endTime = System.nanoTime()
    val pdeElapsed = (endTime - startTime) * 1e-9
    var formulaVol = 0.0
    for (i <- 0 until 100) {
      formulaVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    }
    println(formulaVol)
    startTime = System.nanoTime()
    for (i <- 0 until 1000) {
      formulaVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    }
    endTime = System.nanoTime()
    val formulaElapsed = (endTime - startTime) * 1e-9

    val ratio = pdeElapsed / formulaElapsed
    println(f"$pdeElapsed%2.4f $formulaElapsed%2.4f $ratio%2.1f")
  }

  test("pde-vs-formula-performance-ah") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val strike = forward

    val FmaxTruncation = 10 * forward
    //warmup
    var vol = 0.0
    for (i <- 0 until 100) {
      val pde =
      //      new HaganSABRDensitySolver(spec, forward, tte, 125, 5, FmaxTruncation)
      //    pde.useRannacher=true
        new HaganLMG3SABRDensitySolver(spec, forward, tte, 250, 2, FmaxTruncation)
      val price = pde.price(true, strike)
      vol = solver.impliedVolatility(true, strike, price, forward, df, tte)
    }
    println(vol)
    var startTime = System.nanoTime()
    for (i <- 0 until 1000) {
      val pde =
      //           new HaganSABRDensitySolver(spec, forward, tte, 125, 5, FmaxTruncation)
      //    pde.useRannacher=true

        new HaganLMG3SABRDensitySolver(spec, forward, tte, 250, 2, FmaxTruncation)
      val price = pde.price(true, strike)
      vol = solver.impliedVolatility(true, strike, price, forward, df, tte)
    }
    var endTime = System.nanoTime()
    val pdeElapsed = (endTime - startTime) * 1e-9
    var formulaVol = 0.0
    for (i <- 0 until 100) {
      formulaVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    }
    println(formulaVol)
    startTime = System.nanoTime()
    for (i <- 0 until 1000) {
      formulaVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
    }
    endTime = System.nanoTime()
    val formulaElapsed = (endTime - startTime) * 1e-9

    val ratio = pdeElapsed / formulaElapsed
    println(f"$pdeElapsed%2.4f $formulaElapsed%2.4f $ratio%2.1f")
  }

  test("pde-performance-hagan") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    var isCall = false
    //    val strikeLow = forward*0.1
    //    val strikeHigh = forward*2.9
    //    val strikeSize = 100

    val strike = forward

    val FmaxTruncation = 5
    //    for (i <- 0 until strikeSize) {

    //    val strike = strikeLow + (strikeHigh - strikeLow)/strikeSize*i
    val refPDE = new HaganSABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 64, FmaxTruncation)
    val refPrice = refPDE.price(true, strike)
    //1280 * 4, 1280 * 64, 5 => 0.14962173392511555
    //    val refPrice = 0.14962173392511555
    val refVol = solver.impliedVolatility(isCall, strike, refPrice, forward, df, tte)
    println(refPDE.computeCourantNumber + " " + refPrice + " " + refVol)
    for (i <- 0 until 5) {
      println("i=" + i)
      //println("SpaceSteps TimeSteps Scheme Price PriceError Vol Error Time")
      println("SpaceSteps & TimeSteps & Scheme & Price & PriceError & Vol &  Error & Time\\\\")
      val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
      val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
      for (spaceStep <- spaceSteps) {
        var buffer = ""
        for (timeStep <- timeSteps) {
          var startTime = System.nanoTime()
          val cn = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          cn.useSmoothing = false
          cn.useRannacher = false
          var price = cn.price(isCall, strike)
          var endTime = System.nanoTime()
          var elapsed = (endTime - startTime) * 1e-9
          var error = math.abs(price - refPrice)
          var vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          var volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep CN $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & CN & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          //          startTime = System.nanoTime()
          //          var scn = new HaganBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          //          price = scn.price(isCall, strike)
          //          endTime = System.nanoTime()
          //          elapsed = (endTime - startTime) * 1e-9
          //          error = math.abs(price - refPrice)
          //          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          //          volError = math.abs(vol - refVol)
          //          buffer += f"$spaceStep $timeStep BDF2 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          startTime = System.nanoTime()
          var ran = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          ran.useRannacher = true
          price = ran.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep RAN $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & RAN & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          startTime = System.nanoTime()
          val lmg2 = new HaganLMG2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = lmg2.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep LMG2 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & LMG2 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          startTime = System.nanoTime()
          val lmg3 = new HaganLMG3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = lmg3.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep LMG3 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & LMG3 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"

          //          startTime = System.nanoTime()
          //          val bdf3 = new HaganBDF3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          //          price = bdf3.price(isCall, strike)
          //          endTime = System.nanoTime()
          //          elapsed = (endTime - startTime) * 1e-9
          //          error = math.abs(price - refPrice)
          //          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          //          volError = math.abs(vol - refVol)
          //          buffer += f"$spaceStep $timeStep BDF3 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //          //         buffer += f"$spaceStep & $timeStep & LMG3 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"

          startTime = System.nanoTime()
          val trbdf2 = new HaganTRBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = trbdf2.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep TRBDF2 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & LMG3 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          startTime = System.nanoTime()
          val re = new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = re.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep RE $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & RE & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
        }
        print(buffer)
      }
      //    }
    }
  }

  test("hagan-vol-simple") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 40
    val FmaxTruncation = 5 * forward
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    var pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = false
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = true
    pdeMap += "RAN" -> pde
    pdeMap += "LMG2" -> new HaganLMG2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LMG3" -> new HaganLMG3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)

    for (pde <- pdeMap.values) {
      pde.solve()
      println(pde.computeCourantNumber + " " + pde.h)
    }
    println("vols")
    println("Scheme Strike Vol")
    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 10.0)
      val isCall = strike >= forward
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      println(f"Hagan $strike%2.4f $haganVol%2.3f")
      //      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      //      val normalVolScaled = normalVol * math.log(forward / strike) / (forward - strike)
      //      var normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      //      normalVol = solver.impliedVolatility(isCall, strike, normalPrice, forward, df, tte)
      //      val ahVolApprox = SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, strike, tte)
      //      val benPrice = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
      //      var benVol = Double.NaN
      //      try {
      //        benVol = solver.impliedVolatility(isCall, strike, benPrice, forward, df, tte)
      //      } catch {
      //        case e: RuntimeException =>
      //      }
      for ((name, pde) <- pdeMap) {
        val pdePrice = pde.price(isCall, strike)
        var pdeVol = Double.NaN
        try {
          pdeVol = solver.impliedVolatility(isCall, strike, pdePrice, forward, df, tte)
        } catch {
          case e: RuntimeException =>
        }
        println(f"$name $strike%2.4f $pdeVol%2.3f")
      }
    }

  }

  test("hagan-key-numbers") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 5
    val FmaxTruncation = 5 * forward
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = false
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = true
    pdeMap += "RAN" -> pde
    pdeMap += "LMG2" -> new HaganLMG2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LMG3" -> new HaganLMG3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    println("Scheme & h & Q(f) & QL & QR\\\\")
    for ((name, pde) <- pdeMap) {
      pde.solve()
      val h = pde.h
      val j0 = pde.indexForward
      val Qforward = pde.Q0(j0)
      val QL = pde.QL
      val QR = pde.QR
      val dt = pde.dt
      println(f"$name & $h%2.12f & $dt%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")

      //      println(f"$name & $h%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")
    }

  }
  test("hagan-density-simple") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 5
    val FmaxTruncation = 5 * forward
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    var pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = false
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = true
    pdeMap += "RAN" -> pde
    pdeMap += "DOUGLAS" -> new HaganDouglasSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LMG2" -> new HaganLMG2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LMG3" -> new HaganLMG3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)

    for ((name, pde) <- pdeMap) {
      pde.solve()
      println(name + " " + pde.computeCourantNumber + " " + pde.h)
    }
    println("densities")
    //density
    println("Scheme Strike Density")
    val strikeList = pde.F
    for (i <- 1 until strikeList.length) {
      val strike = strikeList(i)
      val eps = 1e-3 * strike
      val isCall = strike >= forward
      val haganVolUp = SABRVanilla.impliedVolatilityHagan(spec, forward, strike + eps, tte)
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      val haganVolDown = SABRVanilla.impliedVolatilityHagan(spec, forward, strike - eps, tte)
      var variance = haganVolUp * haganVolUp * tte
      val callUp = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike + eps, forward, variance, 1.0, df)
      variance = haganVol * haganVol * tte
      val call = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, variance, 1.0, df)
      variance = haganVolDown * haganVolDown * tte
      val callDown = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike - eps, forward, variance, 1.0, df)
      val densityHagan = (callUp - 2 * call + callDown) / (eps * eps)

      println(f"Hagan $strike%2.4f $densityHagan%2.2e")

      var normalVolUp = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike + eps, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      var normalVolDown = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike - eps, tte)
      val callUpNormal = BachelierVanillaEuropean.price(isCall, strike + eps, forward, normalVolUp, tte)
      val callNormal = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      val callDownNormal = BachelierVanillaEuropean.price(isCall, strike - eps, forward, normalVolDown, tte)
      val densityNormal = (callUpNormal - 2 * callNormal + callDownNormal) / (eps * eps)

      for ((name, pde) <- pdeMap) {
        val Q = pde.Q0(i)

        println(f"$name $strike%2.4f $Q%2.2e")
      }
    }
  }

  test("pde-performance-andreasen") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val isCall = false
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val strike = forward * 1.0
    val FmaxTruncation = 5 * forward
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 8, FmaxTruncation)
    //    val refPrice = refPDE.price(isCall, strike)
    //1280 * 4, 1280 * 64, 5 => 0.007979119446512159 
    val refPrice = 0.007979119446512159
    val refVol = solver.impliedVolatility(isCall, strike, refPrice, forward, df, tte)
    println(refPDE.computeCourantNumber + " " + refPrice + " " + refVol)
    for (i <- 0 until 5) {
      println("i=" + i)
      println("SpaceSteps & TimeSteps & Scheme & Price & PriceError & Vol &  Error & Time\\\\")
      val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
      val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
      for (spaceStep <- spaceSteps) {
        var buffer = ""
        for (timeStep <- timeSteps) {
          var startTime = System.nanoTime()
          val cn = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          cn.useSmoothing = false
          cn.useRannacher = false
          var price = cn.price(isCall, strike)
          var endTime = System.nanoTime()
          var elapsed = (endTime - startTime) * 1e-9
          var error = math.abs(price - refPrice)
          var vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          var volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep CN $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & CN & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          //          startTime = System.nanoTime()
          //          var scn = new HaganBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          //          price = scn.price(isCall, strike)
          //          endTime = System.nanoTime()
          //          elapsed = (endTime - startTime) * 1e-9
          //          error = math.abs(price - refPrice)
          //          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          //          volError = math.abs(vol - refVol)
          //          buffer += f"$spaceStep $timeStep BDF2 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          startTime = System.nanoTime()
          var ran = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          ran.useRannacher = true
          price = ran.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep RAN $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & RAN & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          startTime = System.nanoTime()
          val lmg2 = new HaganLMG2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = lmg2.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep LMG2 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & LMG2 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          startTime = System.nanoTime()
          val lmg3 = new HaganLMG3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = lmg3.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep LMG3 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & LMG3 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          //          startTime = System.nanoTime()
          //          val bdf3 = new HaganBDF3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          //          price = bdf3.price(isCall, strike)
          //          endTime = System.nanoTime()
          //          elapsed = (endTime - startTime) * 1e-9
          //          error = math.abs(price - refPrice)
          //          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          //          volError = math.abs(vol - refVol)
          //          buffer += f"$spaceStep $timeStep BDF3 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //          //         buffer += f"$spaceStep & $timeStep & LMG3 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"

          startTime = System.nanoTime()
          val trbdf2 = new HaganTRBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = trbdf2.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep TRBDF2 $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & LMG3 & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
          startTime = System.nanoTime()
          val re = new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          price = re.price(isCall, strike)
          endTime = System.nanoTime()
          elapsed = (endTime - startTime) * 1e-9
          error = math.abs(price - refPrice)
          vol = solver.impliedVolatility(isCall, strike, price, forward, df, tte)
          volError = math.abs(vol - refVol)
          buffer += f"$spaceStep $timeStep RE $price%2.4e $error%2.1e $vol%2.5f $volError%2.1e  $elapsed%2.1e\n"
          //         buffer += f"$spaceStep & $timeStep & RE & $price%2.4e & $error%2.1e & $vol%2.5f & $volError%2.1e & $elapsed%2.1e\\\\\n"
        }
        print(buffer)
      }
    }
  }

  def createHaganSABRSolver(name: String, spec: SABRModelSpec, forward: Double, tte: Double, spaceStep: Int = 1000, timeStep: Int = 1000 / 5, FmaxTruncation: Double = 5.0): HaganSABRDensitySolver = {
    var pde: HaganSABRDensitySolver = null
    name match {
      case "CN" => {
        pde = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
        pde.useRannacher = false
        pde.useSmoothing = false
      }
      case "RAN" => {
        pde = new HaganSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
        pde.useRannacher = true
        pde.useSmoothing = false
      }
      case "LMG2" => pde = new HaganLMG2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "LMG3" => pde = new HaganLMG3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "LS" => pde = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "TRBDF2" => pde = new HaganTRBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "CNBDF2" => pde = new HaganCNBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "TRBDF3" | "CNBDF3" => pde = new HaganCNBDF3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "RE" => pde = new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "BDF2" => pde = new HaganBDF2SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "BDF3" => pde = new HaganBDF3SABRDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
    }
    return pde
  }


  def computeMaxError(forward: Double, F: Array[Double], Q: Array[Double], Qold: Array[Double]): (Double, Double, Double) = {
    var maxError = 0.0
    var Fmax = 0.0
    var value = 0.0
    for (i <- 1 until F.size - 2) {
      if (Qold != null) {
        var error = Q(i) - Qold(i)
        if (math.abs(error) > maxError) {
          maxError = error
          Fmax = F(i)
        }
      }
      if (math.abs(F(i) - forward) < Epsilon.MACHINE_EPSILON_SQRT) {
        value = Q(i)
      }
    }
    return (value, maxError, Fmax)
  }

  test("pde-density-convergence-hagan") {
    val isLatex = false
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "CNBDF3", "RE")
    val spaceSteps = Array[Int](500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
    for (l <- 0 to 8) {
      println("TimeSteps & ATM Value & Max Difference & Ratio & Time(s)\\\\")
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {
          var buffer = ""
          println("\\hline \\multicolumn{5}{|c|}{" + scheme + "} \\\\ \\hline")
          var previousQ: Array[Double] = null
          var previousDiff = Double.NaN
          for (timeStep <- timeSteps) {
            var startTime = System.nanoTime()
            val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
            cn.solve()
            var endTime = System.nanoTime()
            var elapsed = (endTime - startTime) * 1e-9
            var err = computeMaxError(forward, cn.F, cn.Q0, previousQ)
            val value: Double = err._1
            val diff: Double = err._2
            val ratio = previousDiff / diff
            buffer += f" $timeStep & $value%2.8f & $diff%2.1e & $ratio%2.1f & $elapsed%2.1e\\\\\n"
            previousQ = cn.Q0
            previousDiff = diff

          }
          print(buffer)

        }
      }
    }
  }

  test("pde-vol-convergence-hagan") {
    val isLatex = false
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "CNBDF3", "RE")
    val spaceSteps = Array[Int](500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
    for (l <- 0 to 8) {
      println("TimeSteps & ATM Value & Max Difference & Ratio & Time(s)\\\\")
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {
          var buffer = ""
          println("\\hline \\multicolumn{5}{|c|}{" + scheme + "} \\\\ \\hline")
          var previousValue = Double.NaN
          var previousDiff = Double.NaN
          for (timeStep <- timeSteps) {
            var startTime = System.nanoTime()
            val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
            val price = cn.price(isCall, forward)
            var endTime = System.nanoTime()
            var elapsed = (endTime - startTime) * 1e-9
            val value = solver.impliedVolatility(isCall, forward, price, forward, df, tte)
            val diff = value - previousValue
            val ratio = previousDiff / diff
            buffer += f" $timeStep & $value%2.8f & $diff%2.1e & $ratio%2.1f & $elapsed%2.1e\\\\\n"
            previousValue = value
            previousDiff = diff

          }
          print(buffer)

        }
      }
    }
  }

}
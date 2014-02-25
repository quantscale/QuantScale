package test.quantscale.analytic

import org.junit.runner.RunWith
import quantscale.analytic._
import quantscale.fdm.mesh._
import quantscale.fdm.sabr._
import quantscale.math.CubicPP
import quantscale.math.CubicSpline
import quantscale.vol.{Li2011ImpliedVolatilityGuess, Li2011SORBlackVolatilitySolver, Li2011SORDRBlackVolatilitySolver, Li2011SORTSBlackVolatilitySolver}
import org.scalatest.junit.{JUnitSuite, JUnitRunner}
import quantscale.fdm._
import _root_.org.junit.Test
import org.scalatest.{FunSuite, Suite}
import quantscale.fdm.sabr.RannacherSmoothing
import quantscale.fdm.payoff.{SimpsonIntegralSmoother, VanillaFDPayoff}
import quantscale.fdm.method._
import quantscale.fdm.sabr.RannacherSmoothing
import scala.Array
import quantscale.fdm.sabr.RannacherSmoothing
import quantscale.fdm.sabr.RannacherSmoothing

@RunWith(classOf[JUnitRunner])
class SABRSuite extends FunSuite {
  test("PutCallParity") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 10
    val nDeviation = 3.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.useRannacher = false
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.useRannacher = true
    pdeMap += "RAN" -> pde

    pdeMap += "LMG2" -> new HaganLMG2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LMG3" -> new HaganLMG3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF3" -> new HaganCNBDF3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "BDF2" -> new HaganBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "BDF3" -> new HaganBDF3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)

    for ((name, pde) <- pdeMap) {
      pde.solve()
      val pricePut = pde.price(false, forward)
      val priceCall = pde.price(true, forward)
      println(name + " " + pricePut + " " + priceCall)
      assert(pricePut - priceCall < 1e-8, "put=" + pricePut + " call=" + priceCall)

      // println(f"$name & $h%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")
    }
  }

  @Test def testLiSorDRIV {
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

  @Test def testLiSORTSIV() {
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

  @Test def testLiSORIV() {
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

  @Test def testLiHagan() {
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

  @Test def testHagan() {
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

  test("HaganFDPrice") {
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

  test("localvol-trbdf2") {
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

    var spaceSize: Int = 400;
    var timeSize: Int = 100;
    val mu = 0.0;
    val r = 0.0;
    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      new MeshBoundaries(0, tte, 0, 4 * forward),
      strike);

    val fdSpec = new SABRLocalVolFDSpec(
      grid,
      spec,
      forward, tte,
      mu,
      r)

    val payoff = new VanillaFDPayoff(isCall, strike, tte);
    val solver = new ThomasTridiagonalSolver();
    val method = new TRBDF2Parabolic1DMethod(payoff);
    val pricer = new FDMSolver1D(fdSpec, method, solver);
    pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricer.solve(payoff);

    val price = pricer.price(forward);
    println("local vol price=" + price);

    val adPricer = new DensityGaussianApproxPricer(spec, forward, tte, 100, 5 * forward)
    val adPrice = adPricer.price(isCall, strike)
    println("density approx price=" + adPrice);


  }


  test("HaganBenhamou") {
    //    val alpha = 0.0758194;
    //    val nu = 0.1;
    //    val beta = 0.5;
    //    val rho = -0.1;
    //    val forward = 0.02;
    //    val tte = 2.0;
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.0325;
    val tte = 10.0;
    val strikes = Array(0.01, 0.02, 0.025, 0.03, 0.04)
    for (strike <- strikes) {
      val spec = new SABRModelSpec(alpha, beta, nu, rho)
      val blackVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      //println(blackVol)
      val isCall = true
      val blackPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, forward, blackVol * blackVol * tte, 1.0, 1.0)
      val benPrice = SABRVanilla.priceBenhamou(spec, isCall, strike, forward, tte)
      println("BEN " + blackPrice + " " + benPrice + " " + math.abs(blackPrice - benPrice))
      val normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      val normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      println("HAG " + blackPrice + " " + normalPrice + " " + (blackPrice - normalPrice))
      val ahPrice = SABRVanilla.priceAndreasenHuge(isCall, strike, spec, forward, tte)
      println("AHS " + blackPrice + " " + ahPrice + " " + (blackPrice - ahPrice))
      val ahzPrice = SABRVanilla.priceAndreasenHugeZABR(isCall, strike, new ZABRModelSpec(alpha, beta, nu, rho, 1.0), forward, tte)
      println("AHZ " + blackPrice + " " + ahzPrice + " " + (blackPrice - ahzPrice))
    }

  }

  @Test def testHaganNormal() {
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

  @Test def testAndreasenHuge() {
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

  @Test def testAndreasenHugePrice() {
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

  test("AndreasenHugeDensity") {
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
    //val pde = new HaganTruncatedSABRDensitySolver(spec, forward, tte, 50, 10)
    val pde = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, 50, 10, 5.0 * forward)

    var spaceSize: Int = 50;
    var timeSize: Int = 10;
    val mu = 0.0;
    val r = 0.0;
    val fdsolver = new ThomasTridiagonalSolver();
    val payoffCall = new VanillaFDPayoff(false, forward, tte);
    val payoffPut = new VanillaFDPayoff(true, forward, tte);
    val timeMesh = new UniformMesh1D(timeSize, new Mesh1DBoundaries(0, tte))
    val spaceBoundaries = new Mesh1DBoundaries(0, 5 * forward);

    val methodCall = new TRBDF2Parabolic1DMethod(payoffCall);
    //     method.lowerBoundary  = ForwardPartialOrder2Parabolic1DBoundaryFactory
    //     method.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
    val grid =
      new StaticAdaptiveMesh2D(new UniformMesh1D(spaceSize, spaceBoundaries, forward, true), timeMesh)

    val fdSpec = new SABRLocalVolFDSpec(grid, spec, forward, tte, mu, r)

    val pricerCall = new FDMSolver1D(fdSpec, methodCall, fdsolver);
    // pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricerCall.solve(payoffCall);
    val callLVSplineCall = CubicSpline.makeCubicSpline(grid.spaceVector, pricerCall.price)
    val methodPut = new TRBDF2Parabolic1DMethod(payoffPut);
    val pricerPut = new FDMSolver1D(fdSpec, methodPut, fdsolver);
    // pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricerPut.solve(payoffPut);
    val callLVSplinePut = CubicSpline.makeCubicSpline(grid.spaceVector, pricerPut.price)

    println("vols")
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeApprox AndreasenHugeRaw AndreasenHuge HaganPDE HaganLocalVol")
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
      val splineLV = if (isCall) callLVSplineCall else callLVSplinePut
      val lvPrice = splineLV.value(strike)
      var lvVol = Double.NaN
      try {
        lvVol = solver.impliedVolatility(isCall, strike, lvPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      println(f"$strike%2.4f $haganVol%2.3f $normalVol%2.3f $normalVolScaled%2.3f $ahVolApprox%2.3f $ahVolRaw%2.3f $ahVol%2.3f $pdeVol%2.3f $lvVol%2.3f")
    }
    println("densities")
    //density
    // println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeRaw AndreasenHuge Benhamou")
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeRaw AndreasenHuge HaganPDE HaganLocalVol")

    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 20.0)
      val eps = 1e-2 * strike
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
      val splineLV = if (isCall) callLVSplineCall else callLVSplinePut
      val densityLV = splineLV.secondDerivative(strike);
      println(f"$strike%2.4f $densityHagan%2.2e $densityNormal%2.2e $densityNormalS%2.2e $densityAhRaw%2.2e $densityAh%2.2e $densityPDE%2.2e $densityLV%2.2e")
      //       println(f"$strike%2.4f $call%2.2e $callNormal%2.2e $callBen%2.2e")
    }
  }

  test("HaganDensity") {
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
    val pde = new HaganSABRDensitySolver(spec, forward, tte, 500, 10, 5.0)
    pde.useSmoothing = false
    pde.useRannacher = false
    pde.solve()
    println(pde.computeCourantNumber + " " + pde.h)
    val pde2 = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, 500, 10, 5.0)

    var spaceSize: Int = 500;
    var timeSize: Int = 10;
    val mu = 0.0;
    val r = 0.0;
    val fdsolver = new ThomasTridiagonalSolver();
    val payoffCall = new VanillaFDPayoff(false, forward, tte);
    val payoffPut = new VanillaFDPayoff(true, forward, tte);
    val timeMesh = new UniformMesh1D(timeSize, new Mesh1DBoundaries(0, tte))
    val spaceBoundaries = new Mesh1DBoundaries(0, 5 * forward);

    val methodCall = new TRBDF2Parabolic1DMethod(payoffCall);
    methodCall.lowerBoundary = ForwardPartialOrder2Parabolic1DBoundaryFactory
    methodCall.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
    val grid =
      new StaticAdaptiveMesh2D(new UniformMesh1D(spaceSize, spaceBoundaries, forward, true), timeMesh)

    val fdSpec = new SABRLocalVolFDSpec(grid, spec, forward, tte, mu, r)

    val pricerCall = new FDMSolver1D(fdSpec, methodCall, fdsolver);
    // pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricerCall.solve(payoffCall);
    val callLVSplineCall = CubicSpline.makeBesselSpline(grid.spaceVector, pricerCall.price)
    val methodPut = new TRBDF2Parabolic1DMethod(payoffPut);
    methodPut.lowerBoundary = ForwardPartialOrder2Parabolic1DBoundaryFactory
    methodPut.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
    val pricerPut = new FDMSolver1D(fdSpec, methodPut, fdsolver);
    // pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricerPut.solve(payoffPut);
    val callLVSplinePut = CubicSpline.makeBesselSpline(grid.spaceVector, pricerPut.price)
    val adPricer = new DensityGaussianApproxPricer(spec, forward, tte, 1000, 5 * forward)

    println("vols")
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeApprox AndreasenHugeRaw AndreasenHuge HaganPDE HaganEulerPDE HaganLocalVol")
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
      val splineLV = if (isCall) callLVSplineCall else callLVSplinePut
      val lvPrice = splineLV.value(strike)
      val callPrice = callLVSplineCall.value(strike)
      val putPrice = callLVSplinePut.value(strike)
      var lvVol = Double.NaN
      try {
        lvVol = solver.impliedVolatility(isCall, strike, lvPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      val adPrice = adPricer.price(isCall, strike)
      var adVol = Double.NaN
      try {
        adVol = solver.impliedVolatility(isCall, strike, adPrice, forward, df, tte)
      } catch {
        case e: RuntimeException =>
      }
      println(f"$strike%2.4f $haganVol%2.3f $normalVol%2.3f $normalVolScaled%2.3f $ahVolApprox%2.3f $ahVolRaw%2.3f $ahVol%2.3f $pdeVol%2.3f $pde2Vol%2.3f $lvVol%2.3f $adVol%2.3f")
    }
    println("densities")
    //density
    println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeRaw AndreasenHuge HaganPDE HaganEulerPDE HaganLocalVol")
    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 20.0)
      val eps = 1e-2 * strike
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
      val splineLV = if (isCall) callLVSplineCall else callLVSplinePut
      val densityLV = splineLV.secondDerivative(strike);

      val adPriceUp = adPricer.price(isCall, strike + eps)
      val adPrice = adPricer.price(isCall, strike)
      val adPriceDown = adPricer.price(isCall, strike - eps)

      val adDensity = (adPriceUp - 2 * adPrice + adPriceDown) / (eps * eps)
      println(f"$strike%2.4f $densityHagan%2.2e $densityNormal%2.2e $densityNormalS%2.2e $densityAhRaw%2.2e $densityAh%2.2e $densityPDE%2.2e $densityPDE2%2.2e $densityLV%2.2e $adDensity%2.2e")
    }

  }

  test("Alan") {
    val forward = 1.0
    val K = 0.8//0.8
    val tte = 10.0
    val nu = 1.00
    val beta = 0.3//0.30
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
    val pdeTransformed = new HaganSABRTransformedDensitySolver(spec, forward, tte, 1000, 100, 4)
    pdeTransformed.smoothing = new RannacherSmoothing()
    val volHagan =  SABRVanilla.impliedVolatilityHagan(spec, forward, K, tte)
    val volAHShort =  SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, K, tte)
    val priceAH = splinePut.value(K)
    val priceAHRaw = splinePutRaw.value(K)
    val pricePDE = pde.price(false, K)
    val pricePDET = pdeTransformed.price(false, K)
    val normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, K, tte)
    val normalPrice = BachelierVanillaEuropean.price(false, K, forward, normalVol, tte)
    //    val volNormal = 100*solver.impliedVolatility(true, K, normalPrice, forward, df, tte)
    val volAH =  solver.impliedVolatility(false, K, priceAH, forward, df, tte)
    val volAHRaw =  solver.impliedVolatility(false, K, priceAHRaw, forward, df, tte)
    val volPDE =  solver.impliedVolatility(false, K, pricePDE, forward, df, tte)
    val volPDET =  solver.impliedVolatility(false, K, pricePDET, forward, df, tte)
    val priceBen = SABRVanilla.priceBenhamou(spec, false, K, forward, tte) + K - forward
    println("BS      " + BlackScholesVanillaEuropean.priceEuropeanVanilla(false, K, forward, volHagan * volHagan * tte, 1.0, df))
    println("AHShort " + BlackScholesVanillaEuropean.priceEuropeanVanilla(false, K, forward, volAHShort * volAHShort * tte, 1.0, df))
    println("AH      "+priceAH)
    println("PDE     " + pricePDE)
    println("PDET    " + pricePDET)
    println(priceBen)
    println(normalPrice)
    val volBen = Double.NaN //100* solver.impliedVolatility(false, K, priceBen, forward, df, tte)

    println(f"$volHagan%2.3f $volBen%2.3f $volAHShort%2.3f $volAH%2.3f $volPDE%2.3f $volPDET%2.3f")
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
      var error = (Q0(i) - Qref)
      if (error < 0) error = -error
      if (error > maxError) {
        maxError = error
        Fmax = F(i)
      }
    }
    return (maxError, Fmax)
  }

  @Test def testPDEDensityAccuracyAH() {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0

    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 5.0 * forward

    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 500, 1280 * 4, FmaxTruncation)

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


  test("AndreasenHugeGamma") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.0325;
    val tte = 10.0;
    val gamma = Array(0.0, 0.5, 1.0, 1.5, 1.7)
    //normal vol
    for (g <- gamma) {
      val spec = new ZABRModelSpec(alpha, beta, nu, rho, g)

      val strikePriceVol = SABRVanilla.priceAndreasenHugeZABRStrikePriceVol(true, spec, forward, tte, true)
      val strikes = strikePriceVol._1
      val xs = strikePriceVol._3
      val vols = Array.ofDim[Double](xs.length)
      for (i <- 0 until strikes.length) {
        val normalvol = if (math.abs(forward - strikes(i)) < 1e-6) {
          alpha * math.pow(forward, beta) //   alpha*math.pow(forward, beta-1)
        } else {
          (forward - strikes(i)) / xs(i) //(forward-strikes(i))/xs(i)
        }
        val isCall = forward <= strikes(i)
        val price = BachelierVanillaEuropean.price(isCall, strikes(i), forward, normalvol, tte)
        val blackvol = new Li2011SORTSBlackVolatilitySolver().impliedVolatility(isCall, strikes(i), price, forward, 1.0, tte)
        vols(i) = blackvol
        println("Normal " + g + " " + strikes(i) + " " + vols(i))
      }
    }
    //density
    //    for (g <- gamma) {
    //      val spec = new ZABRModelSpec(alpha, beta, nu, rho, g)
    //
    //      val strikePriceVol = SABRVanilla.priceAndreasenHugeZABRStrikePriceVol(true, spec, forward, tte, true)
    //      val strikes = strikePriceVol._1
    //      val xs = strikePriceVol._3
    //      val vols = Array.ofDim[Double](xs.length)
    //      for (i <- 0 until strikes.length) {
    //        if (math.abs(forward-strikes(i))<1e-6) {
    //           vols(i) = alpha*math.pow(forward, beta-1)
    //        } else {
    //        vols(i) = Math.log(forward/strikes(i))/xs(i)
    //        }
    //      }
    //      var call0 = BlackScholesVanillaEuropean.priceEuropeanVanilla(true, strikes(0), forward, vols(0)*vols(0)*tte, 1.0, 1.0)
    //      var call1 = BlackScholesVanillaEuropean.priceEuropeanVanilla(true, strikes(1), forward, vols(1)*vols(1)*tte, 1.0, 1.0)
    //
    //      for (i <- 2 until strikes.length) {
    //        val call2 = BlackScholesVanillaEuropean.priceEuropeanVanilla(true, strikes(i), forward, vols(i)*vols(i)*tte, 1.0, 1.0)
    //        val density = (call2 - 2*call1 +call0)/(strikes(i)-strikes(i-1))/(strikes(i-1)-strikes(i-2))
    //        println(g+" "+strikes(i-1)+" "+density)
    //        call0 = call1
    //        call1 = call2
    //      }
    //    }

    for (g <- gamma) {
      val spec = new ZABRModelSpec(alpha, beta, nu, rho, g)

      val splineCall = SABRVanilla.priceAndreasenHugeZABR(true, spec, forward, tte, false)
      val splinePut = SABRVanilla.priceAndreasenHugeZABR(false, spec, forward, tte, false)
      for (i <- 0 to 1950) {
        val strike = 0.005 + i * 1e-4
        val isCall = strike >= forward
        val ahPrice = if (isCall) splineCall.value(strike) else splinePut.value(strike)
        var solver = new Li2011SORDRBlackVolatilitySolver(1e-8)
        val ahVol = solver.impliedVolatility(isCall, strike, ahPrice, forward, 1.0, tte)
        println(g + " " + strike + " " + ahPrice + " " + ahVol)
      }
    }

    //    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    //    val splineCall = SABRVanilla.priceAndreasenHuge(true, spec, forward, tte, true)
    //    val splinePut = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte, true)
    //    for (i <- 1 to 100) {
    //      val strike = 0.0 + i*0.2/100
    //      val isCall = strike >= forward
    //      val ahPrice = if (isCall) splineCall.value(strike) else splinePut.value(strike)
    //      var solver = new Li2011SORDRBlackVolatilitySolver(1e-8)
    //      val ahVol = solver.impliedVolatility(isCall, strike, ahPrice, forward, 1.0, tte)
    //      println("SAB "+strike+" "+ahPrice+" " + ahVol)
    //    }
  }

  @Test def testPDEVolAccuracyAH() {
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

    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "TRBDF3")
    for (l <- 0 to 7) {
      println("SpaceSteps TimeSteps Scheme MaxError Time")
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {
          var buffer = ""
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

  private def createErrorString(isLatex: Boolean, scheme: String, spaceStep: Int, timeStep: Int, error: Double, errorLocation: Double, elapsed: Double): String = {
    return if (isLatex) f"$spaceStep & $timeStep & $scheme & $error%2.1e & $elapsed%2.1e\\\\\n"
    else f"$spaceStep $timeStep $scheme $error%2.1e $errorLocation%2.4f $elapsed%2.1e\n"
  }

  test("PDEDensityAccuracyHagan") {
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
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    val Qspline = CubicSpline.makeCubicSpline(refPDE.F, refPDE.Q0)

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "BDF2", "BDF3", "TRBDF3", "RE")
    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 40, 80, 160, 320, 640, 1280)
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

  @Test def testPDEVolAccuracyHagan() {
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
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 40, 80, 160, 320, 640, 1280)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "TRBDF3", "RE")
    for (l <- 0 until 9) {
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

  @Test def testPDEvsFormulaPerformanceHagan() {
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

  @Test def testPDEvsFormulaPerformanceAH() {
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

  @Test def testPDEPerformanceHagan() {
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

  @Test def testHaganVolSimple() {
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

  @Test def testHaganKeyNumbers() {
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
    val nDeviation = 3.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = false
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useSmoothing = false
    pde.useRannacher = true
    pdeMap += "RAN" -> pde


    //pdeMap += "LMG2" -> new HaganLMG2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "LMG3" -> new HaganLMG3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "TRBDF2" -> new HaganTRBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "RE" -> new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "TRBDF3" -> new HaganCNBDF3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    val h = pde.h
    val dt = pde.dt
    println(f"h=$h%2.12f")
    println(f"dt=$dt%2.12f")
    println("Scheme & ATM price & Q(f) & QL & QR\\\\")
    for (i <- 0 until 100) {
      val startTime = System.nanoTime()
      for ((name, pde) <- pdeMap) {
        pde.solve()
        val j0 = pde.indexForward
        val Qforward = pde.Q0(j0)
        val QL = pde.QL
        val QR = pde.QR

        val price = pde.price(false, forward)
        println(f"$name & $price%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")

        //      println(f"$name & $h%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")
      }
      println("elapsed " + (System.nanoTime() - startTime) * 1e-9)
    }

  }

  @Test def testHaganDensitySimple() {
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

  @Test def testPDEPerfomanceAndreason() {
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
        val error = Q(i) - Qold(i)
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

  @Test def testPDEDensityConvergenceHagan() {
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

  @Test def testPDEVolConvergenceHagan() {
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
package quantscale.analytic

import quantscale.fdm.sabr._
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.scalatest.{Tag, FunSuite}
import quantscale.fdm.mesh._
import quantscale.vol.Li2011SORDRBlackVolatilitySolver
import quantscale.fdm.{FDMSolver1D, SABRLocalVolFDSpec, ThomasTridiagonalSolver, Epsilon}
import quantscale.math.{CubicPP, CubicSpline}
import org.junit.Test
import quantscale.fdm.payoff.VanillaFDPayoff
import quantscale.fdm.method._
import quantscale.fdm.sabr.RannacherSmoothing
import quantscale.fdm.sabr.BDF2Smoothing
import quantscale.fdm.sabr.TRBDF2BDF3Smoothing
import quantscale.fdm.sabr.TRBDF2BDF3Smoothing
import quantscale.fdm.sabr.RannacherSmoothing
import quantscale.fdm.sabr.BDF2Smoothing

object PerformanceTest extends Tag("PerformanceTest")

//in SBT, test-only quantscale.analytic.TransformedSABRSuite -- -n "PerformanceTest"
@RunWith(classOf[JUnitRunner])
class TransformedSABRSuite extends FunSuite {

  test("Alan") {
    val forward = 1.0
    val K = 0.8
    val tte = 10.0
    val nu = 1.00
    val beta = 0.30
    val alpha = 1.00
    val rho = 0.90
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val splinePut = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte)
    val splinePutRaw = SABRVanilla.priceAndreasenHuge(false, spec, forward, tte, false)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val volHagan = 100 * SABRVanilla.impliedVolatilityHagan(spec, forward, K, tte)
    val volAHShort = 100 * SABRVanilla.impliedVolatilityAndreasenHuge(spec, forward, K, tte)
    val priceAH = splinePut.value(K)
    val priceAHRaw = splinePutRaw.value(K)

    val normalVol = 100 * SABRVanilla.normalVolatilityHagan2013(spec, forward, K, tte)
    val normalPrice = BachelierVanillaEuropean.price(false, K, forward, normalVol / 100, tte)
    val normalPriceParity = K - forward + BachelierVanillaEuropean.price(true, K, forward, normalVol / 100, tte)
    //    val volNormal = 100*solver.impliedVolatility(true, K, normalPrice, forward, df, tte)
    val volAH = 100 * solver.impliedVolatility(false, K, priceAH, forward, df, tte)
    val volAHRaw = 100 * solver.impliedVolatility(false, K, priceAHRaw, forward, df, tte)

    val lognormalPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(false, K, forward, volHagan / 100 * volHagan / 100 * tte, 1.0, df)
    println(f"HaganLN $lognormalPrice%2.3f $volHagan%2.3f")
    println(f"HaganN $normalPrice%2.3f $normalVol%2.3f")
    println(normalPriceParity)
    println(f"priceAH $priceAH%2.3f $volAH%2.3f")

    val Fmax = Array(5, 50, 500, 5000, 10000)
    val ndev = Array(3, 4, 5, 10, 100)
    val size = Array(10, 20, 100, 1000, 10000)
    val timesteps = Array(5, 10, 20, 40, 160)

    //    val Fmax = Array(50)
    //    val ndev = Array(10)
    //    val size = Array(10)
    //    val timesteps = Array(5)

    for (tsteps <- timesteps) {
      for (xsteps <- size) {
        for (i <- 0 until Fmax.length) {
          val pde = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xsteps, tsteps, Fmax(i))
          val pdeTransformed = new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xsteps, tsteps, ndev(i))
          val pricePDE = pde.price(false, K)
          val pricePDET = pdeTransformed.price(false, K)

          val volPDE = 100 * solver.impliedVolatility(false, K, pricePDE, forward, df, tte)
          val volPDET = 100 * solver.impliedVolatility(false, K, pricePDET, forward, df, tte)

          print(Fmax(i) + " & " + xsteps + " & " + tsteps + f" & $pricePDE%2.5f & $volPDE%2.3f && ")
          println(ndev(i) + " & " + xsteps + " & " + tsteps + f" & $pricePDET%2.5f & $volPDET%2.3f")
        }
      }
    }
  }


  test("F_z") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 50
    val tSteps = 10
    val nDeviation = 4.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.solve()
    val zmin = pde.zmin
    val Fm = pde.Fm
    val h = pde.h
    var j = 0
    val size = Fm.length
    val zm = Array.ofDim[Double](size)
    println("z F")
    while (j < size) {
      zm(j) = zmin + (j - 0.5) * h
      println(zm(j) + " " + Fm(j))
      j += 1
    }
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

    val pdeTPrice = new HaganSABRTransformedDensitySolver(spec, forward, tte).price(isCall, strike)
    println(pdeTPrice)

    assert(math.abs(pdePrice - blackPrice) < 1e-4)
    assert(math.abs(pdeTPrice - blackPrice) < 1e-4)

  }

  test("courant") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 40
    val nDeviation = 4.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.computeCourantNumber();
    var pdeu = new HaganSABRDensitySolver(spec, forward, tte, xSteps, tSteps, 5.0)
    println(pdeu.computeCourantNumber());
  }
  test("put-call-parity") {
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
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.smoothing = new RannacherSmoothing()
    pdeMap += "RAN" -> pde
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.smoothing = new TRBDF2BDF3Smoothing()
    pdeMap += "LEF" -> pde

    pdeMap += "LMG2" -> new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LMG2a" -> new HaganLMG2aSABRTransformedDensitySolver(spec, forward, tte, xSteps, 1e-1, nDeviation)
    pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "BDF2" -> new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    val pde2 = new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde2.useTRinit = true
    pdeMap += "BDF2t" -> pde2
    pdeMap += "BDF3" -> new HaganBDF3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    val pde3 = new HaganBDF3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde3.useTRinit = true
    pdeMap += "BDF3t" -> pde3

    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "TRBDF2a" -> new HaganTRBDF2aSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "RE3" -> new HaganRichardsonEuler3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "REN" -> new HaganRichardsonEulerNSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "Bathe" -> new HaganBathe3SubstepsSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "IMBDF3" -> new HaganIMBDF3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)

    for ((name, pde) <- pdeMap) {
      pde.solve()
      val pricePut = pde.price(false, forward)
      val priceCall = pde.price(true, forward)
      println(name + " " + pricePut + " " + priceCall)
      assert(pricePut - priceCall < 1e-8, "put=" + pricePut + " call=" + priceCall)

      // println(f"$name & $h%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")
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
    val nDeviation = 4.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.smoothing = new RannacherSmoothing()
    pdeMap += "RAN" -> pde
    pdeMap += "BDF2" -> new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LMG2" -> new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    //pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "Bathe" -> new HaganBathe3SubstepsSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    //pdeMap += "TRBDF3" -> new HaganTRBDF3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    val h = pde.h
    val dt = pde.dt
    println(f"h=$h%2.12f")
    println(f"dt=$dt%2.12f")
    println("courant=" + pde.computeCourantNumber())
    println("Scheme & ATM price & Q(f) & QL & QR\\\\")
    for (i <- 0 until 1) {
      val startTime = System.nanoTime()
      for ((name, pde) <- pdeMap) {
        pde.solve()
        val j0 = pde.indexForward
        val Qforward = pde.P(j0)
        val QL = pde.PL
        val QR = pde.PR

        val price = pde.price(true, forward)
        println(f"$name & $price%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")

        // println(f"$name & $h%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")
      }
      println("elapsed " + (System.nanoTime() - startTime) * 1e-9)
    }

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
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    val Qspline = CubicSpline.makeCubicSpline(refPDE.Fm, refPDE.P)

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "BDF2", "LMG2", "LMG2a", "LS", "TRBDF2", "Bathe3", "RE3", "RE")
    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 40, 80, 160, 320, 640, 1280)
    var l = 0
    while (l >= 0) {
      //            println("SpaceSteps TimeSteps Scheme MaxError FError Time")
      println("SpaceSteps & TimeSteps & Scheme & MaxError & Time\\\\")
      var buffer = ""
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {

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

        }

      }
      print(buffer)
      l -= 1
    }
  }

  test("pde-density-accuracy-ah") {
    val isLatex = false
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    val Qspline = CubicSpline.makeCubicSpline(refPDE.Fm, refPDE.P)

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "Bathe", "RE", "RE3")
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

  private def createErrorString(isLatex: Boolean, scheme: String, spaceStep: Int, timeStep: Int, error: Double, errorLocation: Double, elapsed: Double): String = {
    return if (isLatex) f"$spaceStep & $timeStep & $scheme & $error%2.1e & $elapsed%2.1e\\\\\n"
    else f"$spaceStep $timeStep $scheme $error%2.1e $errorLocation%2.4f $elapsed%2.1e\n"
  }

  def computeMaxError(solver: HaganSABRTransformedDensitySolver, Qspline: CubicPP): (Double, Double) = {
    var maxError = 0.0
    var Fmax = 0.0
    val F = solver.Fm
    val Q0 = solver.P

    for (i <- 2 until F.size - 2) {
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

  test("pde-vol-accuracy-hagan") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 40, 80, 160, 320, 640, 1280)
    val schemes = Array("CN", "RAN", "LMG2", "BDF2", "LS", "TRBDF2", "Bathe", "RE")
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

  test("pde-vol-accuracy-ah", Tag("PerformanceTest")) {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 40, 80, 160, 320, 640, 1280)
    val schemes = Array("CN", "RAN", "BDF2", "LMG2", "LS", "TRBDF2", "Bathe", "RE")
    var l = 9
    while (l >= 0) {
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

        }
        print(buffer)
      }
      l -= 1
    }

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
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2.0)).x

    var isCall = false

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "BDF2", "LMG2", "LS", "TRBDF2", "Bathe", "RE")
    val spaceSteps = Array[Int](500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(5, 10, 20, 40, 80, 160, 320, 640)
    for (l <- 0 until 1) {
      println("TimeSteps & ATM Value & Max Difference & Ratio & Time(s)\\\\")
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {
          var buffer = ""
          println("\\hline \\multicolumn{5}{|c|}{" + scheme + "} \\\\ \\hline")
          var previousQ: Array[Double] = Array.ofDim(spaceStep)
          var previousDiff = Double.NaN
          for (timeStep <- timeSteps) {
            val startTime = System.nanoTime()
            val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
            cn.solve()
            val endTime = System.nanoTime()
            val elapsed = (endTime - startTime) * 1e-9
            val err = computeMaxError(forward, cn.Fm, cn.P, previousQ)
            val value: Double = err._1
            val diff: Double = err._2
            val ratio = previousDiff / diff
            buffer += f" $timeStep & $value%2.8f & $diff%2.1e & $ratio%2.1f & $elapsed%2.1e\\\\\n"
            Array.copy(cn.P, 0, previousQ, 0, cn.P.length)
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
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "BDF2", "LMG2", "LMG3", "LS", "TRBDF2", "Bathe", "RE")
    val spaceSteps = //Array(1000, 1000, 1000, 1000, 1000, 1000)
      Array(80, 160, 320, 640, 1280, 1280 * 2)
    val timeSteps = Array(5, 10, 20, 40, 80, 160)
    for (l <- 0 until 10) {
      println("$J$ & $N$ & ATM Value & Max Difference & Ratio & Time(s)\\\\")
      for (scheme <- schemes) {

        var buffer = ""
        println("\\hline \\multicolumn{6}{|c|}{" + scheme + "} \\\\ \\hline")
        var previousValue = Double.NaN
        var previousDiff = Double.NaN
        var index = 0
        for (timeStep <- timeSteps) {
          val spaceStep = spaceSteps(index)
          index += 1
          var startTime = System.nanoTime()
          val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
          val price = cn.price(isCall, forward)
          var endTime = System.nanoTime()
          var elapsed = (endTime - startTime) * 1e-9
          val value = solver.impliedVolatility(isCall, forward, price, forward, df, tte) * 100
          val diff = value - previousValue
          val ratio = previousDiff / diff
          buffer += f" $spaceStep & $timeStep & $value%2.8f & $diff%2.1e & $ratio%2.1f & $elapsed%2.1e\\\\\n"
          previousValue = value
          previousDiff = diff

        }
        print(buffer)


      }
    }
  }

  test("pde-density-convergence-ah") {
    val isLatex = false
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "BDF2", "LMG2", "LMG3", "LS", "TRBDF2", "Bathe", "RE", "RE3")
    val spaceSteps = Array[Int](500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(5, 10, 20, 40, 80, 160, 320, 640)
    for (l <- 0 to 5) {
      println("TimeSteps & ATM Value & Max Difference & Ratio & Time(s)\\\\")
      for (scheme <- schemes) {
        for (spaceStep <- spaceSteps) {
          var buffer = ""
          println("\\hline \\multicolumn{5}{|c|}{" + scheme + "} \\\\ \\hline")
          var previousQ: Array[Double] = null
          var previousDiff = Double.NaN
          for (timeStep <- timeSteps) {
            val startTime = System.nanoTime()
            val cn = createHaganSABRSolver(scheme, spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
            cn.solve()
            val endTime = System.nanoTime()
            val elapsed = (endTime - startTime) * 1e-9
            val err = computeMaxError(forward, cn.Fm, cn.P, previousQ)
            val value: Double = err._1
            val diff: Double = err._2
            val ratio = previousDiff / diff
            buffer += f" $timeStep & $value%2.8f & $diff%2.1e & $ratio%2.1f & $elapsed%2.1e\\\\\n"
            previousQ = cn.P
            previousDiff = diff

          }
          print(buffer)

        }
      }
    }
  }

  test("pde-vol-convergence-ah") {
    val isLatex = false
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val FmaxTruncation = 4.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "BDF2", "LMG2", "LMG3", "LS", "TRBDF2", "Bathe", "RE")
    val spaceSteps = Array[Int](500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(2, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
    for (l <- 0 until 1) {
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

  test("HaganDensitySimple") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 80
    val FmaxTruncation = 4
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    var pdeMap = Map("CN" -> pde)
    //    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //    pde.useRannacher = true
    //    pdeMap += "RAN" -> pde
    //    pdeMap += "LMG2" -> new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //    pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //    pdeMap += "RE" -> new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)

    for ((name, pde) <- pdeMap) {
      pde.solve()
      println(name + " " + pde.computeCourantNumber + " " + pde.h)
    }
    println("densities")
    //density
    println("Scheme Strike Density")
    val strikeList = pde.Fm
    for (i <- 1 until strikeList.length - 1) {
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

      println(f"Hagan $strike%2.4f $densityHagan%2.4e")

      var normalVolUp = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike + eps, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      var normalVolDown = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike - eps, tte)
      val callUpNormal = BachelierVanillaEuropean.price(isCall, strike + eps, forward, normalVolUp, tte)
      val callNormal = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      val callDownNormal = BachelierVanillaEuropean.price(isCall, strike - eps, forward, normalVolDown, tte)
      val densityNormal = (callUpNormal - 2 * callNormal + callDownNormal) / (eps * eps)

      for ((name, pde) <- pdeMap) {
        val Q = pde.Q(i)

        println(f"$name $strike%2.4f $Q%2.4e")
      }
    }
  }

  def computeMaxError(solver: SABRDensitySolver, strikes: Array[Double], forward: Double, tte: Double, refSolver: SABRDensitySolver): Double = {
    val ivsolver = new Li2011SORDRBlackVolatilitySolver(1e-12)
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

  def computeMaxError(forward: Double, F: Array[Double], Q: Array[Double], Qold: Array[Double]): (Double, Double, Double) = {
    var maxError = 0.0
    var Fmax = 0.0
    var value = 0.0
    for (i <- F.size / 10 until F.size - F.size / 10) {
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


  def createHaganSABRSolver(name: String, spec: SABRModelSpec, forward: Double, tte: Double, spaceStep: Int = 1000, timeStep: Int = 100, nDeviation: Double = 3.0): HaganSABRTransformedDensitySolver = {
    var pde: HaganSABRTransformedDensitySolver = null
    name match {
      case "CN" => {
        val pdeH = new HaganSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
        pde = pdeH
      }
      case "RAN" => {
        val pdeH = new HaganSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
        pdeH.smoothing = new RannacherSmoothing()
        pde = pdeH
      }
      case "LEF" => {
        val pdeH = new HaganSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
        pdeH.smoothing = new BDF2Smoothing()
        pde = pdeH
      }
      case "LMG2" => pde = new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "LMG2a" => pde = new HaganLMG2aSABRTransformedDensitySolver(spec, forward, tte, spaceStep, 1e-4, nDeviation)
      case "LMG3" => pde = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "LS" => pde = new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "BDF2" => pde = new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "BDF2t" => {
        val pde2 = new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
        pde2.useTRinit = true
        pde = pde2
      }
      case "BDF3" => pde = new HaganBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "BDF3t" => {
        val pde2 = new HaganBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
        pde2.useTRinit = true
        pde = pde2
      }
      case "TRBDF2" => pde = new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "TRBDF2a" => pde = new HaganTRBDF2aSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      //case "CNBDF2" => pde = new HaganCNBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "Bathe" | "Bathe3" => pde = new HaganBathe3SubstepsSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "IMBDF3" => pde = new HaganIMBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "TRBDF3" | "CNBDF3" => pde = new HaganTRBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "RE" => pde = new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "RE3" => pde = new HaganRichardsonEuler3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "REN" => pde = new HaganRichardsonEulerNSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      //case "BDF2" => pde = new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      //case "BDF3" => pde = new HaganBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
    }
    return pde
  }

  test("AndreasenHugeDensityLocalVolTransform") {
    var alpha = 0.0873;
    //-0.01 translates AH formula to Hagan
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val spaceSize: Int = 50;
    val timeSize: Int = 5;
    val mu = 0.0;
    val r = 0.0;

    val pde = //new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, spaceSize, timeSize, 5.0*forward)
      new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, spaceSize, timeSize, 4.0)

    val fdsolver = new ThomasTridiagonalSolver();
    val payoffCall = new VanillaFDPayoff(false, forward, tte);
    val payoffPut = new VanillaFDPayoff(true, forward, tte);
    val timeMesh = new UniformMesh1D(timeSize, new Mesh1DBoundaries(0, tte))
    val spaceBoundaries = new Mesh1DBoundaries(0, 4 * forward);

    val methodCall = //new RannacherParabolic1DMethod()
      new LawsonSwayneParabolic1DMethod(payoffCall)
    //new TRBDF2Parabolic1DMethod(payoffCall)
    //new ThetaParabolic1DMethod()
    methodCall.smearingReducer = null
    val methodPut = //new RannacherParabolic1DMethod()
      new LawsonSwayneParabolic1DMethod(payoffPut)
    //new ThetaParabolic1DMethod()
    //new TRBDF2Parabolic1DMethod(payoffPut);
    methodPut.smearingReducer = null
    new TRBDF2Parabolic1DMethod(payoffCall);
    //         methodCall.lowerBoundary  = ForwardPartialOrder2Parabolic1DBoundaryFactory
    //         methodCall.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
    //    methodPut.lowerBoundary  = ForwardPartialOrder2Parabolic1DBoundaryFactory
    //    methodPut.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
    val grid =
      new StaticAdaptiveMesh2D(new UniformMesh1D(spaceSize - 1, spaceBoundaries, forward, true), timeMesh)

    val fdSpec = new SABRLocalVolFDSpec(grid, spec, forward, tte, mu, r)

    val pricerCall = new FDMSolver1D(fdSpec, methodCall, fdsolver);
    pricerCall.solve(payoffCall);
    val callLVSplineCall = CubicSpline.makeCubicSpline(grid.spaceVector, pricerCall.price)

    val pricerPut = new FDMSolver1D(fdSpec, methodPut, fdsolver);
    pricerPut.solve(payoffPut);
    val callLVSplinePut = CubicSpline.makeCubicSpline(grid.spaceVector, pricerPut.price)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val transformation = new SABRTransformation(spec, forward)
    val zb = transformation.inverseValue(0)
    val zmax = 4.0 * math.sqrt(tte)
    val zmin = math.max(-4.0 * math.sqrt(tte), zb)
    val zforward = transformation.inverseValue(forward)
    val gridT = new StaticAdaptiveMesh2D(new TransformedUniformMesh1D(transformation, spaceSize, new Mesh1DBoundaries(zmin, zmax), new Point(zforward, true), true), timeMesh)
    val fdSpecT = new SABRLocalVolFDSpec(gridT, spec, forward, tte, mu, r)
    val fdsolverT = new ThomasTridiagonalSolver();
    val pricerCallT = new FDMSolver1D(fdSpecT, methodCall, fdsolverT);
    // pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricerCallT.solve(payoffCall);
    val callLVSplineCallT = Array[CubicPP](
      CubicSpline.makeHarmonicSpline(gridT.spaceVector, pricerCallT.price),
      CubicSpline.makeBesselSpline(gridT.spaceVector, pricerCallT.price),
      CubicSpline.makeCubicSpline(gridT.spaceVector, pricerCallT.price),
      CubicSpline.makeHymanCubicSpline(gridT.spaceVector, pricerCallT.price))

    val pricerPutT = new FDMSolver1D(fdSpecT, methodPut, fdsolver);
    // pricer.smoother = new SimpsonIntegralSmoother(Array(strike))
    pricerPutT.solve(payoffPut);
    val callLVSplinePutT = Array[CubicPP](CubicSpline.makeHarmonicSpline(gridT.spaceVector, pricerPutT.price),
      CubicSpline.makeBesselSpline(gridT.spaceVector, pricerPutT.price),
      CubicSpline.makeCubicSpline(gridT.spaceVector, pricerPutT.price),
      CubicSpline.makeHymanCubicSpline(gridT.spaceVector, pricerPutT.price))

    println("black vols")
    println("Strike Hagan NormalHagan HaganTransormedDensity HaganLocalVol HaganTransformedLocalVol")
    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 10.0)
      val isCall = strike >= forward
      val haganVol = SABRVanilla.impliedVolatilityHagan(spec, forward, strike, tte)
      var normalVol = SABRVanilla.normalVolatilityHagan2013(spec, forward, strike, tte)
      var normalPrice = BachelierVanillaEuropean.price(isCall, strike, forward, normalVol, tte)
      normalVol = solver.impliedVolatility(isCall, strike, normalPrice, forward, df, tte)

      val splineLVT = if (isCall) callLVSplineCallT else callLVSplinePutT
      val lvtVol = Array.ofDim[Double](splineLVT.length)
      var j = 0
      for (spline <- splineLVT) {
        val lvtPrice = spline.value(strike)
        lvtVol(j) = Double.NaN
        try {
          lvtVol(j) = solver.impliedVolatility(isCall, strike, lvtPrice, forward, df, tte)
        } catch {
          case e: RuntimeException =>
        }
        j += 1
      }


      val splineLV = if (isCall) callLVSplineCall else callLVSplinePut
      val lvPrice = splineLV.value(strike)
      var lvVol = Double.NaN
      try {
        lvVol = solver.impliedVolatility(isCall, strike, lvPrice, forward, df, tte)
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
      val lvtVolHarm = lvtVol(0)
      val lvtVolBess = lvtVol(1)
      val lvtVolCubi = lvtVol(2)
      val lvtVolHyma = lvtVol(3)
      println(f"$strike%2.4f $haganVol%2.3f $normalVol%2.3f $pdeVol%2.4f $lvVol%2.4f $lvtVolHarm%2.4f $lvtVolBess%2.4f $lvtVolCubi%2.4f $lvtVolHyma%2.4f")

    }

    println("densities")
    //density
    // println("Strike Hagan NormalHagan NormalHaganScaled AndreasenHugeRaw AndreasenHuge Benhamou")
    println("Strike Hagan NormalHagan HaganDensity HaganUniformLV HaganLVHarmonic HaganLVBessel HaganLVCubic HaganLVHyman")

    for (i <- 1 to 100) {
      val strike = forward * (0.0 + i / 20.0)
      val eps = 1e-4 * strike
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

      val callPDE = pde.price(isCall, strike)
      val callUpPDE = pde.price(isCall, strike + eps)
      val callDownPDE = pde.price(isCall, strike - eps)
      val densityPDE = (callUpPDE - 2 * callPDE + callDownPDE) / (eps * eps)

      val splineLV = if (isCall) callLVSplineCall else callLVSplinePut
      val densityLV = splineLV.secondDerivative(strike);

      val splineLVT = if (isCall) callLVSplineCallT else callLVSplinePutT
      val densityLVT = Array.ofDim[Double](splineLVT.length)
      var j = 0
      for (spline <- splineLVT) {
        densityLVT(j) = spline.secondDerivative(strike)
        j += 1
      }
      val densityLVTHarm = densityLVT(0)
      val densityLVTBess = densityLVT(1)
      val densityLVTCubi = densityLVT(2)
      val densityLVTHyma = densityLVT(3)

      println(f"$strike%2.4f $densityHagan%2.2e $densityNormal%2.2e $densityPDE%2.3e $densityLV%2.3e $densityLVTHarm%2.3e $densityLVTBess%2.3e $densityLVTCubi%2.3e $densityLVTHyma%2.3e")
      //       println(f"$strike%2.4f $call%2.2e $callNormal%2.2e $callBen%2.2e")
    }

  }
  test("AndreasenHugeDensityTransform") {
    val alpha = 0.0873;
    val nu = 0.47;
    val beta = 0.7;
    val rho = -0.48;
    val forward = 0.025;
    val tte = 10.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val spaceSize: Int = 50;
    val timeSize: Int = 10;
    val mu = 0.0;
    val r = 0.0;

    val pde = new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, spaceSize, timeSize, 5.0)
    val pdeU = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, spaceSize, timeSize, 5.0 * forward)

    println("refined call price around forward")
    println("strike price priceUniform density")
    // forward*0.1
    val strikeRef = forward
    val isCall = false
    for (i <- 0 to 100) {
      val strike = strikeRef + (-strikeRef * 0.1 + i * (strikeRef * 0.1 * 2) / 100)
      val callPrice = pde.price(isCall, strike)
      val callPriceU = pdeU.price(isCall, strike)
      var dens = 0.0
      if (i > 0 && i < 100) {
        val up = pde.price(isCall, strike * 1.01)
        val down = pde.price(isCall, strike * 0.99)
        dens = (-2 * callPrice + up + down) / (strike * 0.01) / (strike * 0.01)
      }
      var densU = 0.0
      if (i > 0 && i < 100) {
        val up = pdeU.price(isCall, strike * 1.01)
        val down = pdeU.price(isCall, strike * 0.99)
        densU = (-2 * callPriceU + up + down) / (strike * 0.01) / (strike * 0.01)
      }
      println(f"$strike%2.6f $callPrice%2.6f $callPriceU%2.6f $dens%f $densU%f")
    }
  }

  test("HaganDensityTransform") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    val spaceSize: Int = 50;
    val timeSize: Int = 50;
    val mu = 0.0;
    val r = 0.0;

    val pde = new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, spaceSize, timeSize, 5.0)
    val pdeU = new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, spaceSize, timeSize, 5.0 * forward)

    println("refined call price around forward")
    println("strike price priceUniform density")
    // forward*0.1
    val strikeRef = forward
    val isCall = false
    for (i <- 0 to 100) {
      val strike = strikeRef + (-strikeRef * 0.1 + i * (strikeRef * 0.1 * 2) / 100)
      val callPrice = pde.price(isCall, strike)
      val callPriceU = pdeU.price(isCall, strike)
      var dens = 0.0
      if (i > 0 && i < 100) {
        val up = pde.price(isCall, strike * 1.01)
        val down = pde.price(isCall, strike * 0.99)
        dens = (-2 * callPrice + up + down) / (strike * 0.01) / (strike * 0.01)
      }
      var densU = 0.0
      if (i > 0 && i < 100) {
        val up = pdeU.price(isCall, strike * 1.01)
        val down = pdeU.price(isCall, strike * 0.99)
        densU = (-2 * callPriceU + up + down) / (strike * 0.01) / (strike * 0.01)
      }
      println(f"$strike%2.6f $callPrice%2.6f $callPriceU%2.6f $dens%f $densU%f")
    }
  }
}


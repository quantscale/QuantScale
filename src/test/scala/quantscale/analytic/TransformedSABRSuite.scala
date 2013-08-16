package quantscale.analytic

import quantscale.fdm.sabr._
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.scalatest.{Tag, FunSuite}
import quantscale.fdm.mesh.{Mesh1DBoundaries, UniformMesh1D}
import quantscale.vol.Li2011SORDRBlackVolatilitySolver
import quantscale.fdm.Epsilon
import quantscale.math.{CubicPP, CubicSpline}
import org.junit.Test

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

    val normalVol = 100*SABRVanilla.normalVolatilityHagan2013(spec, forward, K, tte)
    val normalPrice = BachelierVanillaEuropean.price(false, K, forward, normalVol/100, tte)
    val normalPriceParity = K-forward +  BachelierVanillaEuropean.price(true, K, forward, normalVol/100, tte)
    //    val volNormal = 100*solver.impliedVolatility(true, K, normalPrice, forward, df, tte)
    val volAH = 100 * solver.impliedVolatility(false, K, priceAH, forward, df, tte)
    val volAHRaw = 100 * solver.impliedVolatility(false, K, priceAHRaw, forward, df, tte)

    val lognormalPrice = BlackScholesVanillaEuropean.priceEuropeanVanilla(false, K, forward, volHagan/100 * volHagan/100 * tte, 1.0, df)
    println(f"HaganLN $lognormalPrice%2.3f $volHagan%2.3f")
    println(f"HaganN $normalPrice%2.3f $normalVol%2.3f")
    println(normalPriceParity)
    println(f"priceAH $priceAH%2.3f $volAH%2.3f")

//    val Fmax = Array(5, 50, 500, 5000, 10000)
//    val ndev = Array(3, 4, 5, 10, 100)
//    val size = Array(10, 20, 100, 1000, 10000)
//    val timesteps = Array(10, 40, 160)

    val Fmax = Array(50)
    val ndev = Array(10)
    val size = Array(10)
    val timesteps = Array(4)

    for (tsteps <- timesteps) {
       for (xsteps <- size) {
         for (i <- 0 until Fmax.length) {
           val pde = new HaganLMG3SABRDensitySolver(spec, forward, tte, xsteps, tsteps, Fmax(i))
           val pdeTransformed = new HaganLMG3SABRTransformedDensitySolver(spec,forward,tte, xsteps ,tsteps,ndev(i))
           val pricePDE = pde.price(false, K)
           val pricePDET = pdeTransformed.price(false, K)

           val volPDE = 100 * solver.impliedVolatility(false, K, pricePDE, forward, df, tte)
           val volPDET = 100 * solver.impliedVolatility(false, K, pricePDET, forward, df, tte)

           print(Fmax(i) +" & "+xsteps+" & "+tsteps+ f" & $pricePDE%2.5f & $volPDE%2.3f && ")
           println(ndev(i) +" & "+xsteps+" & "+tsteps+ f" & $pricePDET%2.5f & $volPDET%2.3f")
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
      println(zm(j) + " "+Fm(j))
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
    pde.useRannacher = false
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.useRannacher = true
    pdeMap += "RAN" -> pde


    pdeMap += "LMG2" -> new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF3" -> new HaganCNBDF3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)

    for ((name, pde) <- pdeMap) {
      pde.solve()
      val pricePut = pde.price(false,forward)
      val priceCall = pde.price(true, forward)
      println(name+" "+pricePut+" "+priceCall)
      assert(pricePut-priceCall < 1e-8, "put="+pricePut+" call="+priceCall)

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
    pde.useRannacher = false
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.useRannacher = true
    pdeMap += "RAN" -> pde


    pdeMap += "LMG2" -> new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF3" -> new HaganCNBDF3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    val h = pde.h
    val dt = pde.dt
    println(f"h=$h%2.12f")
    println(f"dt=$dt%2.12f")
    println("courant="+pde.computeCourantNumber())
    println("Scheme & ATM price & Q(f) & QL & QR\\\\")
    for (i <- 0 until 1) {
      val startTime = System.nanoTime()
      for ((name, pde) <- pdeMap) {
        pde.solve()
        val j0 = pde.indexForward
        val Qforward = pde.P(j0)
        val QL = pde.PL
        val QR = pde.PR

        val price = pde.price(false,forward)
        println(f"$name & $price%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")

        // println(f"$name & $h%2.12f & $Qforward%2.12f & $QL%2.12f & $QR%2.12f\\\\")
      }
      println("elapsed "+(System.nanoTime()-startTime)*1e-9)
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
    val refPDE = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, 1280*4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    val Qspline = CubicSpline.makeCubicSpline(refPDE.Fm, refPDE.P)

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "TRBDF3", "RE")
    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1,2,3,4, 5,6,7,8,9, 10, 20, 40, 80, 160, 320, 640, 1280)
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
    val refPDE = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, 1280*4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    val Qspline = CubicSpline.makeCubicSpline(refPDE.Fm, refPDE.P)

    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    val schemes = Array("CN", "RAN", "LMG2", "LMG3", "LS", "TRBDF2", "TRBDF3", "RE")
    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1,2,3,4, 5,6,7,8,9, 10, 20, 40, 80, 160, 320, 640, 1280)
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
    val FmaxTruncation = 5.0
    val strikes = new UniformMesh1D(100, new Mesh1DBoundaries(forward * 0.2, forward * 2)).x

    var isCall = false
    val refPDE = new HaganLMG3SABRDensitySolver(spec, forward, tte, 1280 * 4, 1280 * 4, FmaxTruncation)

    refPDE.solve()
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)

    val spaceSteps = Array(500) // Array(80, 160, 320, 640, 1280)
    val timeSteps = Array(1, 2, 3, 4, 5, 6,7,8,9, 10, 20, 40, 80, 160, 320, 640, 1280)
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
    val timeSteps = Array(1, 2, 3, 4, 5, 6,7,8,9, 10, 20, 40, 80, 160, 320, 640, 1280)
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
    val FmaxTruncation = 3.0
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
    val FmaxTruncation = 3.0
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

  test("HaganDensitySimple") {
    val alpha = 0.35;
    val nu = 1.0;
    val beta = 0.25;
    val rho = -0.1;
    val forward = 1.0;
    val tte = 1.0;
    val df = 1.0
    val xSteps = 500
    val tSteps = 5
    val FmaxTruncation = 3
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var solver = new Li2011SORDRBlackVolatilitySolver(1e-12)
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useRannacher = false
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pde.useRannacher = true
    pdeMap += "RAN" -> pde
    pdeMap += "LMG2" -> new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "TRBDF2" -> new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "RE" -> new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)

    for ((name, pde) <- pdeMap) {
      pde.solve()
      println(name + " " + pde.computeCourantNumber + " " + pde.h)
    }
    println("densities")
    //density
    println("Scheme Strike Density")
    val strikeList = pde.Fm
    for (i <- 1 until strikeList.length-1) {
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

        println(f"$name $strike%2.4f $Q%2.2e")
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
    for (i <- F.size/10 until F.size - F.size/10) {
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
        pdeH.useRannacher = false
        pde = pdeH
      }
      case "RAN" => {
        val pdeH = new HaganSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
        pdeH.useRannacher = true
        pde = pdeH
      }
      case "LMG2" => pde = new HaganLMG2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "LMG3" => pde = new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "LS" => pde = new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "TRBDF2" => pde = new HaganTRBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      //case "CNBDF2" => pde = new HaganCNBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      case "TRBDF3" | "CNBDF3" => pde = new HaganCNBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      case "RE" => pde = new HaganRichardsonEulerSABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, nDeviation)
      //case "BDF2" => pde = new HaganBDF2SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
      //case "BDF3" => pde = new HaganBDF3SABRTransformedDensitySolver(spec, forward, tte, spaceStep, timeStep, FmaxTruncation)
    }
    return pde
  }
}


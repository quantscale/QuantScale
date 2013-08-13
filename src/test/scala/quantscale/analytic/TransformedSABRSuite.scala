package quantscale.analytic

import quantscale.fdm.sabr._
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.scalatest.FunSuite
import quantscale.fdm.mesh.{Mesh1DBoundaries, UniformMesh1D}
import quantscale.vol.Li2011SORDRBlackVolatilitySolver
import quantscale.fdm.Epsilon
import quantscale.math.{CubicPP, CubicSpline}

@RunWith(classOf[JUnitRunner])
class TransformedSABRSuite extends FunSuite {

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
        val Qforward = pde.Q0(j0)
        val QL = pde.QL
        val QR = pde.QR

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
    val Qspline = CubicSpline.makeCubicSpline(refPDE.Fm, refPDE.Q0)

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
    val Q0 = solver.Q0

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
            val err = computeMaxError(forward, cn.Fm, cn.Q0, previousQ)
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

  def computeMaxError(solver: SABRDensitySolver, strikes: Array[Double], forward: Double, tte: Double, refSolver: HaganSABRDensitySolver): Double = {
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


package quantscale.analytic

import quantscale.fdm.sabr.{HaganLMG3SABRTransformedDensitySolver, HaganLawsonSwayneSABRTransformedDensitySolver, HaganSABRTransformedDensitySolver, HaganSABRDensitySolver}
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.scalatest.FunSuite

@RunWith(classOf[JUnitRunner])
class TransformedSABRSuite extends FunSuite {
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
    val nDeviation = 3.0
    val spec = new SABRModelSpec(alpha, beta, nu, rho)
    var pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.useRannacher = false
    //var pdeMap = Map("LS" -> new HaganLawsonSwayneSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation))
    var pdeMap = Map("CN" -> pde)
    pde = new HaganSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pde.useRannacher = true
    pdeMap += "RAN" -> pde


    //pdeMap += "LMG2" -> new HaganLMG2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    pdeMap += "LMG3" -> new HaganLMG3SABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    pdeMap += "LS" -> new HaganLawsonSwayneSABRTransformedDensitySolver(spec, forward, tte, xSteps, tSteps, nDeviation)
    //pdeMap += "TRBDF2" -> new HaganTRBDF2SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "RE" -> new HaganRichardsonEulerSABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    // pdeMap += "ADE" -> new HaganADESABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
    //pdeMap += "TRBDF3" -> new HaganCNBDF3SABRDensitySolver(spec, forward, tte, xSteps, tSteps, FmaxTruncation)
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
}


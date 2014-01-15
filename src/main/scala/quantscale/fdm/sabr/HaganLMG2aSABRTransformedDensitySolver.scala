package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.{Epsilon, TridiagonalMatrix}

class HaganLMG2aSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, tolerance: Double, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, 1, nDeviations) {

  private var timeStepCounter_ = 0
  private val Q1Half = Array.ofDim[Double](size)

  private[sabr] def buildEm(t: Double) {
    var j = size - 2
    val efactor = spec.rho * spec.nu * spec.alpha
    while (j > 0) {
      Em_(j) = math.exp(efactor * Gammam_(j) * t)
      j -= 1
    }
    Em_(0) = Em_(1)
    Em_(size - 1) = Em_(size - 2)
  }

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    Em_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0

    var ddt = T / 1000
    //evaluate first step size
    //val error = math.abs(solveStep(T, dt))
    //var ddt = math.max(T / 10000, math.sqrt(tolerance / error))
    //println("ddt init="+ddt)
    P0_ = computeP()
    PL_ = 0.0
    PR_ = 0.0
    val dtmin = T / 10000
    var tIndex = 0
    //var ddt = dt
    var t = T
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //             if (tIndex < 4) {
      //              printTimeStep(t)
      //            }

      ddt = math.min(ddt, t)
      val PLold = PL_
      val PRold = PR_
      val error = math.abs(solveStep(t, ddt))

      tIndex += 1
      val r = tolerance * math.abs(P0_(j0)) / error
      var newddt = math.pow(r, 1.0 / 3.0) * ddt
      if (r > 0.5 || ddt == dtmin) {
        t = t - ddt
      } else {
        newddt = 0.9 * newddt
        val Qtmp = P0_
        P0_ = P1_
        P1_ = Qtmp
        PL_ = PLold
        PR_ = PRold
      }
      ddt = math.max(dtmin, newddt) //error*dt^2 = tolerance
      println("t=" + t + " ddt=" + ddt + " error=" + error + " r=" + r)

    }
    println("tCounter=" + tIndex)
  }

  private def solveStep(t0: Double, dt0: Double): Double = {
    var t = t0 - dt0 / 2
    //      System.arraycopy(M0, 0, M1, 0, size)
    buildEm(T - t)
    computeSystem(dt0 / 2, Em_, tri1)
    P0_(0) = 0
    P0_(size - 1) = 0
    solver.solve(tri1, P0_, Q1Half)
    var QL_Part = dt0 / 2 * computedPLdt(Em_, Q1Half)
    var QR_Part = dt0 / 2 * computedPRdt(Em_, Q1Half)
    val Q0_Init = P0_
    P0_ = Q1Half

    t -= dt0 / 2
    buildEm(T - t)
    computeSystem(dt0 / 2, Em_, tri1)
    P0_(0) = 0
    P0_(size - 1) = 0
    solver.solve(tri1, P0_, Q1Half)
    QL_Part += dt0 / 2 * computedPLdt(Em_, Q1Half)
    QR_Part += dt0 / 2 * computedPRdt(Em_, Q1Half)

    P0_ = Q0_Init
    computeSystem(dt0, Em_, tri1)
    P0_(0) = 0
    P0_(size - 1) = 0
    solver.solve(tri1, P0_, P1_)

    PL_ += 2 * QL_Part - dt0 * computedPLdt(Em_, P1_)
    PR_ += 2 * QR_Part - dt0 * computedPRdt(Em_, P1_)


    val error = 2 * (P1_(j0) - Q1Half(j0))
    var i = 0
    while (i < size) {
      P1_(i) = 2 * Q1Half(i) - P1_(i)
      i += 1
    }
    //printSumQF("LMG2",t,Q1,QL_,QR_)
    //check sum
    //      var sum = QR_+QL
    //      i = size-2
    //      while (i>0) {
    //        sum += h*Q1(i)
    //        i-=1
    //      }
    //      println("t="+t+" LMG2 Q="+sum)
    //      M0 = M1


    val Qtmp = P0_
    P0_ = P1_
    P1_ = Qtmp
    return error * error
  }


}
package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganTRBDF2SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = 0.5 * alpha;
  private val frac = 1.0 / (alpha * (2 - alpha));

  override def solve() {
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val M0 = Array.ofDim[Double](size)
    val bdt = 2 * dt * backcoeff
    buildEmCache(bdt, dt - bdt)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    val rhs = Array.ofDim[Double](size)
    val cdth = backcoeff * dt / (h)
    var tIndex = 0

    while (tIndex < timeSteps) {
//                        if (tIndex < 4) {
//                    printTimeStep(t)
//                  }

      t -= dt
      Array.copy(Em_, 0, M0, 0, size)
      advanceEm(bdt, Em_)
      computeSystem(bdt, Em_, M0, tri1, tri0)
      tri0.multiply(P0_, rhs)
      solver.solve(tri1, rhs, P1_)
      val QL_init = PL_
      val QR_init = PR_
      PL_ += 0.5 * bdt * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * P1_(1) + M0(1) * P0_(1))
      PR_ += 0.5 * bdt * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * P1_(size - 2) + M0(size - 2) * P0_(size - 2))
      //      printSumQF("TR", t, Q1, QL_, QR_)

      var j = size - 2
      while (j >= 1) {
        rhs(j) = frac * (P1_(j) - (1 - alpha) * (1 - alpha) * P0_(j))
        j -= 1
      }
      advanceEm(dt - bdt, Em_)
      computeSystem(bdt, Em_, null, tri1, tri0)
      solver.solve(tri1, rhs, P1_)

      PL_ = frac * (PL - (1 - alpha) * (1 - alpha) * QL_init) + 0.5 * bdt * computedPLdt(Em_, P1_)
      PR_ = frac * (PR - (1 - alpha) * (1 - alpha) * QR_init) + 0.5 * bdt * computedPRdt(Em_, P1_)

      //      printSumQF("BDF2 ", t, Q1, QL_, QR_)

      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp

      tIndex += 1
    }
  }

}
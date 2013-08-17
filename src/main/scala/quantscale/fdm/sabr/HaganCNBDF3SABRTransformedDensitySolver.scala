package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.{TridiagonalMatrix}

class HaganCNBDF3SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {


  override def solve() {
    val dt3 = dt / 3.0
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val M0 = Array.ofDim[Double](size)
    val Q2 = Array.ofDim[Double](size)
    buildEmCache(dt3, 0)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0

    var t = T
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (tIndex < timeSteps) {
//                        if (tIndex < 4) {
//                    printTimeStep(t)
//                  }

      t -= dt
      Array.copy(Em_, 0, M0, 0, size)
      advanceEm(dt3, Em_)
      computeSystem(dt3, Em_, M0, tri1, tri0)
      tri0.multiply(P0_, rhs)
      solver.solve(tri1, rhs, P1_)
      val QL1 = PL_
      val QR1 = PR_
      PL_ += 0.5 * dt3 * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * P1_(1) + M0(1) * P0_(1))
      PR_ += 0.5 * dt3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * P1_(size - 2) + M0(size - 2) * P0_(size - 2))

      Array.copy(Em_, 0, M0, 0, size)
      advanceEm(dt3, Em_)
      computeSystem(dt3, Em_, M0, tri1, tri0)
      tri0.multiply(P1_, rhs)
      solver.solve(tri1, rhs, Q2)
      val QL2 = PL_
      val QR2 = PR_
      PL_ += 0.5 * dt3 * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * Q2(1) + M0(1) * P1_(1))
      PR_ += 0.5 * dt3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * Q2(size - 2) + M0(size - 2) * P1_(size - 2))

      var j = size - 2
      while (j >= 1) {
        rhs(j) = (18 * Q2(j) - 9 * P1_(j) + 2 * P0_(j)) / 11
        j -= 1
      }
      advanceEm(dt3, Em_)
      computeSystem(dt3 * 6.0 / 11, Em_, tri1)
      solver.solve(tri1, rhs, P1_)

      PL_ = (18 * PL_ - 9 * QL2 + 2 * QL1) / 11 + 6 * dt3 / 11 * computedPLdt(Em_, P1_)
      PR_ = (18 * PR_ - 9 * QR2 + 2 * QR1) / 11 + 6 * dt3 / 11 * computedPRdt(Em_, P1_)

      //printSumQF("TRBDF3", t, Q1, QL_, QR_)
      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }
  }

}

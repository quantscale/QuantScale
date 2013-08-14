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
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0

    var t = T
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (tIndex < timeSteps) {
      //                  if (tIndex < 4) {
      //              printTimeStep(t)
      //              tIndex += 1
      //            }

      t -= dt
      Array.copy(Em_, 0, M0, 0, size)
      advanceEm(dt3, Em_)
      computeSystem(dt3, Em_, M0, tri1, tri0)
      tri0.multiply(Q0_, rhs)
      solver.solve(tri1, rhs, Q1)
      val QL1 = QL_
      val QR1 = QR_
      QL_ += 0.5 * dt3 * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * Q1(1) + M0(1) * Q0_(1))
      QR_ += 0.5 * dt3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * Q1(size - 2) + M0(size - 2) * Q0_(size - 2))

      Array.copy(Em_, 0, M0, 0, size)
      advanceEm(dt3, Em_)
      computeSystem(dt3, Em_, M0, tri1, tri0)
      tri0.multiply(Q1, rhs)
      solver.solve(tri1, rhs, Q2)
      val QL2 = QL_
      val QR2 = QR_
      QL_ += 0.5 * dt3 * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * Q2(1) + M0(1) * Q1(1))
      QR_ += 0.5 * dt3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * Q2(size - 2) + M0(size - 2) * Q1(size - 2))

      var j = size - 2
      while (j >= 1) {
        rhs(j) = (18 * Q2(j) - 9 * Q1(j) + 2 * Q0_(j)) / 11
        j -= 1
      }
      advanceEm(dt3, Em_)
      computeSystem(dt3 * 6.0 / 11, Em_, tri1)
      solver.solve(tri1, rhs, Q1)

      QL_ = (18 * QL_ - 9 * QL2 + 2 * QL1) / 11 + 6 * dt3 / 11 * computedQLdt(Em_, Q1)
      QR_ = (18 * QR_ - 9 * QR2 + 2 * QR1) / 11 + 6 * dt3 / 11 * computedQRdt(Em_, Q1)

      //printSumQF("TRBDF3", t, Q1, QL_, QR_)
      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
      tIndex += 1
    }
  }

}

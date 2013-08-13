package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganRichardsonEulerSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  override def solve() {
    solve(1)
    val Q0_full = Array.ofDim[Double](size)
    Array.copy(Q0_, 0, Q0_full, 0, size) // a new array will be create later
    val QL_full = QL_
    val QR_full = QR_
    dt_ = dt_ / 2
    solve(2)
    var i = 0
    while (i < size) {
      Q0_(i) = 2 * Q0_(i) - Q0_full(i)
      i += 1
    }
    dt_ = dt_ * 2
    QL_ = 2 * QL_ - QL_full
    QR_ = 2 * QR_ - QR_full
  }

  def solve(divisor : Int) {
    tri1 = new TridiagonalMatrix(size)
    buildEmCache(dt_, 0)

    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    var tIndex = 0
    while (tIndex < timeSteps*divisor) {
      t -= dt
      advanceEm(dt_, Em_)
      computeSystem(dt_, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1)
      QL_ += dt_ * computedQLdt(Em_, Q1)
      QR_ += dt_ * computedQRdt(Em_, Q1)
   //   printSumQF("RE",t,Q1, QL_, QR_)
      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
      tIndex += 1
    }
  }


}
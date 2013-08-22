package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganRichardsonEulerSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  override def solve() {
    buildEmCache(dt_, dt_ / 2)
    val EmInit = Em_.clone()
    solve(1)
    val Q0_full = Array.ofDim[Double](size)
    Array.copy(P0_, 0, Q0_full, 0, size) // a new array will be create later
    val QL_full = PL_
    val QR_full = PR_
    dt_ = dt_ / 2
    Array.copy(EmInit, 0, Em_, 0, size)
    solve(2)
    var i = 0
    while (i < size) {
      P0_(i) = 2 * P0_(i) - Q0_full(i)
      i += 1
    }
    dt_ = dt_ * 2
    PL_ = 2 * PL_ - QL_full
    PR_ = 2 * PR_ - QR_full
  }

  def solve(divisor: Int) {
    tri1 = new TridiagonalMatrix(size)


    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    var tIndex = 0
    while (tIndex < timeSteps * divisor) {
      t -= dt
      advanceEm(dt_, Em_)
      computeSystem(dt_, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, P1_)
      PL_ += dt_ * computedPLdt(Em_, P1_)
      PR_ += dt_ * computedPRdt(Em_, P1_)
      //   printSumQF("RE",t,Q1, QL_, QR_)
      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }
  }


}
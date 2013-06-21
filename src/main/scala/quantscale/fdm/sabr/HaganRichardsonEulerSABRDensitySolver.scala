package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganRichardsonEulerSABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  override def solve() {
    val Q1full = Array.ofDim[Double](size)
    solve(dt_)
    System.arraycopy(Q1, 0, Q1full, 0, size)
    val QL_full = QL_
    val QR_full = QR_
    dt_ = dt_ / 2
    solve(dt_)
    var i = 0
    while (i < size) {
      Q1(i) = 2 * Q1(i) - Q1full(i)
      i += 1
    }
    dt_ = dt_ * 2
    val Qtmp = Q0_
    Q0_ = Q1
    Q1 = Qtmp
    QL_ = 2 * QL_ - QL_full
    QR_ = 2 * QR_ - QR_full
  }

  def solve(dt: Double) {
    tri1 = new TridiagonalMatrix(size)
    //    tri0 = new TridiagonalMatrix(size)
    M1 = Array.ofDim[Double](size)
    buildMcache(dt, 0, M1)

    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    //    var rhs = Array.ofDim[Double](size)
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      t -= dt
      advance(dt, M1)
      computeSystem(dt)
      //      System.arraycopy(Q0_, 0, rhs,0, size)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1)
      QL_ += dt / (h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
      QR_ -= dt / (h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))
//      M0 = M1
      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
    }
  }

  private[sabr] override def computeSystem(dt: Double) {

    var j = 1
    val frac = dt / (h * h)
    while (j < size - 1) {
      tri1.lower(j) = -frac * M1(j - 1)
      tri1.middle(j) = 1 + 2 * frac * M1(j)
      tri1.upper(j) = -frac * M1(j + 1)
      //      tri0.lower(j) = frac * M0(j - 1)
      //      tri0.middle(j) = 1 - 2 * frac * M0(j)
      //      tri0.upper(j) = frac * M0(j + 1)
      j += 1
    }
    tri1.upper(0) = M1(1)
    tri1.middle(0) = M1(0)
    tri1.lower(size - 1) = M1(size - 2)
    tri1.middle(size - 1) = M1(size - 1)
  }
}
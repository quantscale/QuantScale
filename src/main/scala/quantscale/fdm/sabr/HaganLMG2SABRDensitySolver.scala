package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganLMG2SABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    //    val tri1Half = new TridiagonalMatrix(size) //Half
    //    M0 = Array.ofDim[Double](size)
    M1 = Array.ofDim[Double](size)
    buildMcache(dt / 2, 0, M1)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    val Q1Half = Array.ofDim[Double](size)
    QL_ = 0.0
    var QL_Full = 0.0
    var QR_Full = 0.0
    QR_ = 0.0
    var t = T
    //    var rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //       if (tIndex < 4) {
      //        printTimeStep(t)
      //        tIndex += 1
      //      }
      t -= dt / 2
      //      System.arraycopy(M0, 0, M1, 0, size)
      advance(dt / 2, M1)
      computeSystem(tri1, dt / 2)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      var QL_Part = dt / (2 * h) * (M1(1) * Q1Half(1) - M1(0) * Q1Half(0))
      var QR_Part = dt / (2 * h) * (M1(size - 1) * Q1Half(size - 1) - M1(size - 2) * Q1Half(size - 2))
      val Q0_Init = Q0_
      Q0_ = Q1Half

      t -= dt / 2
      advance(dt / 2, M1)
      computeSystem(tri1, dt / 2)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      QL_Part += dt / (2 * h) * (M1(1) * Q1Half(1) - M1(0) * Q1Half(0))
      QR_Part += dt / (2 * h) * (M1(size - 1) * Q1Half(size - 1) - M1(size - 2) * Q1Half(size - 2))

      Q0_ = Q0_Init
      computeSystem(tri1, dt)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1)

      QL_ += 2 * QL_Part - dt / (h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
      QR_ -= 2 * QR_Part - dt / (h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

      var i = 0
      while (i < size) {
        Q1(i) = 2 * Q1Half(i) - Q1(i)
        i += 1
      }

      //check sum
      //      var sum = QR_+QL
      //      i = size-2
      //      while (i>0) {
      //        sum += h*Q1(i)
      //        i-=1
      //      }
      //      println("t="+t+" LMG2 Q="+sum)
      //      M0 = M1
      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
    }

  }

  private def computeSystem(tri1: TridiagonalMatrix, dt: Double) {

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
package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganLMG3SABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    //    M0 = Array.ofDim[Double](size)
    M1 = Array.ofDim[Double](size)
    buildMcache(dt / 3, 0, M1)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    var Q1Third = Array.ofDim[Double](size)
    var Q1ThirdTmp = Array.ofDim[Double](size)
    var Q1Half = Array.ofDim[Double](size)
    QL_ = 0.0
    var QL_Full = 0.0
    var QR_Full = 0.0
    QR_ = 0.0
    var t = T
    var Qtmp = Q0_
    var tIndex = 0

    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //                  if (tIndex < 4) {
      //              printTimeStep(t)
      //              tIndex += 1
      //            }
      t -= dt / 3

      advance(dt / 3, M1)
      computeSystem(tri1, dt / 3)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1ThirdTmp)
      var QL_ThirdPart = QL_ + dt / (3 * h) * (M1(1) * Q1ThirdTmp(1) - M1(0) * Q1ThirdTmp(0))
      var QR_ThirdPart = QR_ - dt / (3 * h) * (M1(size - 1) * Q1ThirdTmp(size - 1) - M1(size - 2) * Q1ThirdTmp(size - 2))
      var QL_HalfPart = QL_ThirdPart
      var QR_HalfPart = QR_ThirdPart
      val Q0_Init = Q0_
      Q0_ = Q1ThirdTmp

      t -= dt / 3
      advance(dt / 3, M1)
      computeSystem(tri1, dt / 3)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Third)
      QL_ThirdPart += dt / (3 * h) * (M1(1) * Q1Third(1) - M1(0) * Q1Third(0))
      QR_ThirdPart -= dt / (3 * h) * (M1(size - 1) * Q1Third(size - 1) - M1(size - 2) * Q1Third(size - 2))
      Q0_ = Q1Third

      t -= dt / 3
      advance(dt / 3, M1)
      computeSystem(tri1, dt / 3)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Third)
      QL_ThirdPart += dt / (3 * h) * (M1(1) * Q1Third(1) - M1(0) * Q1Third(0))
      QR_ThirdPart -= dt / (3 * h) * (M1(size - 1) * Q1Third(size - 1) - M1(size - 2) * Q1Third(size - 2))

      Q0_ = Q1ThirdTmp
      computeSystem(tri1, 2 * dt / 3)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      QL_HalfPart += 2 * dt / (3 * h) * (M1(1) * Q1Half(1) - M1(0) * Q1Half(0))
      QR_HalfPart -= 2 * dt / (3 * h) * (M1(size - 1) * Q1Half(size - 1) - M1(size - 2) * Q1Half(size - 2))

      Q0_ = Q0_Init
      computeSystem(tri1, dt)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1)

      QL_ += 4.5 * QL_ThirdPart - 4.5 * QL_HalfPart + dt / (h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
      QR_ += 4.5 * QR_ThirdPart - 4.5 * QR_HalfPart - dt / (h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

      var i = 0
      while (i < size) {
        Q1(i) = 4.5 * Q1Third(i) - 4.5 * Q1Half(i) + Q1(i)
        i += 1
      }

      //      QL_ = QL_HalfPart
      //      QR_ = QR_HalfPart
      //
      //      var i = 0
      //      while (i < size) {
      //        Q1(i) = Q1Half(i)
      //        i += 1
      //      }

      //      var sumQ = QR_ + QL_
      //      var sumF = Fmin * QL_ + Fmax * QR_
      //      i = size - 2
      //      while (i > 0) {
      //        sumQ += h * Q1(i)
      //        sumF += h * Q1(i) * F(i)
      //        i -= 1
      //      }
      //      println("t=" + t + " LMG3 Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0)+" QL_="+QL+" QR_="+QR)

      //      M0 = M1
      Qtmp = Q0_
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
package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganBDF3SABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    M1 = Array.ofDim[Double](size)
    buildMcache(dt, 0, M1)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    var Q2 = Array.ofDim[Double](size)
    var Q3 = Array.ofDim[Double](size)
    val rhs = Array.ofDim[Double](size)
    QL_ = 0.0
    var QL_0, QR_0, QL_1, QR_1 = 0.0
    QR_ = 0.0
    var t = T
    //    var rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //       if (tIndex < 4) {
      //        printTimeStep(t)
      //      }
      if (tIndex == 0) {
        //Euler
        t -= dt
        advance(dt, M1)
        computeSystem(tri1, dt)
        Q0_(0) = 0
        Q0_(size - 1) = 0
        solver.solve(tri1, Q0_, Q1)
        QL_0 = QL_
        QR_0 = QR_
        QL_ += dt / (h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
        QR_ -= dt / (h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))
      } else if (tIndex == 1) {
        //BDF2
        t -= dt
        advance(dt, M1)
        computeSystem(tri1, 2 * dt / 3)

        var j = size - 2
        while (j > 0) {
          rhs(j) = (4 * Q1(j) - Q0_(j)) / 3
          j -= 1
        }
        solver.solve(tri1, rhs, Q2)
        QL_1 = QL_
        QR_1 = QR_
        QL_ = (4 * QL_ - QL_0) / 3 + 2 * dt / 3 / (h) * (M1(1) * Q2(1) - M1(0) * Q2(0))
        QR_ = (4 * QR_ - QR_0) / 3 - 2 * dt / 3 / (h) * (M1(size - 1) * Q2(size - 1) - M1(size - 2) * Q2(size - 2))
        //        var sumQ = QR_ + QL_
        //        var sumF = Fmin * QL_ + Fmax * QR_
        //        var i = size - 2
        //        while (i > 0) {
        //          sumQ += h * Q2(i)
        //          sumF += h * Q2(i) * F(i)
        //          i -= 1
        //        }
        //        println("t=" + t + " BDF2 Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0))

      } else {
        //BDF3
        t -= dt
        advance(dt, M1)
        computeSystem(tri1, 6 * dt / 11)

        var j = size - 2
        while (j > 0) {
          rhs(j) = (18 * Q2(j) - 9 * Q1(j) + 2 * Q0_(j)) / 11
          j -= 1
        }
        solver.solve(tri1, rhs, Q3)
        var QL_temp = QL_
        var QR_temp = QR_
        QL_ = (18 * QL_ - 9 * QL_1 + 2 * QL_0) / 11 + 6 * dt / 11 / (h) * (M1(1) * Q3(1) - M1(0) * Q3(0))
        QR_ = (18 * QR_ - 9 * QR_1 + 2 * QR_0) / 11 - 6 * dt / 11 / (h) * (M1(size - 1) * Q3(size - 1) - M1(size - 2) * Q3(size - 2))

        //        var sumQ = QR_ + QL_
        //        var sumF = Fmin * QL_ + Fmax * QR_
        //        var i = size - 2
        //        while (i > 0) {
        //          sumQ += h * Q3(i)
        //          sumF += h * Q3(i) * F(i)
        //          i -= 1
        //        }
        //        println("t=" + t + " BDF3 Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0))

        QL_0 = QL_1
        QR_0 = QR_1
        QL_1 = QL_temp
        QR_1 = QR_temp
        var Qtmp = Q0_
        Q0_ = Q1
        Q1 = Q2
        Q2 = Q3
        Q3 = Qtmp
      }
      tIndex += 1

    }
    Q0_ = Q2
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
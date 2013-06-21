package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.{Epsilon}


class HaganADESABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {
  override def solve() {

    M1 = Array.ofDim[Double](size)
    M0 = Array.ofDim[Double](size)
    buildMcache(dt, 0, M0)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    var Q2 = Array.ofDim[Double](size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    var rhs = Array.ofDim[Double](size)
    var tIndex = 0
    var r = dt / (h * h)

    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //                        if (tIndex < 4) {
      //                          printTimeStep(t)
      //                          tIndex += 1
      //                        }
      t -= dt
      Array.copy(M0, 0, M1, 0, size)
      advance(dt, M1)

      //LR
      var i = 1
      Q1(0) = -Q0(1) * M0(1) / M1(0)
      while (i < size - 1) {
        Q1(i) = (r * M1(i - 1) * Q1(i - 1) + (1 - M0(i) * r) * Q0_(i) + r * M0(i + 1) * Q0_(i + 1)) / (1 + r * M1(i))
        i += 1
      }
      Q1(size - 1) = -Q0_(size - 2) * M0(size - 2) / M1(size - 1) //or opposite?

      val QL1 = dt / h * (M1(0) * Q1(0) - M0(1) * Q0_(1))
      val QR1 = dt / h * (M0(size - 1) * Q0_(size - 1) - M1(size - 2) * Q1(size - 2))

      var sumQ = -QR1 - QL1 + QR_ + QL_
      var sumF = Fmin * (QL_ - QL1) + Fmax * (QR_ - QR1)
      i = size - 2
      while (i > 0) {
        sumQ += h * Q1(i)
        sumF += h * Q1(i) * F(i)
        i -= 1
      }
      println("t=" + t + " LR Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0) + " QL_=" + (QL_ - QL1) + " QR_=" + (QR_ - QR1))

      //RL
      Q2(size - 1) = -Q0_(size - 2) * M0(size - 2) / M1(size - 1)
      i = size - 2
      while (i > 0) {
        Q2(i) = (r * M1(i + 1) * Q2(i + 1) + (1 - M0(i) * r) * Q0_(i) + r * M0(i - 1) * Q0_(i - 1)) / (1 + r * M1(i))
        i -= 1
      }
      Q2(0) = -Q0(1) * M0(1) / M1(0) //or opposite?

      val QL2 = dt / h * (M0(0) * Q0(0) - M1(1) * Q2(1))
      val QR2 = dt / h * (M1(size - 1) * Q2(size - 1) - M0(size - 2) * Q0(size - 2))

      i = 0
      while (i < size - 1) {
        Q0_(i) = 0.5 * (Q1(i) + Q2(i))
        i += 1
      }
      QL_ -= 0.5 * (QL1 + QL2)
      QR_ -= 0.5 * (QR1 + QR2)

      sumQ = QR_ + QL_
      sumF = Fmin * QL_ + Fmax * QR_
      i = size - 2
      while (i > 0) {
        sumQ += h * Q0(i)
        sumF += h * Q0(i) * F(i)
        i -= 1
      }
      println("t=" + t + " AV Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0) + " QL_=" + QL_ + " QR_=" + QR_)

      val Mtmp = M0
      M0 = M1
      M1 = Mtmp
    }
  }
}

package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.{Epsilon, TridiagonalMatrix}

class HaganCNBDF3SABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  override def solve() {
    val dt3 = dt/3.0
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    M0 = Array.ofDim[Double](size)
    M1 = Array.ofDim[Double](size)
    var Q2 = Array.ofDim[Double](size)
    buildMcache(dt3, 0, M0)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0

    var t = T
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0
    val cdth = 0.5*dt3/(h)
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
//                  if (tIndex < 4) {
//              printTimeStep(t)
//              tIndex += 1
//            }

      t -= dt
      Array.copy(M0, 0, M1, 0, size)
      advance(dt3, M1)
      computeSystem(tri1, 0.5*dt3)
      tri0.multiply(Q0_, rhs)
      //      rhs(0) = 0
      //      rhs(size - 1) = 0
      solver.solve(tri1, rhs, Q1)
      val QL1 = QL_
      val QR1 = QR_
      QL_ += cdth * (M1(1) * Q1(1) - M1(0) * Q1(0) + M0(1) * Q0_(1) - M0(0) * Q0_(0))
      QR_ -= cdth * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2) + M0(size - 1) * Q0_(size - 1) - M0(size - 2) * Q0_(size - 2))


      Array.copy(M1, 0, M0, 0, size)
      advance(dt3, M1)
      computeSystem(tri1, 0.5*dt3)
      tri0.multiply(Q1, rhs)
      //      rhs(0) = 0
      //      rhs(size - 1) = 0
      solver.solve(tri1, rhs, Q2)
      val QL2 = QL_
      val QR2 = QR_
      QL_ += cdth * (M1(1) * Q2(1) - M1(0) * Q2(0) + M0(1) * Q1(1) - M0(0) * Q1(0))
      QR_ -= cdth * (M1(size - 1) * Q2(size - 1) - M1(size - 2) * Q2(size - 2) + M0(size - 1) * Q1(size - 1) - M0(size - 2) * Q1(size - 2))

//                        var sum = QR_ + QL_
//
//                        var i = size-2
//                        while (i>0) {
//                          sum += h*Q2(i)
//                          //sumQ0_ += h*Q0_(i)
//                          i-=1
//                        }
//                        println("t="+t+" TR Q="+sum)

      var j = size - 2
      while (j >= 1) {
        rhs(j) =  (18 * Q2(j) - 9 * Q1(j) + 2 * Q0_(j)) / 11
        j -= 1
      }
      advance(dt3, M1)
      computeSystem(tri1, dt3*6.0/11)
      solver.solve(tri1, rhs, Q1)
      //      QL_Part += dt / (2*h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
      //      QR_Part += dt / (2*h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

      QL_ = (18 * QL_ - 9 * QL2 + 2 * QL1) / 11 + 6 * dt3 / 11 / (h) *(M1(1) * Q1(1) - M1(0) * Q1(0))
      QR_ = (18 * QR_ - 9 * QR2 + 2 * QR1) / 11 - 6 * dt3 / 11 / (h) *  (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

//                        var sumQ = QR_ + QL_
//                        var sumF = Fmin * QL_ + Fmax * QR_
//             i = size - 2
//                        while (i > 0) {
//                          sumQ += h * Q1(i)
//                          sumF += h * Q1(i) * F(i)
//                          i -= 1
//                        }
//                        println("t=" + t + " TRBDF3 Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0))

      val Mtmp = M0
      M0 = M1
      M1 = Mtmp
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
      tri1 middle(j) = 1 + 2 * frac * M1(j)
      tri1.upper(j) = -frac * M1(j + 1)
      tri0.lower(j) = frac * M0(j - 1)
      tri0 middle(j) = 1 - 2 * frac * M0(j)
      tri0.upper(j) = frac * M0(j + 1)
      j += 1
    }
    tri1.upper(0) = M1(1)
    tri1 middle(0) = M1(0)
    tri1.lower(size - 1) = M1(size - 2)
    tri1 middle(size - 1) = M1(size - 1)
  }
}

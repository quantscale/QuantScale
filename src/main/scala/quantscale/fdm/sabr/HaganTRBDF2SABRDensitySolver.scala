package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganTRBDF2SABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = 0.5 * alpha;
  private val frac = 1.0 / (alpha * (2 - alpha));

  override def solve() {
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    M0 = Array.ofDim[Double](size)
    M1 = Array.ofDim[Double](size)
    val bdt = dt * backcoeff
    buildMcache(2 * bdt, dt - 2 * bdt, M0)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    val rhs = Array.ofDim[Double](size)
    val cdth = backcoeff * dt / (h)
    var tIndex = 0

    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //                  if (tIndex < 4) {
      //              printTimeStep(t)
      //              tIndex += 1
      //            }

      t -= dt
      Array.copy(M0, 0, M1, 0, size)
      advance(2 * bdt, M1)
      computeSystem(tri1, bdt)
      tri0.multiply(Q0_, rhs)
      //      rhs(0) = 0
      //      rhs(size - 1) = 0
      solver.solve(tri1, rhs, Q1)
      val QL_init = QL_
      val QR_init = QR_
      QL_ += cdth * (M1(1) * Q1(1) - M1(0) * Q1(0) + M0(1) * Q0_(1) - M0(0) * Q0_(0))
      QR_ -= cdth * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2) + M0(size - 1) * Q0_(size - 1) - M0(size - 2) * Q0_(size - 2))

      //                  var sum = QR_+QL
      //                  var sumQ0_ = QL_init+QRinit
      //                  var i = size-2
      //                  while (i>0) {
      //                    sum += h*Q1(i)
      //                    sumQ0_ += h*Q0_(i)
      //                    i-=1
      //                  }
      //                  println("t="+t+" TR Q="+sum+" "+sumQ0_)

      var j = size - 2
      while (j >= 1) {
        rhs(j) = frac * (Q1(j) - (1 - alpha) * (1 - alpha) * Q0_(j))
        j -= 1
      }
      advance(dt - 2 * bdt, M1)
      computeSystem(tri1, bdt)
      solver.solve(tri1, rhs, Q1)
      //      QL_Part += dt / (2*h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
      //      QR_Part += dt / (2*h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

      QL_ = frac * (QL - (1 - alpha) * (1 - alpha) * QL_init) + cdth * (M1(1) * Q1(1) - M1(0) * Q1(0))
      QR_ = frac * (QR - (1 - alpha) * (1 - alpha) * QR_init) - cdth * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

      //            var sumQ = QR_ + QL_
      //            var sumF = Fmin * QL_ + Fmax * QR_
      // i = size - 2
      //            while (i > 0) {
      //              sumQ += h * Q1(i)
      //              sumF += h * Q1(i) * F(i)
      //              i -= 1
      //            }
      //            println("t=" + t + " TRBDF2 Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0))

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
      tri1.middle(j) = 1 + 2 * frac * M1(j)
      tri1.upper(j) = -frac * M1(j + 1)
      tri0.lower(j) = frac * M0(j - 1)
      tri0.middle(j) = 1 - 2 * frac * M0(j)
      tri0.upper(j) = frac * M0(j + 1)
      j += 1
    }
    tri1.upper(0) = M1(1)
    tri1.middle(0) = M1(0)
    tri1.lower(size - 1) = M1(size - 2)
    tri1.middle(size - 1) = M1(size - 1)
  }
}
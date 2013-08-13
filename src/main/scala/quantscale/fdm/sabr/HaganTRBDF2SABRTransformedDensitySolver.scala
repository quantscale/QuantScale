package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganTRBDF2SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = 0.5 * alpha;
  private val frac = 1.0 / (alpha * (2 - alpha));

  override def solve() {
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val M0 = Array.ofDim[Double](size)
    val bdt = 2*dt * backcoeff
    buildEmCache(bdt, dt - bdt)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    val rhs = Array.ofDim[Double](size)
    val cdth = backcoeff * dt / (h)
    var tIndex = 0

    while (tIndex < timeSteps) {
//                  if (tIndex < 4) {
//              printTimeStep(t)
//              tIndex += 1
//            }

      t -= dt
      Array.copy(Em_, 0, M0, 0, size)
      advanceEm(bdt, Em_)
      computeSystem(bdt, Em_, M0, tri1, tri0)
      tri0.multiply(Q0_, rhs)
      solver.solve(tri1, rhs, Q1)
      val QL_init = QL_
      val QR_init = QR_
      QL_ += 0.5*bdt * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * Q1(1) + M0(1) * Q0_(1))
      QR_ += 0.5*bdt * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * Q1(size - 2) + M0(size - 2) * Q0_(size - 2))
//      printSumQF("TR", t, Q1, QL_, QR_)

      var j = size - 2
      while (j >= 1) {
        rhs(j) = frac * (Q1(j) - (1 - alpha) * (1 - alpha) * Q0_(j))
        j -= 1
      }
      advanceEm(dt - bdt, Em_)
      computeSystem(bdt, Em_, null, tri1, tri0)
      solver.solve(tri1, rhs, Q1)

      QL_ = frac * (QL - (1 - alpha) * (1 - alpha) * QL_init) + 0.5*bdt * computedQLdt(Em_, Q1)
      QR_ = frac * (QR - (1 - alpha) * (1 - alpha) * QR_init) + 0.5*bdt * computedQRdt(Em_, Q1)

//      printSumQF("BDF2 ", t, Q1, QL_, QR_)

      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp

      tIndex+=1
    }
  }

}
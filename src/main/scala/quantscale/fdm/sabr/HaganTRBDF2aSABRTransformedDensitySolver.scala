package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.{Epsilon, TridiagonalMatrix}

class HaganTRBDF2aSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = 0.5 * alpha;
  private val frac = 1.0 / (alpha * (2 - alpha));
  var relativeAccuracy = 1e-3
  var absoluteAccuracy = 1e-6
  private var timeStepCounter_ = 0

  private def buildEm(t: Double) {
    var j = size - 2
    val efactor = spec.rho * spec.nu * spec.alpha
    while (j > 0) {
      Em_(j) = math.exp(efactor * Gammam_(j) * (T - t))
      j -= 1
    }
    Em_(0) = Em_(1)
    Em_(size - 1) = Em_(size - 2)
  }

  override def solve() {
    timeStepCounter_ = 0
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    Em_ = Array.fill[Double](size)(1.0)

    val M0 = Array.ofDim[Double](size)

    //buildEmCache(bdt, dt - bdt)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    val L0 = Array.ofDim[Double](size)
    var L1 = Array.ofDim[Double](size)
    var L2 = Array.ofDim[Double](size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    val dtmin = T / 10000
    val dtmax = T / timeSteps //T/2
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0
    val k = (-3 * alpha * alpha + 4 * alpha - 2) / (12 * (2 - alpha))

    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //                        if (tIndex < 4) {
      //                    printTimeStep(t)
      //                  }
      val bdt = 2 * dt_ * backcoeff
      t -= bdt
      Array.copy(Em_, 0, M0, 0, size)
      buildEm(t)
      //advanceEm(bdt, Em_)
      computeSystem(bdt, Em_, M0, tri1, tri0)
      tri0.multiply(P0_, rhs)
      solver.solve(tri1, rhs, P1_)
      var j = size - 2
      while (j >= 1) {
        L0(j) = -(tri0.lower(j) * P0_(j - 1) + (tri0.middle(j) - 1) * P0_(j) + tri0.upper(j) * P0_(j + 1)) / bdt * 2
        L1(j) = -(tri1.lower(j) * P1_(j - 1) + (tri1.middle(j) - 1) * P1_(j) + tri1.upper(j) * P1_(j + 1)) / bdt * 2
        j -= 1
      }
      val QL_init = PL_
      val QR_init = PR_
      PL_ += 0.5 * bdt * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * P1_(1) + M0(1) * P0_(1))
      PR_ += 0.5 * bdt * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * P1_(size - 2) + M0(size - 2) * P0_(size - 2))
      //      printSumQF("TR", t, Q1, QL_, QR_)

      j = size - 2
      while (j >= 1) {
        rhs(j) = frac * (P1_(j) - (1 - alpha) * (1 - alpha) * P0_(j))
        j -= 1
      }
      t -= dt_ - bdt
      buildEm(t)

      //advanceEm(dt - bdt, Em_)
      computeSystem(bdt, Em_, null, tri1, tri0)
      solver.solve(tri1, rhs, P1_)
      j = size - 2
      while (j >= 1) {
        L2(j) = -(tri1.lower(j) * P1_(j - 1) + (tri1.middle(j) - 1) * P1_(j) + tri1.upper(j) * P1_(j + 1)) / bdt * 2
        j -= 1
      }
      PL_ = frac * (PL - (1 - alpha) * (1 - alpha) * QL_init) + 0.5 * bdt * computedPLdt(Em_, P1_)
      PR_ = frac * (PR - (1 - alpha) * (1 - alpha) * QR_init) + 0.5 * bdt * computedPRdt(Em_, P1_)

      //      printSumQF("BDF2 ", t, Q1, QL_, QR_)

      var tau = 0.0
      var error = 0.0
      var rsq = 0.0
      j = size - 2
      while (j >= 1) {
        val tauj = 2 * k * dt_ * (1 / alpha * L0(j) - 1 / (alpha * (1 - alpha)) * L1(j) + 1 / (1 - alpha) * L2(j))
        val errorj = relativeAccuracy * math.abs(P1_(j)) + absoluteAccuracy
        tau += tauj * tauj
        error += errorj * errorj
        rsq += tauj * tauj / (errorj * errorj)
        j -= 1
      }
      val ratio = math.sqrt(rsq) / (size - 2)
      val r = math.abs(ratio)
      var newdt = dt_ * math.pow(r, -1.0 / 3.0)
      if (r < 2 || dt_ == dtmin) {
        val Qtmp = P0_
        P0_ = P1_
        P1_ = Qtmp
      } else {
        t = t + dt_ //rollback
        newdt *= 0.9
        Array.copy(M0, 0, Em_, 0, size)
        PL_ = QL_init
        PR_ = QR_init
      }
      dt_ = math.max(math.min(newdt, dtmax), dtmin)
      if (dt_ > t) dt_ = t
      println("t=" + t + " dt_=" + dt_ + " tau=" + tau + " error=" + error + " r=" + r + " Pj0=" + P1_(j0))


      tIndex += 1
      timeStepCounter_ += 1
    }
    println("-- timeStepCounter=" + timeStepCounter_)
  }

}
package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix

class HaganRichardsonEuler3SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  private var Em3Cache: Array[Double] = null

  override def solve() {
    buildEmCache()
    tri1 = new TridiagonalMatrix(size)
    val EmInit = Em_.clone()
    solve(1)
    val P1_1 = P0_
    val QL_1 = PL_
    val QR_1 = PR_
    Array.copy(EmInit, 0, Em_, 0, size)
    solve(2)
    val P1_2 = P0_
    val QL_2 = PL_
    val QR_2 = PR_
    Array.copy(EmInit, 0, Em_, 0, size)
    solve(3)


    var i = 0
    while (i < size) {
      P0_(i) = 4.5 * P0_(i) - 4 * P1_2(i) + 0.5 * P1_1(i)
      i += 1
    }
    PL_ = 4.5 * PL_ - 4 * QL_2 + 0.5 * QL_1
    PR_ = 4.5 * PR_ - 4 * QR_2 + 0.5 * QR_1

  }

  private[sabr] def buildEmCache() {
    Em_ = Array.ofDim[Double](size)
    Em1Cache = Array.ofDim[Double](size)
    Em2Cache = Array.ofDim[Double](size)
    Em3Cache = Array.ofDim[Double](size)
    var j = size - 2
    val efactor = spec.rho * spec.nu * spec.alpha
    while (j > 0) {
      Em2Cache(j) = math.exp(efactor * Gammam_(j) * dt / 2)
      Em3Cache(j) = math.exp(efactor * Gammam_(j) * dt / 3)

      Em1Cache(j) = Em2Cache(j) * Em2Cache(j)

      Em_(j) = 1.0
      j -= 1
    }
    Em_(0) = Em_(1)
    Em_(size - 1) = Em_(size - 2)
    Em1Cache(0) = Em1Cache(1)
    Em1Cache(size - 1) = Em1Cache(size - 2)
    Em2Cache(0) = Em2Cache(1)
    Em2Cache(size - 1) = Em2Cache(size - 2)
    Em3Cache(0) = Em3Cache(1)
    Em3Cache(size - 1) = Em3Cache(size - 2)

  }

  private[sabr] def advanceEm(divisor: Int, M: Array[Double]) {
    var j = 0
    var Mcache: Array[Double] = null
    divisor match {
      case 1 => Mcache = Em1Cache
      case 2 => Mcache = Em2Cache
      case 3 => Mcache = Em3Cache
    }
    while (j < size) {
      M(j) = M(j) * Mcache(j)
      j += 1
    }
  }

  def solve(divisor: Int) {
    val ddt = dt_ / divisor
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    var tIndex = 0
    while (tIndex < timeSteps * divisor) {
      t -= ddt
      advanceEm(divisor, Em_)
      computeSystem(ddt, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, P1_)
      PL_ += ddt * computedPLdt(Em_, P1_)
      PR_ += ddt * computedPRdt(Em_, P1_)
      //   printSumQF("RE",t,Q1, QL_, QR_)
      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }

  }


}
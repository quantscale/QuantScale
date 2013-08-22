package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix

class HaganBDF3SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  var useTRinit = false

  override def solve() {
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val M0 = Array.ofDim[Double](size)
    buildEmCache(dt, 0)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    var P2 = Array.ofDim[Double](size)
    var P3 = Array.ofDim[Double](size)
    PL_ = 0.0
    PR_ = 0.0
    var PL0, PL1 = PL_
    var PR0, PR1 = PR_

    var t = T
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0

    while (tIndex < timeSteps) {
      //                        if (tIndex < 4) {
      //                    printTimeStep(t)
      //                  }

      t -= dt
      if (tIndex == 0) {
        if (useTRinit) {
          Array.copy(Em_, 0, M0, 0, size)
          advanceEm(dt, Em_)
          computeSystem(dt, Em_, M0, tri1, tri0)
          tri0.multiply(P0_, rhs)
          solver.solve(tri1, rhs, P1_)
          PL_ += 0.5 * dt * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * P1_(1) + M0(1) * P0_(1))
          PR_ += 0.5 * dt * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * P1_(size - 2) + M0(size - 2) * P0_(size - 2))
        } else {
          advanceEm(dt, Em_)
          computeSystem(dt, Em_, tri1)
          solver.solve(tri1, P0_, P1_)
          PL_ += dt * computedPLdt(Em_, P1_)
          PR_ += dt * computedPRdt(Em_, P1_)
        }
        //printSumPF("BDF1", t, P1_, PL_, PR_)
      } else if (tIndex == 1) {
        advanceEm(dt, Em_)
        computeSystem(2 * dt / 3.0, Em_, tri1)
        var j = size - 2
        while (j >= 1) {
          rhs(j) = (4 * P1_(j) - P0_(j)) / 3
          j -= 1
        }
        solver.solve(tri1, rhs, P2)
        PL1 = PL_
        PR1 = PR_
        PL_ = (4 * PL_ - PL0) / 3 + 2 / 3.0 * dt * computedPLdt(Em_, P2)
        PR_ = (4 * PR_ - PR0) / 3 + 2 / 3.0 * dt * computedPRdt(Em_, P2)

        //printSumPF("BDF2 ", t, P2, PL_, PR_)
      } else {
        advanceEm(dt, Em_)
        computeSystem(6 * dt / 11.0, Em_, tri1)
        var j = size - 2
        while (j >= 1) {
          rhs(j) = (18 * P2(j) - 9 * P1_(j) + 2 * P0_(j)) / 11
          j -= 1
        }
        solver.solve(tri1, rhs, P3)
        val PL_temp = PL_
        val PR_temp = PR_
        PL_ = (18 * PL_ - 9 * PL1 + 2 * PL0) / 11 + 6 / 11.0 * dt * computedPLdt(Em_, P3)
        PR_ = (18 * PR_ - 9 * PR1 + 2 * PR0) / 11 + 6 / 11.0 * dt * computedPRdt(Em_, P3)
        PL0 = PL1
        PR0 = PR1
        PL1 = PL_temp
        PR1 = PR_temp
        //printSumPF("BDF2 ", t, P2, PL_, PR_)
        val Ptmp = P0_
        P0_ = P1_
        P1_ = P2
        P2 = P3
        P3 = Ptmp
      }
      tIndex += 1
    }
    if (timeSteps == 1) {
      P0_ = P1_
    } else {
      P0_ = P2
    }
  }

}
package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix

class HaganIMBDF3SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  val gamma = 0.158983899989
  val beta20 = -1.644972541469
  val beta21 = 2.644972541469
  val beta30 = 3.911869324372
  val beta31 = -6.012776528027
  val beta32 = 3.100907203655

  override def solve() {
    val gdt = gamma * dt
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val Q2 = Array.ofDim[Double](size)
    buildEmCache(gdt, dt - 3 * gdt)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0

    var t = T
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (tIndex < timeSteps) {
      //                        if (tIndex < 4) {
      //                    printTimeStep(t)
      //                  }

      t -= dt
      advanceEm(gdt, Em_)
      computeSystem(gdt, Em_, tri1)
      solver.solve(tri1, P0_, P1_)
      val QL1 = PL_
      val QR1 = PR_
      PL_ += gdt * computedPLdt(Em_, P1_)
      PR_ += gdt * computedPRdt(Em_, P1_)
      //    printSumPF("IMBDF1/3", t, P1_, PL_, PR_)

      var j = size - 2
      while (j >= 1) {
        rhs(j) = beta21 * P1_(j) + beta20 * P0_(j)
        j -= 1
      }
      advanceEm(gdt, Em_)
      computeSystem(gdt, Em_, tri1)
      solver.solve(tri1, rhs, Q2)
      val QL2 = PL_
      val QR2 = PR_
      PL_ = (beta21 * PL_ + beta20 * QL1) + gdt * computedPLdt(Em_, Q2)
      PR_ = (beta21 * PR_ + beta20 * QR1) + gdt * computedPRdt(Em_, Q2)

      //printSumPF("IMBDF2/3", t, Q2, PL_, PR_)

      j = size - 2
      while (j >= 1) {
        rhs(j) = (beta32 * Q2(j) + beta31 * P1_(j) + beta30 * P0_(j))
        j -= 1
      }
      advanceEm(gdt, Em_)
      computeSystem(gdt, Em_, tri1)
      solver.solve(tri1, rhs, P1_)
      P1_(0) = 0
      P1_(size - 1) = 0
      PL_ = (beta32 * PL_ + beta31 * QL2 + beta30 * QL1) + gdt * computedPLdt(Em_, P1_)
      PR_ = (beta32 * PR_ + beta31 * QR2 + beta30 * QR1) + gdt * computedPRdt(Em_, P1_)
      advanceEm(dt - 3 * gdt, Em_)
      //printSumPF("IMBDF3", t, P1_, PL_, PR_)
      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }
  }

}

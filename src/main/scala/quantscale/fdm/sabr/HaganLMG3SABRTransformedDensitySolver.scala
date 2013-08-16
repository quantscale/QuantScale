package quantscale.fdm.sabr

import quantscale.fdm.{Epsilon, TridiagonalMatrix}
import quantscale.analytic.SABRModelSpec

/**
 * Date: 8/12/13
 */
class HaganLMG3SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    buildEmCache(dt / 3, 0)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    val Q1Third = Array.ofDim[Double](size)
    val Q1ThirdTmp = Array.ofDim[Double](size)
    val Q1Half = Array.ofDim[Double](size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    var Qtmp = P0_
    var tIndex = 0

    while (tIndex < timeSteps) {
//                        if (tIndex < 4) {
//                    printTimeStep(t)
//                  }
      t -= dt / 3

      advanceEm(dt / 3, Em_)
      computeSystem(dt / 3, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1ThirdTmp)
      var QL_ThirdPart = PL_ + dt / 3 * computedPLdt(Em_, Q1ThirdTmp)
      var QR_ThirdPart = PR_ + dt / 3 * computedPRdt(Em_, Q1ThirdTmp)
      var QL_HalfPart = QL_ThirdPart
      var QR_HalfPart = QR_ThirdPart
      val Q0_Init = P0_
      P0_ = Q1ThirdTmp

      t -= dt / 3
      advanceEm(dt / 3, Em_)
      computeSystem(dt / 3, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1Third)
      QL_ThirdPart += dt / 3 * computedPLdt(Em_, Q1Third)
      QR_ThirdPart += dt / 3 * computedPRdt(Em_, Q1Third)
      P0_ = Q1Third

      t -= dt / 3
      advanceEm(dt / 3, Em_)
      computeSystem(dt / 3, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1Third)
      QL_ThirdPart += dt / 3 * computedPLdt(Em_, Q1Third)
      QR_ThirdPart += dt / 3 * computedPRdt(Em_, Q1Third)

      P0_ = Q1ThirdTmp
      computeSystem(2 * dt / 3, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1Half)
      QL_HalfPart += 2 * dt / 3 * computedPLdt(Em_, Q1Half)
      QR_HalfPart += 2 * dt / 3 * computedPRdt(Em_, Q1Half)

      P0_ = Q0_Init
      computeSystem(dt, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, P1_)

      PL_ += 4.5 * QL_ThirdPart - 4.5 * QL_HalfPart + dt * computedPLdt(Em_, P1_)
      PR_ += 4.5 * QR_ThirdPart - 4.5 * QR_HalfPart + dt * computedPRdt(Em_, P1_)

      var i = 0
      while (i < size) {
        P1_(i) = 4.5 * Q1Third(i) - 4.5 * Q1Half(i) + P1_(i)
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
      Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }
  }

}

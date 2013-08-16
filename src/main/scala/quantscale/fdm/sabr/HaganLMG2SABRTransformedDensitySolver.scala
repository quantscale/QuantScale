package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganLMG2SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    buildEmCache(dt / 2, 0)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    val Q1Half = Array.ofDim[Double](size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    //    var rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (tIndex < timeSteps) {
//             if (tIndex < 4) {
//              printTimeStep(t)
//            }
      t -= dt / 2
      //      System.arraycopy(M0, 0, M1, 0, size)
      advanceEm(dt / 2, Em_)
      computeSystem(dt / 2, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1Half)
      var QL_Part = dt / 2 * computedPLdt(Em_, Q1Half)
      var QR_Part = dt / 2 * computedPRdt(Em_, Q1Half)
      val Q0_Init = P0_
      P0_ = Q1Half

      t -= dt / 2
      advanceEm(dt / 2, Em_)
      computeSystem(dt / 2, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1Half)
      QL_Part += dt / 2 * computedPLdt(Em_, Q1Half)
      QR_Part += dt / 2 * computedPRdt(Em_, Q1Half)

      P0_ = Q0_Init
      computeSystem(dt, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, P1_)

      PL_ += 2 * QL_Part - dt * computedPLdt(Em_, P1_)
      PR_ += 2 * QR_Part - dt * computedPRdt(Em_, P1_)

      var i = 0
      while (i < size) {
        P1_(i) = 2 * Q1Half(i) - P1_(i)
        i += 1
      }
      //printSumQF("LMG2",t,Q1,QL_,QR_)
      //check sum
      //      var sum = QR_+QL
      //      i = size-2
      //      while (i>0) {
      //        sum += h*Q1(i)
      //        i-=1
      //      }
      //      println("t="+t+" LMG2 Q="+sum)
      //      M0 = M1
      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }

  }

}
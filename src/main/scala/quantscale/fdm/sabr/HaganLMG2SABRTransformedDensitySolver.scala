package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganLMG2SABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    buildEmCache(dt / 2, 0)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    val Q1Half = Array.ofDim[Double](size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    //    var rhs = Array.ofDim[Double](size)
    var tIndex = 0
    while (tIndex < timeSteps) {
      //       if (tIndex < 4) {
      //        printTimeStep(t)
      //        tIndex += 1
      //      }
      t -= dt / 2
//      System.arraycopy(M0, 0, M1, 0, size)
      advanceEm(dt / 2, Em_)
      computeSystem(dt / 2, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      var QL_Part = dt / 2  * computedQLdt(Em_, Q1Half)
      var QR_Part = dt / 2 * computedQRdt(Em_, Q1Half)
      val Q0_Init = Q0_
      Q0_ = Q1Half

      t -= dt / 2
      advanceEm(dt / 2, Em_)
      computeSystem(dt / 2, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      QL_Part += dt / 2 * computedQLdt(Em_, Q1Half)
      QR_Part += dt / 2 * computedQRdt(Em_, Q1Half)

      Q0_ = Q0_Init
      computeSystem(dt, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1)

      QL_ += 2 * QL_Part - dt * computedQLdt(Em_, Q1)
      QR_ += 2 * QR_Part - dt * computedQRdt(Em_, Q1)

      var i = 0
      while (i < size) {
        Q1(i) = 2 * Q1Half(i) - Q1(i)
        i += 1
      }

      //check sum
      //      var sum = QR_+QL
      //      i = size-2
      //      while (i>0) {
      //        sum += h*Q1(i)
      //        i-=1
      //      }
      //      println("t="+t+" LMG2 Q="+sum)
//      M0 = M1
      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
      tIndex += 1
    }

  }

}
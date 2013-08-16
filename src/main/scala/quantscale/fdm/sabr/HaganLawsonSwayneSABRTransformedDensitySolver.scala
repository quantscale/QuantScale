package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix

/* Date: 8/12/13
 */
class HaganLawsonSwayneSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  private val b = 1 - math.sqrt(2) / 2
  private val sqrt2 = math.sqrt(2)

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    buildEmCache(dt * b, dt * (1 - 2 * b))
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    val Q1Half = Array.ofDim[Double](size)
    PL_ = 0.0
    var QL_Full, QL_Half, Q1HalfStart, Q1HalfEnd = 0.0
    var QR_Full, QR_Half = 0.0
    PR_ = 0.0
    var t = T
    var tIndex = 0
    var Qtmp = P0_

    while (tIndex < timeSteps) {
//                   if (tIndex < 4) {
//                    printTimeStep(t)
//                  }
      t -= dt
      advanceEm(dt * b, Em_)
      computeSystem(dt * b, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, Q1Half)
      QL_Full = PL_ + dt * b * computedPLdt(Em_, Q1Half)
      QR_Full = PR_ + dt * b * computedPRdt(Em_, Q1Half)
      QL_Half = QL_Full
      QR_Half = QR_Full
      // //check sum
      //                        var sum = QR_Half + QL_Half
      //                  var i = size - 2
      //                  while (i > 0) {
      //                    sum += h * Q1Half(i)
      //                    i -= 1
      //                  }
      //                  println("t=" + t + " LawsonSwayne-1 Q=" + sum+ "Q(f)="+Q1Half(j0)+" QL="+QL_Half+" QR="+QR_Half)

      //      t -= dt *b
      advanceEm(dt * b, Em_)
      computeSystem(dt * b, Em_, tri1)
      Q1Half(0) = 0
      Q1Half(size - 1) = 0
      solver.solve(tri1, Q1Half, P1_)
      QL_Full += dt * b * computedPLdt(Em_, P1_)
      QR_Full += dt * b * computedPRdt(Em_, P1_)

      //      sum = QR_Full + QL_Full
      //      i = size - 2
      //      while (i > 0) {
      //        sum += h * Q1(i)
      //        i -= 1
      //      }
      //      println("t=" + t + " LawsonSwayne-1.5 Q=" + sum)

      var i = 0
      while (i < size) {
        P1_(i) = (sqrt2 + 1) * P1_(i) - sqrt2 * Q1Half(i)
        i += 1
      }
      PL_ = (sqrt2 + 1) * QL_Full - sqrt2 * QL_Half
      PR_ = (sqrt2 + 1) * QR_Full - sqrt2 * QR_Half
      advanceEm(dt * (1 - 2 * b), Em_)
      //printSumQF(t)
      Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }

  }

}

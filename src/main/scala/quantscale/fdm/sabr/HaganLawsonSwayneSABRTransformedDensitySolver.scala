package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix

/* Date: 8/12/13
 */
class HaganLawsonSwayneSABRTransformedDensitySolver (spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  private val b = 1 - math.sqrt(2) / 2
  private val sqrt2 = math.sqrt(2)

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    buildEmCache(dt *b, dt*(1-2*b))
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    val Q1Half = Array.ofDim[Double](size)
    QL_ = 0.0
    var QL_Full, QL_Half, Q1HalfStart, Q1HalfEnd = 0.0
    var QR_Full, QR_Half = 0.0
    QR_ = 0.0
    var t = T
    var tIndex = 0
    var Qtmp = Q0_

    while (tIndex < timeSteps) {
      //             if (tIndex < 4) {
      //              printTimeStep(t)
      //              tIndex += 1
      //            }
      t -= dt
      advanceEm(dt*b, Em_)
      computeSystem(dt * b, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      QL_Full = QL_ + dt * b * computedQLdt(Em_, Q1Half)
      QR_Full = QR_ + dt * b * computedQRdt(Em_, Q1Half)
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
      advanceEm(dt*b, Em_)
      computeSystem(dt * b, Em_, tri1)
      Q1Half(0) = 0
      Q1Half(size - 1) = 0
      solver.solve(tri1, Q1Half, Q1)
      QL_Full += dt * b * computedQLdt(Em_, Q1)
      QR_Full += dt * b * computedQRdt(Em_, Q1)

      //      sum = QR_Full + QL_Full
      //      i = size - 2
      //      while (i > 0) {
      //        sum += h * Q1(i)
      //        i -= 1
      //      }
      //      println("t=" + t + " LawsonSwayne-1.5 Q=" + sum)

      var i = 0
      while (i < size) {
        Q1(i) = (sqrt2 + 1) * Q1(i) - sqrt2 * Q1Half(i)
        i += 1
      }
      QL_ = (sqrt2 + 1) * QL_Full - sqrt2 * QL_Half
      QR_ = (sqrt2 + 1) * QR_Full - sqrt2 * QR_Half
      advanceEm(dt*(1-2*b),Em_)
      //printSumQF(t)
      Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
      tIndex += 1
    }

  }

}

package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Epsilon

class HaganLawsonSwayneSABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  private val b = 1 - math.sqrt(2) / 2
  private val sqrt2 = math.sqrt(2)

  override def solve() {
    tri1 = new TridiagonalMatrix(size) //Full
    //    val tri1Half = new TridiagonalMatrix(size) //Half
      M1 = Array.ofDim[Double](size)
      M0=M1
    buildMcache(dt *b, dt*(1-2*b), M1)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    var Q1Half = Array.ofDim[Double](size)
    QL_ = 0.0
    var QL_Full, QL_Half, Q1HalfStart, Q1HalfEnd = 0.0
    var QR_Full, QR_Half = 0.0
    QR_ = 0.0
    var t = T
    //    var rhs = Array.ofDim[Double](size)
    var tIndex = 0
    var Qtmp = Q0_

    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
//             if (tIndex < 4) {
//              printTimeStep(t)
//              tIndex += 1
//            }
      t -= dt
       advance(dt*b, M1)
      computeSystem(tri1, dt * b)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      QL_Full = QL_ + dt * b / (h) * (M1(1) * Q1Half(1) - M1(0) * Q1Half(0))
      QR_Full = QR_ - dt * b / (h) * (M1(size - 1) * Q1Half(size - 1) - M1(size - 2) * Q1Half(size - 2))
      QL_Half = QL_Full
      QR_Half = QR_Full
      // //check sum
//                  var sum = QR_Half + QL_Half
//            var i = size - 2
//            while (i > 0) {
//              sum += h * Q1Half(i)
//              i -= 1
//            }
//            println("t=" + t + " LawsonSwayne-1 Q=" + sum+ "Q(f)="+Q1Half(j0)+" QL="+QL_Half+" QR="+QR_Half)

      //      t -= dt *b
      advance(dt*b, M1)
      computeSystem(tri1, dt * b)
      Q1Half(0) = 0
      Q1Half(size - 1) = 0
      solver.solve(tri1, Q1Half, Q1)
      QL_Full += dt * b / (h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
      QR_Full -= dt * b / (h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

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
      advance(dt*(1-2*b),M1)
      //      //check sum
      //      sum = QR_ + QL_
      //      i = size - 2
      //      while (i > 0) {
      //        sum += h * Q1(i)
      //        i -= 1
      //      }
      //      println("t=" + t + " LawsonSwayne-2 Q=" + sum)
//      M0 = M1
      Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
    }

  }

  private def computeSystem(tri1: TridiagonalMatrix, dt: Double) {

    var j = 1
    val frac = dt / (h * h)
    while (j < size - 1) {
      tri1.lower(j) = -frac * M1(j - 1)
      tri1.middle(j) = 1 + 2 * frac * M1(j)
      tri1.upper(j) = -frac * M1(j + 1)
      //      tri0.lower(j) = frac * M0(j - 1)
      //      tri0.middle(j) = 1 - 2 * frac * M0(j)
      //      tri0.upper(j) = frac * M0(j + 1)
      j += 1
    }
    tri1.upper(0) = M1(1)
    tri1.middle(0) = M1(0)
    tri1.lower(size - 1) = M1(size - 2)
    tri1.middle(size - 1) = M1(size - 1)
  }
}
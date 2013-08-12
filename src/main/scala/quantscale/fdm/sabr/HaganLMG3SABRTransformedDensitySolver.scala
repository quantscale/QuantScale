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
    //    M0 = Array.ofDim[Double](size)
    buildEmCache(dt / 3, 0)
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    val Q1Third = Array.ofDim[Double](size)
    val Q1ThirdTmp = Array.ofDim[Double](size)
    val Q1Half = Array.ofDim[Double](size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    var Qtmp = Q0_
    var tIndex = 0

    while (tIndex < timeSteps) {
      //                  if (tIndex < 4) {
      //              printTimeStep(t)
      //              tIndex += 1
      //            }
      t -= dt / 3

      advanceEm(dt / 3, Em_)
      computeSystem(dt / 3, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1ThirdTmp)
      var QL_ThirdPart =  QL_ + dt/3 *  Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1ThirdTmp(1)
      var QR_ThirdPart = QR_ + dt/3  * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1ThirdTmp(size - 2)
      var QL_HalfPart = QL_ThirdPart
      var QR_HalfPart = QR_ThirdPart
      val Q0_Init = Q0_
      Q0_ = Q1ThirdTmp

      t -= dt / 3
      advanceEm(dt / 3, Em_)
      computeSystem(dt/3, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Third)
      QL_ThirdPart += dt/3  *  Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1Third(1)
      QR_ThirdPart += dt/3  * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1Third(size - 2)
      Q0_ = Q1Third

      t -= dt / 3
      advanceEm(dt / 3, Em_)
      computeSystem(dt / 3, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Third)
      QL_ThirdPart += dt/3 *  Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1Third(1)
      QR_ThirdPart += dt/3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1Third(size - 2)

      Q0_ = Q1ThirdTmp
      computeSystem(2 * dt / 3, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1Half)
      QL_HalfPart += 2*dt/3 *  Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1Half(1)
      QR_HalfPart += 2*dt/3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1Half(size - 2)

      Q0_ = Q0_Init
      computeSystem(dt, Em_, tri1)
      Q0_(0) = 0
      Q0_(size - 1) = 0
      solver.solve(tri1, Q0_, Q1)

      QL_ += 4.5 * QL_ThirdPart - 4.5 * QL_HalfPart + dt *  Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1(1)
      QR_ += 4.5 * QR_ThirdPart - 4.5 * QR_HalfPart + dt * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1Half(size - 2)

      var i = 0
      while (i < size) {
        Q1(i) = 4.5 * Q1Third(i) - 4.5 * Q1Half(i) + Q1(i)
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
      Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
      tIndex +=1
    }
  }

  private def computeSystem(dt: Double, M1: Array[Double], tri1: TridiagonalMatrix) {

    var j = 1
    val frac = dt / (2*h)
    while (j < size - 1) {
      val a = Cm_(j - 1) / (Fm_(j) - Fm_(j - 1));
      val b = Cm_(j) * (1 / (Fm_(j + 1) - Fm_(j)) + 1 / (Fm_(j) - Fm_(j - 1)))
      val c = Cm_(j + 1) / (Fm_(j + 1) - Fm_(j));
      tri1.lower(j) = -a * frac * M1(j - 1)
      tri1.middle(j) = 1 + b * frac * M1(j)
      tri1.upper(j) = -c * frac * M1(j + 1)
      j += 1
    }
    tri1.upper(0) = Cm_(1) / (Fm_(1) - Fm_(0)) * M1(1)
    tri1.middle(0) = Cm_(0) / (Fm_(1) - Fm_(0)) * M1(0)
    tri1.lower(size - 1) = Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * M1(size - 2)
    tri1.middle(size - 1) = Cm_(size - 1) / (Fm_(size - 1) - Fm_(size - 2)) * M1(size - 1)
  }

}

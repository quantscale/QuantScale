package quantscale.fdm.sabr

import quantscale.fdm.TridiagonalMatrix
import quantscale.analytic.SABRModelSpec
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.Epsilon
import quantscale.analytic.SABRVanilla

class HaganDouglasSABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0)
  extends HaganSABRDensitySolver(spec, forward, T, size, timeSteps, FmaxTruncation) {

  var theta = 0.5

  override def solve() {

    //smooth payoff: take Tsmooth = T-0.01 or T-dt whichever is greater 

    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    M0 = Array.ofDim[Double](size)
    M1 = Array.ofDim[Double](size)
    if (useRannacher) {
      buildMcache(dt / 2, dt, M0)
    } else {
      buildMcache(dt, 0, M0)
    }
    theta = math.max(0,0.5 - h*h/(12*M0(j0)*dt))
    println("theta=" + theta + " " + computeCourantNumber + " " + (dt * M0(size-1) / (h * h)))
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    var indexRannacher = 0
    var rhs = Array.ofDim[Double](size)
    var tIndex = 0

    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      //                        if (tIndex < 4) {
      //                          printTimeStep(t)
      //                          tIndex += 1
      //                        }
      if (useRannacher && indexRannacher < 2) {
        val thetaTmp = theta
        theta = 1.0
        t -= dt / 2
        System.arraycopy(M0, 0, M1, 0, size)
        advance(dt / 2, M1)
        computeSystem(dt / 2)
        Q0_(0) = 0
        Q0_(size - 1) = 0
        solver.solve(tri1, Q0_, Q1)
        QL_ += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
        QR_ -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))
        val Qtmp = Q0_
        Q0_ = Q1
        Q1 = Qtmp
        t -= dt / 2
        advance(dt / 2, M1)
        computeSystem(dt / 2)
        Q0_(0) = 0
        Q0_(size - 1) = 0
        solver.solve(tri1, Q0_, Q1)
        QL_ += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
        QR_ -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))
        theta = thetaTmp
        indexRannacher += 1
      } else {
        t -= dt
        System.arraycopy(M0, 0, M1, 0, size)
        advance(dt, M1)
        computeSystem(dt)
        tri0.multiply(Q0_, rhs)
        solver.solve(tri1, rhs, Q1)
        QL_ += dt / (h) * (theta * (M1(1) * Q1(1) - M1(0) * Q1(0)) + (1 - theta) * (M0(1) * Q0_(1) - M0(0) * Q0_(0)))
        QR_ -= dt / (h) * (theta * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2)) + (1 - theta) * (M0(size - 1) * Q0_(size - 1) - M0(size - 2) * Q0_(size - 2)))

        //        QL_ += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0) + M0(1) * Q0_(1) - M0(0) * Q0_(0))
        //        QR_ -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2) + M0(size - 1) * Q0_(size - 1) - M0(size - 2) * Q0_(size - 2))
      }
      //            var sumQ = QR_ + QL_
      //            var sumF = Fmin * QL_ + Fmax * QR_
      //            var i = size - 2
      //            while (i > 0) {
      //              sumQ += h * Q1(i)
      //              sumF += h * Q1(i) * F(i)
      //              i -= 1
      //            }
      //            println("t=" + t + " CN Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0) + " QL_=" + QL_ + " QR_=" + QR_)

      val Mtmp = M0
      M0 = M1
      M1 = Mtmp
      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
    }
  }

  override private[sabr] def computeSystem(dt: Double) {

    var j = 1
    val frac = dt / (h * h)
    while (j < size - 1) {
      tri1.lower(j) = -frac * theta * M1(j - 1)
      tri1.middle(j) = 1 + 2 * theta * frac * M1(j)
      tri1.upper(j) = -frac * theta * M1(j + 1)
      tri0.lower(j) = frac * (1 - theta) * M0(j - 1)
      tri0.middle(j) = 1 - 2 * frac * (1 - theta) * M0(j)
      tri0.upper(j) = frac * (1 - theta) * M0(j + 1)
      j += 1
    }
//    tri1.upper(0) = 0
//    tri1.middle(0) = 1
//      tri1.lower(size - 1) = 0
//      tri1.middle(size - 1) = 1
//        tri0.upper(0) = 0
//    tri0.middle(0) = 1
//      tri0.lower(size - 1) = 0
//      tri0.middle(size - 1) = 1
    
    tri1.upper(0) = M1(1)
    tri1.middle(0) = M1(0)
    if (useApproxBoundary) {
      tri1.lower(size - 1) = 0
      tri1.middle(size - 1) = 1
    } else {
      tri1.lower(size - 1) = M1(size - 2)
      tri1.middle(size - 1) = M1(size - 1)
    }
  }

}
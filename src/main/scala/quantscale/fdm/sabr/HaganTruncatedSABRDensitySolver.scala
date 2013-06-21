package quantscale.fdm.sabr

import quantscale.fdm.Epsilon
import quantscale.fdm.TridiagonalMatrix
import quantscale.analytic.SABRModelSpec
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.analytic.SABRVanilla

class HaganTruncatedSABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5) {
  private[sabr] var M0, M1: Array[Double] = null
  private[sabr] var Q0, Q1: Array[Double] = null
  private[sabr] var Fmin = 0
  private[sabr] var Fmax = math.min(forward * 10, SABRVanilla.computeFmax(spec, forward, T))
  private[sabr] var h = 0.0
  private[sabr] var j0 = 0
  private[sabr] var QL, QR = 0.0
  private[sabr] var dt = T / timeSteps
  private[sabr] var F: Array[Double] = computeFAndh()
  private[sabr] var tri1, tri0: TridiagonalMatrix = null
  private[sabr] var gammaCache, factorCache: Array[Double] = null
  private[sabr] var solver = new ThomasTridiagonalSolver()
  solver.init(size)

  var useRannacher: Boolean = false

  def solve() {
    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    M0 = computeM(T)
    Q0 = computeQ()
    Q1 = Array.ofDim(size)
    QL = 0.0
    QR = 0.0
    var t = T
    var indexRannacher = 0
    var rhs = Array.ofDim[Double](size)
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      if (useRannacher && indexRannacher < 2) {
        t -= dt / 2
        M1 = computeM(t)
        computeSystem() //will use dt/2 because of Crank
        //                Q0(0) = 0
        //                Q0(size - 1) = 0
        solver.solve(tri1, Q0, Q1)
        //          QL+=dt/(4*h)*(-3*M1(0)*Q1(0)+4*M1(1)*Q1(1)-M1(2)*Q1(2))
        //        QR -= dt / (4 * h) * (3*M1(size - 1) * Q1(size - 1) -4* M1(size - 2) * Q1(size - 2)+M1(size - 3) * Q1(size - 3))

                QL += dt / (2*h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
                QR -= dt / (2*h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

        val Qtmp = Q0
        Q0 = Q1
        Q1 = Qtmp
        t -= dt / 2
        M1 = computeM(t)
        computeSystem() //will use dt/2 because of Crank
        //                 Q0(0) = 0
        //                Q0(size - 1) = 0
        solver.solve(tri1, Q0, Q1)
        //          QL+=dt/(4*h)*(-3*M1(0)*Q1(0)+4*M1(1)*Q1(1)-M1(2)*Q1(2))
        //        QR -= dt / (4 * h) * (3*M1(size - 1) * Q1(size - 1) -4* M1(size - 2) * Q1(size - 2)+M1(size - 3) * Q1(size - 3))

        QL += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
        QR -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))

        indexRannacher += 1
      } else {
        t -= dt
        M1 = computeM(t)
        computeSystem()
        tri0.multiply(Q0, rhs)
        //        rhs(0) = 0
        //        rhs(size-1) = 0
        solver.solve(tri1, rhs, Q1)
        //alpha3=1/(2h), alpha2 = -2/h, alpha1= 2/h-1/2h=3/2h
        //(3*Q1(size-1)-4*Q1(size-2)+1*Q1(size-3))/(2*h)
        //        QL+=dt/(4*h)*(-3*M1(0)*Q1(0)+4*M1(1)*Q1(1)-M1(2)*Q1(2)-3*M0(0)*Q0(0)+4*M0(1)*Q0(1)-M0(2)*Q0(2))
        //        QR -= dt / (4 * h) * (3*M1(size - 1) * Q1(size - 1) -4* M1(size - 2) * Q1(size - 2)+M1(size - 3) * Q1(size - 3) +
        //            3*M0(size - 1) * Q0(size - 1) - 4*M0(size - 2) * Q0(size - 2)+M0(size - 3) * Q0(size - 3))

                QL += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0) + M0(1) * Q0(1) - M0(0) * Q0(0))
                QR -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2) + M0(size - 1) * Q0(size - 1) - M0(size - 2) * Q0(size - 2))

//        QL += dt / (2 * h) * (2 * M1(1) * Q1(1) + 2 * M0(1) * Q0(1))
//        QR -= dt / (2 * h) * (2 * M1(size - 1) * Q1(size - 1) + 2 * M0(size - 1) * Q0(size - 1))

      }

      var sumQ = QR + QL
      var sumF = Fmin * QL + Fmax * QR
      var i = size - 1
      while (i >= 0) {
        sumQ += h * Q1(i)
        sumF += h * Q1(i) * F(i)
        i -= 1
      }
      println("t=" + t + " TruncCN Q=" + sumQ + " F=" + sumF + " Fj0=" + F(j0) + " QL=" + QL + " QR=" + QR)

      M0 = M1
      val Qtmp = Q0
      Q0 = Q1
      Q1 = Qtmp
    }
  }

  def price(isCall: Boolean, strike: Double): Double = {
    if (Q1 == null) {
      solve()
    }
    if (isCall) {
      if (strike < Fmin) {
        return forward - strike
      } else if (strike > Fmax) {
        return 0.0
      } else {
        var k = 0
        while (strike > Fmin + (k + 0.5) * h) {
          k += 1
        }
        var term = (Fmin + (k + 0.5) * h - strike)
        var price = 0.5 * term * term * Q1(k) + (Fmax - strike) * QR
        k += 1
        while (k < size - 1) {
          price += (Fmin + (k) * h - strike) * h * Q1(k)
          k += 1
        }
        return price
      }
    } else {
      if (strike < Fmin) {
        return 0
      } else if (strike > Fmax) {
        return strike - forward
      } else {
        var k = 0
        while (strike > Fmin + (k + 0.5) * h) {
          k += 1
        }
        var term = strike - Fmin - (k - 0.5) * h
        var price = 0.5 * term * term * Q1(k) + (strike - Fmin) * QL
        k -= 1
        while (k > 0) {
          price += (strike - Fmin - (k) * h) * h * Q1(k)
          k -= 1
        }
        return price
      }
    }
  }

  private[sabr] def computeSystem() {

    var j = 1
    val frac = dt / (2 * h * h)
    while (j < size - 1) {
      tri1.lower(j) = -frac * M1(j - 1)
      tri1.middle(j) = 1 + 2 * frac * M1(j)
      tri1.upper(j) = -frac * M1(j + 1)
      tri0.lower(j) = frac * M0(j - 1)
      tri0.middle(j) = 1 - 2 * frac * M0(j)
      tri0.upper(j) = frac * M0(j + 1)
      j += 1
    }
    tri1.upper(0) = 0
    tri1.middle(0) = 1
    tri1.lower(size - 1) = 0
    tri1.middle(size - 1) = 1
  }

  private[sabr] def computeQ(): Array[Double] = {
    var Q = Array.ofDim[Double](size)
    Q(j0) = 1.0 / h
    return Q
  }

  private[sabr] def computeM(t: Double): Array[Double] = {
    var M = Array.ofDim[Double](size)
    var j = 0
    val efactor = spec.rho * spec.nu * spec.alpha * t
    if (gammaCache == null) {
      gammaCache = Array.ofDim[Double](size)
      factorCache = Array.ofDim[Double](size)
      val C0 = math.pow(forward, spec.beta)
      val forwardonebeta = math.pow(forward, 1 - spec.beta)
      while (j < size) {
        val f = F(j)
        if (f <= Epsilon.MACHINE_EPSILON_SQRT) {
          gammaCache(j) = 0.0
          factorCache(j) = 0.0
        } else {
          val C = math.pow(f, spec.beta)
          val fonebeta = f / C
          val z = (fonebeta - forwardonebeta) / (spec.alpha * (1 - spec.beta))
          val factor = C * C * 0.5 * spec.alpha * spec.alpha * (1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z)
          val gamma = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) 1 / forwardonebeta else (C - C0) / (f - forward)
          gammaCache(j) = gamma
          factorCache(j) = factor
        }
        j += 1
      }
    }

    j = 0
    while (j < size) {
      M(j) = factorCache(j) * math.exp(efactor * gammaCache(j))
      j += 1
    }
    return M
  }
  private[sabr] def computeFAndh(): Array[Double] = {
    var j = 0
    var F = Array.ofDim[Double](size)
    var h0 = (Fmax - Fmin) / size
    //    if (h0 > Fmin) {
    //      Fmin = h0 * 2
    //      h0 = (Fmax - Fmin) / size
    //    }
    j0 = (((forward - Fmin) / h0)).toInt
    h = (forward - Fmin) / (j0)
    while (j < size) {
      F(j) = Fmin + (j) * h
      j += 1
    }
    Fmax = Fmin + (j - 1) * h //Fmin+Jh and Fj goes from 0 to J+1=> size-1 = J+1
    return F
  }

}
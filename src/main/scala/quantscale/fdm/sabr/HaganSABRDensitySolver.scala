package quantscale.fdm.sabr

import quantscale.fdm.TridiagonalMatrix
import quantscale.analytic.SABRModelSpec
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.Epsilon
import quantscale.analytic.SABRVanilla

class HaganSABRDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, FmaxTruncation: Double = 5.0) {
  private[sabr] var M0, M1: Array[Double] = null
  private[sabr] var Q0_, Q1: Array[Double] = null
  private[sabr] var Fmin = spec.b
  private[sabr] var Fmax = math.min(FmaxTruncation, SABRVanilla.computeFmax(spec, forward, T))
  private[sabr] var h_ = 0.0
  private[sabr] var j0 = 0
  private[sabr] var QL_, QR_ = 0.0
  private[sabr] var dt_ = T / timeSteps
  private[sabr] var F_ = computeFAndh()
  private[sabr] var tri1, tri0: TridiagonalMatrix = null
  private[sabr] var M1cache, M2cache, gammaCache, factorCache: Array[Double] = null
  private[sabr] var dt1cache, dt2cache = 0.0
  private[sabr] var solver = new ThomasTridiagonalSolver()
  solver.init(size)

  var useRannacher: Boolean = false
  var useSmoothing: Boolean = false
  var useApproxBoundary: Boolean = false
  private val forwardonebeta = math.pow(forward + spec.b, 1 - spec.beta)

  def h: Double = h_

  def Q0: Array[Double] = Q0_
  def QL = QL_
  def QR = QR_

  def dt = dt_

  def indexForward = j0

  def F: Array[Double] = F_

  def printTimeStep(t: Double) {
    var str = ""
    for (i <- 0 until Q0.length) {
      str += t + " " + F(i) + " " + Q0(i) + " " + M0(i) + " " + M1(i) + "\n"
    }
    println(str)
  }

  private def computeApproxQ(j: Int, dt: Double): Double = {
    val f = math.abs(F(j) + spec.b) - spec.b
    val C = math.pow(f + spec.b, spec.beta)
    val fonebeta = f / C
    val z = (fonebeta - forwardonebeta) / (spec.alpha * (1 - spec.beta))
    val y = 1 / spec.nu * math.log((math.sqrt(1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z) + spec.rho + spec.nu * z) / (1 + spec.rho))
    val D = (1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z) * math.exp(spec.rho * spec.nu * spec.alpha * gammaCache(j) * dt) * C * C
    val dy = 1 / (spec.alpha * math.sqrt(D))
    val Q = 1 / math.sqrt(2 * math.Pi * dt) * math.exp(-y * y / (2 * dt)) * dy
    return Q
  }

  def solve() {

    //smooth payoff: take Tsmooth = T-0.01 or T-dt whichever is greater 

    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    M0 = Array.ofDim[Double](size)
    M1 = Array.ofDim[Double](size)
    if (useRannacher || useSmoothing) {
      buildMcache(dt / 2, dt, M0)
    } else {
      buildMcache(dt, 0, M0)
    }
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
      if (useSmoothing && indexRannacher < 1) {
        val dtSmooth = dt / 2
        t -= dtSmooth
        //compute smooth approx
        var j = 0
        while (j < Q0_.length) {
          Q0_(j) = computeApproxQ(j, dtSmooth)
          j += 1
        }
        Array.copy(M0, 0, M1, 0, size)
        advance(dtSmooth, M1)
        QL_ += dtSmooth / (h) * (M1(1) * Q0_(1) - M1(0) * Q0_(0))
        QR_ -= dtSmooth / (h) * (M1(size - 1) * Q0_(size - 1))

        //        printTimeStep(t)
        M0 = M1
        val dtRem = dt - dtSmooth
        t -= dtRem
        //compute CN with step dt-dtSmooth
        advance(dtSmooth, M1)
        computeSystem(dtRem)
        Q0_(0) = 0
        Q0_(size - 1) = 0
        tri0.multiply(Q0_, rhs)
        solver.solve(tri1, rhs, Q1)
        QL_ += dtRem / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0) + M0(1) * Q0_(1) - M0(0) * Q0_(0))
        QR_ -= dtRem / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2) + M0(size - 1) * Q0_(size - 1) - M0(size - 2) * Q0_(size - 2))
        indexRannacher += 1
      } else if (useRannacher && indexRannacher < 2) {
        t -= dt / 2
        Array.copy(M0, 0, M1, 0, size)
        advance(dt / 2, M1)
        computeSystem(dt) //will use dt/2 because of Crank
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
        computeSystem(dt) //will use dt/2 because of Crank
        Q0_(0) = 0
        Q0_(size - 1) = 0
        solver.solve(tri1, Q0_, Q1)
        QL_ += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0))
        QR_ -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2))
        indexRannacher += 1
      } else {
        t -= dt
        Array.copy(M0, 0, M1, 0, size)
        advance(dt, M1)
        computeSystem(dt)
        tri0.multiply(Q0_, rhs)
        if (useApproxBoundary) {
          rhs(size - 1) = computeApproxQ(size - 1, T - t)
        }
        solver.solve(tri1, rhs, Q1)
        //         if (useApproxBoundary) {
        //          Q1(size - 1) = computeApproxQ(size - 1, T - t)
        //        }
        QL_ += dt / (2 * h) * (M1(1) * Q1(1) - M1(0) * Q1(0) + M0(1) * Q0_(1) - M0(0) * Q0_(0))
        QR_ -= dt / (2 * h) * (M1(size - 1) * Q1(size - 1) - M1(size - 2) * Q1(size - 2) + M0(size - 1) * Q0_(size - 1) - M0(size - 2) * Q0_(size - 2))

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

  def price(isCall: Boolean, strike: Double): Double = {
    if (Q0_ == null) {
      solve()
    }
    if (isCall) {
      if (strike < Fmin) {
        return forward - strike
      } else if (strike > Fmax) {
        return 0.0
      } else {
        var k = 0
        while (strike > Fmin + (k) * h) {
          k += 1
        }
        val term = (Fmin + k * h - strike)
        var price = 0.5 * term * term * Q0_(k) + (Fmax - strike) * QR_
        k += 1
        while (k < size - 1) {
          price += (Fmin + (k - 0.5) * h - strike) * h * Q0_(k)
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
        while (strike > Fmin + (k) * h) {
          k += 1
        }
        val term = strike - Fmin - (k - 1) * h
        var price = 0.5 * term * term * Q0_(k) + (strike - Fmin) * QL_
        k -= 1
        while (k > 0) {
          price += (strike - Fmin - (k - 0.5) * h) * h * Q0_(k)
          k -= 1
        }
        return price
      }
    }
  }

  private[sabr] def computeSystem(dt: Double) {

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

  private[sabr] def computeQ(): Array[Double] = {
    var Q = Array.ofDim[Double](size)
    Q(j0) = 1.0 / h
    return Q
  }

  //  private def computeM(t: Double): Array[Double] = {
  //    var M = Array.ofDim[Double](size)
  //    var j = 0
  //    val efactor = spec.rho * spec.nu * spec.alpha * t
  //    if (gammaCache == null) {
  //      gammaCache = Array.ofDim[Double](size)
  //      factorCache = Array.ofDim[Double](size)
  //      val C0 = math.pow(forward, spec.beta)
  //      while (j < size) {
  //        var f = math.abs(F(j)) //abs necessary for Q0
  //        if (f <= Epsilon.MACHINE_EPSILON_SQRT) {
  //          gammaCache(j) = 0.0
  //          factorCache(j) = 0.0
  //        } else {
  //          val C = math.pow(f, spec.beta)
  //          val fonebeta = f / C
  //          val z = (fonebeta - forwardonebeta) / (spec.alpha * (1 - spec.beta))
  //          val factor = C * C * 0.5 * spec.alpha * spec.alpha * (1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z)
  //          val gamma = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) spec.beta / forwardonebeta else (C - C0) / (f - forward)
  //          gammaCache(j) = gamma
  //          factorCache(j) = factor
  //        }
  //        j += 1
  //
  //      }
  //    }
  //
  //    j = 0
  //    while (j < size) {
  //      M(j) = factorCache(j) * math.exp(efactor * gammaCache(j))
  //      j += 1
  //    }
  //    return M
  //  }

  private[sabr] def buildMcache(dt1: Double, dt2: Double, M: Array[Double]) {
    M1cache = Array.ofDim[Double](size)
    if (dt2 != 0) M2cache = Array.ofDim[Double](size)
    dt1cache = dt1
    dt2cache = dt2
    var j = 0
    val efactor = spec.rho * spec.nu * spec.alpha
    val C0 = math.pow(forward + spec.b, spec.beta)
    while (j < size) {
      var gamma = 0.0
      var factor = 0.0
      var f = math.abs(F(j) + spec.b) - spec.b //abs necessary for Q0
      if (f <= Epsilon.MACHINE_EPSILON_SQRT) {
        gamma = 0.0
        factor = 0.0
      } else {
        val C = math.pow(f + spec.b, spec.beta)
        val fonebeta = f / C
        val z = (fonebeta - forwardonebeta) / (spec.alpha * (1 - spec.beta))
        factor = C * C * 0.5 * spec.alpha * spec.alpha * (1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z)
        gamma = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) spec.beta / forwardonebeta else (C - C0) / (f - forward)

      }
      M1cache(j) = math.exp(efactor * gamma * dt1)
      if (dt2 != 0) M2cache(j) = math.exp(efactor * gamma * dt2)
      M(j) = factor
      j += 1
    }
  }

  private[sabr] def advance(dt: Double, M: Array[Double]) {
    var j = 0
    val Mcache = if (dt == dt1cache) M1cache else M2cache
    while (j < size) {
      M(j) = M(j) * Mcache(j)
      j += 1
    }
  }

  def computeCourantNumber(): Double = {
    var cfl = spec.alpha * math.pow(forward, spec.beta)
    cfl = cfl * cfl //* math.exp(spec.rho * spec.nu * spec.alpha * T * spec.beta * math.pow(forward, spec.beta - 1))
    cfl *= 0.5*dt / (h * h)
    return cfl
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
    h_ = (forward - Fmin) / (j0 - 0.5)
    while (j < size) {
      F(j) = Fmin + (j - 0.5) * h
      j += 1
    }
    Fmax = Fmin + (j - 2) * h //Fmin+Jh and Fj goes from 0 to J+1=> size-1 = J+1

    //    println("Fmax = "+Fmax+" h="+h+" CFL="+computeCourantNumber)
    return F
  }

}
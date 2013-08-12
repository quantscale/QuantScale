package quantscale.fdm.sabr

import quantscale.analytic.{SABRModelSpec}
import quantscale.fdm.{Epsilon, ThomasTridiagonalSolver, TridiagonalMatrix}
import scala.Array

/**
 * Date: 8/12/13
 */
class HaganSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, nDeviations: Double = 3.0) extends SABRDensitySolver {
  private[sabr] var Q0_, Q1: Array[Double] = null
  private[sabr] var Fmin = spec.b
  private[sabr] var Fmax = 0.0
  private[sabr] var zmin, zmax = 0.0
  private[sabr] var h_ = 0.0
  private[sabr] var j0 = 0
  private[sabr] var QL_, QR_ = 0.0
  private[sabr] var dt_ = T / timeSteps
  private[sabr] var Cm_, Fm_, Gammam_, Em_, Em1Cache, Em2Cache: Array[Double] = null
  private[sabr] var tri1, tri0: TridiagonalMatrix = null
  private[sabr] var dt1cache, dt2cache = 0.0
  private[sabr] var solver = new ThomasTridiagonalSolver()
  solver.init(size)

  var useRannacher: Boolean = false
  var useApproxBoundary: Boolean = false
  private val forwardonebeta = math.pow(forward + spec.b, 1 - spec.beta)

  initGrid()

  def h: Double = h_

  def Q0: Array[Double] = Q0_

  def QL = QL_

  def QR = QR_

  def dt = dt_

  def indexForward = j0

  def Fm: Array[Double] = Fm_

  def printTimeStep(t: Double) {
    var str = ""
    for (i <- 0 until Q0.length) {
      str += t + " " + Fm(i) + " " + Q0(i) + " " + Em_(i) + "\n"
    }
    println(str)
  }

  def solve() {

    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val M0 = Array.ofDim[Double](size)
    if (useRannacher) {
      buildEmCache(dt / 2, dt)
    } else {
      buildEmCache(dt, 0)
    }
    Q0_ = computeQ()
    Q1 = Array.ofDim(size)
    QL_ = 0.0
    QR_ = 0.0
    var t = T
    var indexRannacher = 0
    var rhs = Array.ofDim[Double](size)
    var tIndex = 0

    while (tIndex < timeSteps) {
      //                        if (tIndex < 4) {
      //                          printTimeStep(t)
      //                          tIndex += 1
      //                        }
      if (useRannacher && indexRannacher < 2) {
        t -= dt / 2
        advanceEm(dt / 2, Em_)
        computeSystem(dt, Em_, null) //will use dt/2 because of Crank
        Q0_(0) = 0
        Q0_(size - 1) = 0
        solver.solve(tri1, Q0_, Q1)
        QL_ += dt/2 * Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1(1)
        QR_ += dt/2 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1(size - 2)
        val Qtmp = Q0_
        Q0_ = Q1
        Q1 = Qtmp
        t -= dt / 2
        advanceEm(dt / 2, Em_)
        computeSystem(dt, Em_, null) //will use dt/2 because of Crank
        Q0_(0) = 0
        Q0_(size - 1) = 0
        solver.solve(tri1, Q0_, Q1)
        QL_ += dt/2 * Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * Q1(1)
        QR_ += dt/2 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * Q1(size - 2)
        indexRannacher += 1
      } else {
        t -= dt
        Array.copy(Em_, 0, M0, 0, size)
        advanceEm(dt, Em_)
        computeSystem(dt, Em_, M0)
        Q0_(0) = 0
        Q0_(size - 1) = 0
        tri0.multiply(Q0_, rhs)
        solver.solve(tri1, rhs, Q1)

        QL_ += dt/2  * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * Q1(1) + M0(1) * Q0(1))
        QR_ += dt/2  * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * Q1(size - 2) + M0(size - 2) * Q0(size - 2))
      }
     // printSumQF(t)

      val Qtmp = Q0_
      Q0_ = Q1
      Q1 = Qtmp
      tIndex += 1
    }
  }


  def printSumQF(t: Double) {
    var sumQ = QR_ + QL_
    var sumF = Fmin * QL_ + Fmax * QR_
    var i = size - 2
    while (i > 0) {
      sumQ += h * Q1(i)
      sumF += h * Q1(i) * Fm_(i)
      i -= 1
    }
    println("t=" + t + " CN Q=" + sumQ + " F=" + sumF + " Fj0=" + Fm_(j0) + " QL_=" + QL_ + " QR_=" + QR_)
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
        val zstrike = computeZFromY(computeYFromF(strike))
        while (zstrike > zmin + (k) * h) {
          k += 1
        }
        val ztilde = zmin + k * h
        val ftilde = computeFFromY(computeYFromZ(ztilde))
        val term = (ftilde - strike)
        var price = 0.5 * term * term * Q0_(k) + (Fmax - strike) * QR_
        k += 1
        while (k < size - 1) {
          price += (Fm_(k) - strike) * h * Q0_(k)
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
        val zstrike = computeZFromY(computeYFromF(strike))
        while (zstrike > zmin + (k) * h) {
          k += 1
        }
        val ztilde = zmin + k * h
        val ftilde = computeFFromY(computeYFromZ(ztilde))
        val term = strike - ftilde
        var price = 0.5 * term * term * Q0_(k) + (strike - Fmin) * QL_
        k -= 1
        while (k > 0) {
          price += (strike - Fm_(k)) * h * Q0_(k)
          k -= 1
        }
        return price
      }
    }
  }

  private[sabr] def computeSystem(dt: Double, M1: Array[Double], M0: Array[Double]) {

    var j = 1
    val frac = dt / (4 * h)
    while (j < size - 1) {
      val a = Cm_(j - 1) / (Fm_(j) - Fm_(j - 1));
      val b = Cm_(j) * (1 / (Fm_(j + 1) - Fm_(j)) + 1 / (Fm_(j) - Fm_(j - 1)))
      val c = Cm_(j + 1) / (Fm_(j + 1) - Fm_(j));
      tri1.lower(j) = -a * frac * M1(j - 1)
      tri1.middle(j) = 1 + b * frac * M1(j)
      tri1.upper(j) = -c * frac * M1(j + 1)
      if (M0 != null) {
        tri0.lower(j) = a * frac * M0(j - 1)
        tri0.middle(j) = 1 - b * frac * M0(j)
        tri0.upper(j) = c * frac * M0(j + 1)
      }
      j += 1
    }
    tri1.upper(0) = Cm_(1) / (Fm_(1) - Fm_(0)) * M1(1)
    tri1.middle(0) = Cm_(0) / (Fm_(1) - Fm_(0)) * M1(0)
    tri1.lower(size - 1) = Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * M1(size - 2)
    tri1.middle(size - 1) = Cm_(size - 1) / (Fm_(size - 1) - Fm_(size - 2)) * M1(size - 1)
    if (M0 != null) {
      tri0.upper(0) = -Cm_(1) / (Fm_(1) - Fm_(0)) * M0(1)
      tri0.middle(0) = -Cm_(0) / (Fm_(1) - Fm_(0)) * M0(0)
      tri0.lower(size - 1) = -Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * M0(size - 2)
      tri0.middle(size - 1) = -Cm_(size - 1) / (Fm_(size - 1) - Fm_(size - 2)) * M0(size - 1)
    }
  }

  private[sabr] def computeQ(): Array[Double] = {
    val Q = Array.ofDim[Double](size)
    Q(j0) = 1.0 / h
    return Q
  }


  private[sabr] def buildEmCache(dt1: Double, dt2: Double) {
    Em_ = Array.ofDim[Double](size)
    Em1Cache = Array.ofDim[Double](size)
    if (dt2 != 0) Em2Cache = Array.ofDim[Double](size)
    dt1cache = dt1
    dt2cache = dt2
    var j = size - 2
    val efactor = spec.rho * spec.nu * spec.alpha
    while (j > 0) {
      Em1Cache(j) = math.exp(efactor * Gammam_(j) * dt1)
      if (dt2 != 0) Em2Cache(j) = math.exp(efactor * Gammam_(j) * dt2)
      Em_(j) = 1.0
      j -= 1
    }
    Em_(0) = Em_(1)
    Em_(size-1) = Em_(size-2)
    Em1Cache(0) = Em1Cache(1)
    Em1Cache(size-1) = Em1Cache(size-2)
    if (dt2 != 0) {
      Em2Cache(0) = Em2Cache(1)
      Em2Cache(size-1) = Em2Cache(size-2)
    }
  }

  private[sabr] def advanceEm(dt: Double, M: Array[Double]) {
    var j = 0
    val Mcache = if (dt == dt1cache) Em1Cache else Em2Cache
    while (j < size) {
      M(j) = M(j) * Mcache(j)
      j += 1
    }
  }

  def computeCourantNumber(): Double = {
    val cfl = Cm_(j0) * (1 / (Fm_(j0 + 1) - Fm_(j0)) + 1 / (Fm_(j0) - Fm_(j0 - 1)))*dt/(2*h)
    return cfl
  }

  private[sabr] final def computeYFromF(f: Double): Double = {
    return (math.pow(f, 1 - spec.beta) - forwardonebeta) / (1 - spec.beta)
  }

  private[sabr] final def computeZFromY(y: Double): Double = {
    val temp = (spec.rho + spec.nu * y / spec.alpha)
    return -1 / spec.nu * math.log((math.sqrt(1 - spec.rho * spec.rho + temp * temp) - spec.rho - spec.nu * y / spec.alpha) / (1 - spec.rho))
  }

  private[sabr] final def computeYFromZ(z: Double): Double = {
//    val sh = math.sinh(spec.nu * z)
//    val ch = math.cosh(spec.nu * z)
    val e = math.exp(spec.nu*z)
    val sh = 0.5*(e-1/e)
    val ch = 0.5*(e+1/e)
    return spec.alpha / spec.nu * (sh + spec.rho * (ch - 1))
  }

  private[sabr] final def computeFFromY(y: Double): Double = {
    val factor = (1 - spec.beta) * y
    return if (factor + forwardonebeta < 0) 0.0 else math.pow(forwardonebeta + factor, 1 / (1 - spec.beta))
  }

  def initGrid() {
    val sqrtT = math.sqrt(T)
    zmin = -nDeviations * sqrtT
    zmax = -zmin
    if (spec.beta < 1) {
      val ybar = -forwardonebeta / (1 - spec.beta)
      val zbar = computeZFromY(ybar)
      if (zbar > zmin) {
        zmin = zbar
      }
    }
    var j = 0
    var zm = Array.ofDim[Double](size)
    var h0 = (zmax - zmin) / size
    j0 = (((0.0 - zmin) / h0)).toInt
    h_ = (0.0 - zmin) / (j0 - 0.5)
    while (j < size) {
      zm(j) = zmin + (j - 0.5) * h
      j += 1
    }
    zmax = zm(size - 1) + 0.5 * h
    j = size - 2
    Fm_ = Array.ofDim[Double](size)
    Cm_ = Array.ofDim[Double](size)
    Gammam_ = Array.ofDim[Double](size)
    val forwardbeta = math.pow(forward, spec.beta)
    while (j > 0) {
      val ym = computeYFromZ(zm(j))
      Fm_(j) = computeFFromY(ym)
      val f = Fm_(j)
      val fpow = math.pow(f, spec.beta)
      Cm_(j) = math.sqrt(spec.alpha * spec.alpha + 2 * spec.rho * spec.alpha * spec.nu * ym + spec.nu * spec.nu * ym * ym) * fpow
      Gammam_(j) = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) spec.beta / forwardonebeta else (fpow - forwardbeta) / (f - forward)
      j -= 1
    }

    Fmax = computeFFromY(computeYFromZ(zmax))
    Fmin = computeFFromY(computeYFromZ(zmin))
    Fm_(0) = 2 * Fmin - Fm_(1)
    Fm_(size - 1) = 2 * Fmax - Fm_(size - 2)
    Cm_(0) = Cm_(1)
    Cm_(size - 1) = Cm_(size - 2)
    //    println("Fmax = "+Fmax+" h="+h+" CFL="+computeCourantNumber)
  }

}
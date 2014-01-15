package quantscale.fdm.sabr

import quantscale.analytic.{SABRModelSpec}
import quantscale.fdm.{Epsilon, ThomasTridiagonalSolver, TridiagonalMatrix}
import scala.Array

/**
 * Date: 8/12/13
 *
 * TODO: rename Q in P, and the Q = P*h/(Fm+1-Fm)
 */

abstract class Smoothing {}

case class RannacherSmoothing extends Smoothing {}

case class TRBDF2BDF3Smoothing extends Smoothing {}

case class LMG3Smoothing extends Smoothing {}

case class BDF2Smoothing extends Smoothing {}

class HaganSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 1000, timeSteps: Int = 1000 / 5, nDeviations: Double = 3.0) extends SABRDensitySolver {
  private[sabr] var P0_, P1_ : Array[Double] = null
  private[sabr] var Fmin = spec.b
  private[sabr] var Fmax = 0.0
  private[sabr] var zmin_, zmax_ = 0.0
  private[sabr] var h_ = 0.0
  private[sabr] var j0 = 0
  private[sabr] var PL_, PR_ = 0.0
  private[sabr] var dt_ = T / timeSteps
  private[sabr] var Cm_, Fm_, Gammam_, Em_, Em1Cache, Em2Cache: Array[Double] = null
  private[sabr] var tri1, tri0: TridiagonalMatrix = null
  private[sabr] var dt1cache, dt2cache = 0.0
  private[sabr] var solver = new ThomasTridiagonalSolver()
  solver.init(size)

  var smoothing: Smoothing = null
  private val forwardonebeta = math.pow(forward + spec.b, 1 - spec.beta)

  initGrid()

  def name: String = {
    smoothing match {
      case RannacherSmoothing() => return "RAN"
      case TRBDF2BDF3Smoothing() => return "CN-SBDF"
      case LMG3Smoothing() => "LEF"
      case _ => return "CN"
    }
  }

  def h: Double = h_

  def P: Array[Double] = P0_

  def Q(i: Integer) = P0_(i) / (Fm_(i + 1) - Fm_(i)) * h

  def PL = PL_

  def PR = PR_

  def dt = dt_

  def indexForward = j0

  def Fm: Array[Double] = Fm_

  def zmin = zmin_

  def zmax = zmax_

  def printTimeStep(t: Double) {
    var str = ""
    for (i <- 1 until P.length - 1) {
      str += t + " " + Fm(i) + " " + P(i) + " " + Em_(i) + "\n"
    }
    println(str)
  }

  private[sabr] def computedPLdt(Em: Array[Double], P: Array[Double]): Double = {
    return Cm_(1) / (Fm_(1) - Fm_(0)) * Em(1) * P(1)
  }

  private[sabr] def computedPRdt(Em: Array[Double], P: Array[Double]): Double = {
    return Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em(size - 2) * P(size - 2)
  }

  def solve() {

    tri1 = new TridiagonalMatrix(size)
    tri0 = new TridiagonalMatrix(size)
    val M0 = Array.ofDim[Double](size)
    smoothing match {
      case RannacherSmoothing() | BDF2Smoothing() => buildEmCache(dt / 2, dt)
      case TRBDF2BDF3Smoothing() | LMG3Smoothing() => buildEmCache(dt / 3, dt)
      case _ => buildEmCache(dt, 0)
    }

    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    var indexRannacher = 0
    val rhs = Array.ofDim[Double](size)
    var tIndex = 0

    while (tIndex < timeSteps) {
      //                              if (tIndex < 4) {
      //                                printTimeStep(t)
      //                              }

      if (smoothing.isInstanceOf[RannacherSmoothing] && indexRannacher < 2) {
        t -= dt / 2
        advanceEm(dt / 2, Em_)
        computeSystem(dt, Em_, null, tri1, tri0) //will use dt/2 because of Crank
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, P1_)
        PL_ += dt / 2 * Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * P1_(1)
        PR_ += dt / 2 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * P1_(size - 2)
        val Qtmp = P0_
        P0_ = P1_
        P1_ = Qtmp
        t -= dt / 2
        advanceEm(dt / 2, Em_)
        computeSystem(dt, Em_, null, tri1, tri0) //will use dt/2 because of Crank
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, P1_)
        PL_ += dt / 2 * computedPLdt(Em_, P1_)
        PR_ += dt / 2 * computedPRdt(Em_, P1_)
        indexRannacher += 1
      } else if (smoothing.isInstanceOf[BDF2Smoothing] && indexRannacher > 2) {
        t -= dt / 2
        advanceEm(dt / 2, Em_)
        computeSystem(dt, Em_, null, tri1, tri0) //will use dt/2 because of Crank
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, P1_)
        val PL1 = PL_
        val PR1 = PR_
        PL_ += dt / 2 * Cm_(1) / (Fm_(1) - Fm_(0)) * Em_(1) * P1_(1)
        PR_ += dt / 2 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * Em_(size - 2) * P1_(size - 2)
        val Qtmp = P0_
        P0_ = P1_
        P1_ = Qtmp
        t -= dt / 2
        advanceEm(dt / 2, Em_)
        var j = size - 2
        while (j >= 1) {
          rhs(j) = (4 * P1_(j) - P0_(j)) / 3.0
          j -= 1
        }

        computeSystem(dt / 3.0, Em_, tri1)
        solver.solve(tri1, rhs, P1_)
        PL_ = (4 * PL_ - PL1) / 3.0 + 1.0 / 3.0 * dt * computedPLdt(Em_, P1_)
        PR_ = (4 * PR_ - PR1) / 3.0 + 1.0 / 3.0 * dt * computedPRdt(Em_, P1_)

        indexRannacher += 1

      } else if (smoothing.isInstanceOf[LMG3Smoothing] && indexRannacher < 2) {
        val Q1ThirdTmp = Array.ofDim[Double](size)
        val Q1Third = Array.ofDim[Double](size)
        val Q1Half = Array.ofDim[Double](size)
        t -= dt / 3

        advanceEm(dt / 3, Em_)
        computeSystem(dt / 3, Em_, tri1)
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, Q1ThirdTmp)
        var QL_ThirdPart = PL_ + dt / 3 * computedPLdt(Em_, Q1ThirdTmp)
        var QR_ThirdPart = PR_ + dt / 3 * computedPRdt(Em_, Q1ThirdTmp)
        var QL_HalfPart = QL_ThirdPart
        var QR_HalfPart = QR_ThirdPart
        val Q0_Init = P0_
        P0_ = Q1ThirdTmp

        t -= dt / 3
        advanceEm(dt / 3, Em_)
        computeSystem(dt / 3, Em_, tri1)
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, Q1Third)
        QL_ThirdPart += dt / 3 * computedPLdt(Em_, Q1Third)
        QR_ThirdPart += dt / 3 * computedPRdt(Em_, Q1Third)
        P0_ = Q1Third

        t -= dt / 3
        advanceEm(dt / 3, Em_)
        computeSystem(dt / 3, Em_, tri1)
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, Q1Third)
        QL_ThirdPart += dt / 3 * computedPLdt(Em_, Q1Third)
        QR_ThirdPart += dt / 3 * computedPRdt(Em_, Q1Third)

        P0_ = Q1ThirdTmp
        computeSystem(2 * dt / 3, Em_, tri1)
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, Q1Half)
        QL_HalfPart += 2 * dt / 3 * computedPLdt(Em_, Q1Half)
        QR_HalfPart += 2 * dt / 3 * computedPRdt(Em_, Q1Half)

        P0_ = Q0_Init
        computeSystem(dt, Em_, tri1)
        P0_(0) = 0
        P0_(size - 1) = 0
        solver.solve(tri1, P0_, P1_)

        PL_ += 4.5 * QL_ThirdPart - 4.5 * QL_HalfPart + dt * computedPLdt(Em_, P1_)
        PR_ += 4.5 * QR_ThirdPart - 4.5 * QR_HalfPart + dt * computedPRdt(Em_, P1_)

        var i = 0
        while (i < size) {
          P1_(i) = 4.5 * Q1Third(i) - 4.5 * Q1Half(i) + P1_(i)
          i += 1
        }
        indexRannacher += 1
      } else if (smoothing.isInstanceOf[TRBDF2BDF3Smoothing] && indexRannacher < 2) {
        val dt3 = dt / 3.0
        val Q2 = Array.ofDim[Double](size)
        t -= dt
        Array.copy(Em_, 0, M0, 0, size)
        advanceEm(dt3, Em_)
        computeSystem(dt3, Em_, M0, tri1, tri0)
        tri0.multiply(P0_, rhs)
        solver.solve(tri1, rhs, P1_)
        val QL1 = PL_
        val QR1 = PR_
        PL_ += 0.5 * dt3 * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * P1_(1) + M0(1) * P0_(1))
        PR_ += 0.5 * dt3 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * P1_(size - 2) + M0(size - 2) * P0_(size - 2))
        //printSumPF("TR/3", t, P1_, PL_, PR_)

        var j = size - 2
        while (j >= 1) {
          rhs(j) = (4 * P1_(j) - P0_(j)) / 3.0
          j -= 1
        }
        advanceEm(dt3, Em_)
        computeSystem(dt3 / 3.0 * 2, Em_, tri1)
        solver.solve(tri1, rhs, Q2)
        val QL2 = PL_
        val QR2 = PR_
        PL_ = (4 * PL_ - QL1) / 3.0 + 2.0 / 3.0 * dt3 * computedPLdt(Em_, Q2)
        PR_ = (4 * PR_ - QR1) / 3.0 + 2.0 / 3.0 * dt3 * computedPRdt(Em_, Q2)

        //printSumPF("TRBDF2/3", t, Q2, PL_, PR_)

        j = size - 2
        while (j >= 1) {
          rhs(j) = (18 * Q2(j) - 9 * P1_(j) + 2 * P0_(j)) / 11
          j -= 1
        }
        advanceEm(dt3, Em_)
        computeSystem(dt3 * 6.0 / 11, Em_, tri1)
        solver.solve(tri1, rhs, P1_)

        PL_ = (18 * PL_ - 9 * QL2 + 2 * QL1) / 11 + 6 * dt3 / 11 * computedPLdt(Em_, P1_)
        PR_ = (18 * PR_ - 9 * QR2 + 2 * QR1) / 11 + 6 * dt3 / 11 * computedPRdt(Em_, P1_)
        indexRannacher += 1
      } else {
        t -= dt
        Array.copy(Em_, 0, M0, 0, size)
        advanceEm(dt, Em_)
        computeSystem(dt, Em_, M0, tri1, tri0)
        P0_(0) = 0
        P0_(size - 1) = 0
        tri0.multiply(P0_, rhs)
        solver.solve(tri1, rhs, P1_)

        PL_ += dt / 2 * Cm_(1) / (Fm_(1) - Fm_(0)) * (Em_(1) * P1_(1) + M0(1) * P(1))
        PR_ += dt / 2 * Cm_(size - 2) / (Fm_(size - 1) - Fm_(size - 2)) * (Em_(size - 2) * P1_(size - 2) + M0(size - 2) * P(size - 2))
      }
      // printSumQF(t)

      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }
  }


  def printSumPF(name: String, t: Double, Q1: Array[Double], PL: Double, PR: Double) {
    var sumQ = PR + PL
    var sumF = Fmin * PL + Fmax * PR
    var i = size - 2
    while (i > 0) {
      sumQ += h * Q1(i)
      sumF += h * Q1(i) * Fm_(i)
      i -= 1
    }
    println("t=" + t + " " + name + " P=" + sumQ + " F=" + sumF + " Fj0=" + Fm_(j0) + " PL_=" + PL + " PR_=" + PR)
  }

  def price(isCall: Boolean, strike: Double): Double = {
    if (P0_ == null) {
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
        while (zstrike > zmin_ + (k) * h) {
          k += 1
        }
        val ztilde = zmin_ + k * h
        val ftilde = computeFFromY(computeYFromZ(ztilde))
        val term = (ftilde - strike)
        val dFdZ = if (math.abs(term) > 1e-5) (ftilde - Fm_(k)) / (h * 0.5) else 1

        //val dFdZ =  if (math.abs( term) > 1e-5) (ftilde - strike)/(ztilde-zstrike) else 1
        // val dFdZ = if (math.abs( term) > 1e-6) (Cm_(k)+Cm_(k-1))/2  else 1  //create discontinuity
        var price = 0.5 * term * term * P0_(k) / dFdZ + (Fmax - strike) * PR_
        k += 1
        while (k < size - 1) {
          price += (Fm_(k) - strike) * h * P0_(k)
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
        while (zstrike > zmin_ + (k) * h) {
          k += 1
        }
        val ztilde = zmin_ + (k) * h
        val ftilde = computeFFromY(computeYFromZ(ztilde))
        val term = strike - ftilde
        //  val dFdZ = if (math.abs(term) > 1e-6) (Cm_(k)+Cm_(k-1))/2 else 1   //create discontinuity
        //val dFdZ = if (math.abs( term) > 1e-5) (ftilde - strike)/(ztilde-zstrike) else 1

        val dFdZ = if (math.abs(term) > 1e-5) (ftilde - Fm_(k)) / (0.5 * h) else 1
        var price = 0.5 * term * term * P0_(k) / dFdZ + (strike - Fmin) * PL_
        while (k > 0) {
          price += (strike - Fm_(k)) * h * P0_(k)
          k -= 1
        }
        return price
      }
    }
  }

  private[sabr] def computeSystem(dt: Double, M1: Array[Double], tri1: TridiagonalMatrix) {

    var j = 1
    val frac = dt / (2 * h)
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

  private[sabr] def computeSystem(dt: Double, M1: Array[Double], M0: Array[Double], tri1: TridiagonalMatrix, tri0: TridiagonalMatrix) {

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
        //-(lower+middle-1+upper)/dt*2
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

  private[sabr] def computeP(): Array[Double] = {
    val Q = Array.ofDim[Double](size)
    Q(j0) = 1.0 / h
    return Q
  }


  private[sabr] def buildEm() {
    Em_ = Array.ofDim[Double](size)
    var j = size - 1
    while (j >= 0) {
      Em_(j) = 1.0
      j -= 1
    }
  }

  private[sabr] def buildEmCache(dt1: Double, dt2: Double) {
    buildEm()
    Em1Cache = Array.ofDim[Double](size)
    if (dt2 != 0) Em2Cache = Array.ofDim[Double](size)
    dt1cache = dt1
    dt2cache = dt2
    var j = size - 2
    val efactor = spec.rho * spec.nu * spec.alpha
    while (j > 0) {
      Em1Cache(j) = math.exp(efactor * Gammam_(j) * dt1)
      if (dt2 != 0) Em2Cache(j) = math.exp(efactor * Gammam_(j) * dt2)
      j -= 1
    }

    Em1Cache(0) = Em1Cache(1)
    Em1Cache(size - 1) = Em1Cache(size - 2)
    if (dt2 != 0) {
      Em2Cache(0) = Em2Cache(1)
      Em2Cache(size - 1) = Em2Cache(size - 2)
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
    val cfl = Cm_(j0) * (1 / (Fm_(j0 + 1) - Fm_(j0)) + 1 / (Fm_(j0) - Fm_(j0 - 1))) * dt / (2 * h)
    val cfl2 = math.pow(forward, spec.beta) * dt_ / (h_ * h_)
    println(cfl + " " + cfl2 + " " + math.abs(cfl - cfl2))
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
    val e = math.exp(spec.nu * z)
    val sh = 0.5 * (e - 1 / e)
    val ch = 0.5 * (e + 1 / e)
    return spec.alpha / spec.nu * (sh + spec.rho * (ch - 1))
  }

  private[sabr] final def computeFFromY(y: Double): Double = {
    val factor = (1 - spec.beta) * y
    return if (factor + forwardonebeta < 0) 0.0 else math.pow(forwardonebeta + factor, 1 / (1 - spec.beta))
  }

  def initGrid() {
    val sqrtT = math.sqrt(T)
    zmin_ = -nDeviations * sqrtT
    zmax_ = -zmin_
    if (spec.beta < 1) {
      val ybar = -forwardonebeta / (1 - spec.beta)
      val zbar = computeZFromY(ybar)
      if (zbar > zmin_) {
        zmin_ = zbar
      }
    }
    var j = 0
    val zm = Array.ofDim[Double](size)
    val h0 = (zmax_ - zmin_) / size
    j0 = math.max(1, (((0.0 - zmin_) / h0)).toInt)
    h_ = (0.0 - zmin_) / (j0 - 0.5)
    while (j < size) {
      zm(j) = zmin_ + (j - 0.5) * h
      j += 1
    }
    zmax_ = zm(size - 1) + 0.5 * h
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
      //if (Cm_(j).isInfinite) throw new RuntimeException("too large numbers: try a smaller number of deviations")
      Gammam_(j) = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) spec.beta / forwardonebeta else (fpow - forwardbeta) / (f - forward)
      j -= 1
    }

    Fmax = computeFFromY(computeYFromZ(zmax_))
    Fmin = computeFFromY(computeYFromZ(zmin_))
    Fm_(0) = 2 * Fmin - Fm_(1)
    Fm_(size - 1) = 2 * Fmax - Fm_(size - 2)
    Cm_(0) = Cm_(1)
    Cm_(size - 1) = Cm_(size - 2)
    //    println("Fmax = "+Fmax+" h="+h+" CFL="+computeCourantNumber)
  }

}
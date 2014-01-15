package quantscale.analytic

import quantscale.fdm.mesh.{Mesh1D, UniformMesh1D, Mesh1DBoundaries}
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.GridUtil
import quantscale.fdm.DifferentialCache
import org.apache.commons.math3.analysis.integration.RombergIntegrator
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator
import org.apache.commons.math3.complex.Complex
import quantscale.fdm.Epsilon
import quantscale.fdm.method.Parabolic1DBoundaryFactory
import quantscale.fdm.method.ForwardPartialOrder2Parabolic1DBoundaryFactory
import quantscale.fdm.OperatorLine
import quantscale.fdm.method.BackwardPartialOrder2Parabolic1DBoundaryFactory
import quantscale.math.{RungeKuttaODESolver, Function2D, CubicSpline, CubicPP}
import quantscale.fdm.TridiagonalSolver
import scala.math
import scala.math

/**
 * dF = v*(F+b)^beta dW1, F(0) = f
 * dv = nu*v*dW2, v(0)=alpha
 */
class SABRModelSpec(val alpha: Double, val beta: Double, val nu: Double, val rho: Double, val b: Double = 0.0) {

}

class ZABRModelSpec(val alpha: Double, val beta: Double, val nu: Double, val rho: Double, val gamma: Double) {

}

object SABRVanilla {
  private val AndreasenHugeSize = 1950

  def priceBenhamou(spec: SABRModelSpec, isCall: Boolean, strike: Double, forward: Double, tte: Double): Double = {
    val k = strike
    val f = forward
    val beta = spec.beta
    val alpha = spec.alpha
    val rho = spec.rho
    val nu = spec.nu
    val f0 = k;
    // projected fwd as suggested by Benhamou et al. makes z0=0
    val f1beta = math.pow(f, 1 - beta)
    val k1beta = math.pow(k, 1 - beta)
    val f01beta = k1beta
    val z = (f1beta - k1beta) / (alpha * (1 - beta));
    val z0 = (f01beta - k1beta) / (alpha * (1 - beta));

    val x = math.log((-rho + nu * z + math.sqrt(1 - 2 * nu * rho * z + nu * nu * z * z)) / (1 - rho)) / nu;
    val b1 = beta / (alpha * (1 - beta) * z0 + k1beta);

    val theta = 0.25 * rho * nu * alpha * b1 * z * z +
      math.log(alpha * math.pow(f * k, 0.5 * beta) * z / (f - k)) +
      math.log(x / z * math.pow(1 - 2 * nu * rho * z + nu * nu * z * z, 0.25))
    val kappa = 0.125 * ((alpha * alpha * (beta - 2) * beta * math.pow(k, 2 * beta - 2)) +
      6 * alpha * beta * math.pow(k, beta - 1) * nu * rho +
      nu * nu * (2 - 3 * rho * rho));

    //    val d1 = (k - alpha * (beta - 1) * math.pow(k, beta) * z0)
    //    val kappa = 0.125 * ((alpha * alpha * (beta - 2) * beta * math.pow(k, 2 * beta)) / (d1 * d1) +
    //      6 * alpha * beta * math.pow(k, beta) * nu * rho / (k - alpha * (beta - 1) * math.pow(k, beta) * z0) +
    //      nu * nu * (2 - 3 * rho * rho + 2 * nu * rho * z0 - nu * nu * z0 * z0) / (1 - 2 * nu * rho * z0 + nu * nu * z0 * z0));
    val sqrtkappa = new Complex(kappa, 0).sqrt().multiply(Complex.I)
    val i1 = new Complex(0.5 * math.sqrt(math.Pi), 0).divide(sqrtkappa)
    val i2 = sqrtkappa.multiply(math.sqrt(2) * x).exp()
    val i3 = new Complex(x / math.sqrt(2 * tte), 0).add(sqrtkappa.multiply(math.sqrt(tte)))
    val i4 = sqrtkappa.multiply(-math.sqrt(2) * x).exp()
    val i5 = new Complex(-x / math.sqrt(2 * tte), 0).add(sqrtkappa.multiply(math.sqrt(tte)))
    val integral = i1.multiply(i2.multiply(erfcomplex(i3).subtract(Complex.ONE)).add(i4.multiply(erfcomplex(i5).add(Complex.ONE))))

    val realIntegral = integral.getReal()
    val price = f - k + (f - k) / (2 * x * math.sqrt(2 * math.Pi)) * math.exp(theta) * realIntegral;
    return price
  }

  private def erfcomplex(z: Complex): Complex = {
    val x = z.getReal()
    val y = z.getImaginary()
    val erf = 2 * CumulativeNormalDistribution.value(x * math.sqrt(2)) - 1
    if (math.abs(y) < Epsilon.MACHINE_EPSILON) {
      return new Complex(erf, 0)
    }
    val firstTerm = new Complex(1 - math.cos(2 * x * y), math.sin(2 * x * y)).multiply(math.exp(-x * x) / (2 * math.Pi * x))
    var sum = Complex.ZERO
    var n = 1
    while (n < 6) {
      val term = new Complex(2 * x - 2 * x * math.cosh(n * y) * math.cos(2 * x * y) + n * math.sinh(n * y) * math.sin(2 * x * y),
        2 * x * math.cosh(n * y) * math.sin(2 * x * y) + n * math.sinh(n * y) * math.cos(2 * x * y)).multiply(math.exp(-0.25 * n * n) / (n * n + 4 * x * x))
      sum = sum.add(term)
      n += 1
    }
    sum.multiply(2 / math.Pi * math.exp(-x * x))
    return sum.add(firstTerm.add(erf))
  }

  def impliedVolatilityHagan(spec: SABRModelSpec, forward: Double, strike: Double, tte: Double): Double = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val logfk = math.log(forward / strike)
    val z = nu / alpha * math.pow(forward * strike, (1 - beta) / 2) * logfk
    //TODO rho = 1? beta = 1?
    val xz = math.log((math.sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho))
    if (math.abs(strike - forward) < 1e-8) {
      //ATM
      val f1beta = math.pow(forward, beta - 1)
      val sigma = alpha * f1beta * (1 + ((1 - beta) * (1 - beta) / 24 * alpha * alpha * f1beta * f1beta + 0.25 * rho * beta * alpha * nu * f1beta + (2 - 3 * rho * rho) / 24 * nu * nu) * tte)
      return sigma
    } else {
      val onebeta = 1 - beta
      val onebeta2 = onebeta * onebeta
      val onebeta4 = onebeta2 * onebeta2
      val fkbetaone = math.pow(forward * strike, -onebeta / 2)
      val firstTerm = alpha * fkbetaone / (1 + onebeta2 / 24 * logfk * logfk + onebeta4 / 1920 * logfk * logfk * logfk * logfk) * z / xz
      val secondTerm = 1 + (onebeta2 / 24 * alpha * alpha * fkbetaone * fkbetaone + 0.25 * rho * beta * nu * alpha * fkbetaone + (2 - 3 * rho * rho) / 24 * nu * nu) * tte
      val sigma = firstTerm * secondTerm
      return sigma
    }
  }

  def normalVolatilityHagan2013(spec: SABRModelSpec, forward: Double, strike: Double, tte: Double): Double = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val fonebeta = math.pow(forward, 1 - beta)
    val konebeta = math.pow(strike, 1 - beta)
    val xsi = nu / alpha * (fonebeta - konebeta) / (1 - beta)
    val x = math.log((math.sqrt(1 - 2 * rho * xsi + xsi * xsi) - rho + xsi) / (1 - rho))
    val g = (1 - beta) * (1 - beta) / ((fonebeta - konebeta) * (fonebeta - konebeta)) *
      math.log(math.pow(forward * strike, beta / 2) * (fonebeta - konebeta) / ((1 - beta) * (forward - strike)))
    val firstTerm = alpha * (1 - beta) * (forward - strike) / (fonebeta - konebeta)
    val secondTerm = 1 + tte * (g * alpha * alpha + rho * nu * alpha / 4 *
      (Math.pow(forward, beta) - Math.pow(strike, beta)) / (forward - strike) + (2 - 3 * rho * rho) / 24 * nu * nu)
    val sigma = firstTerm * xsi / x * secondTerm
    return sigma
  }

  def impliedVolatilityAndreasenHuge(spec: SABRModelSpec, forward: Double, strike: Double, tte: Double): Double = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val onebeta = 1 - beta
    val z = 1.0
    val y = 1.0 / (z * alpha * onebeta) * (math.pow(forward, onebeta) - math.pow(strike, onebeta))
    val jy = math.sqrt((1 + nu * nu * y * y - 2 * rho * nu * y))

    val x = 1.0 / nu * math.log((jy - rho + nu * y) / (1 - rho))
    val normalImpliedVol = (forward - strike) / x
    val lognormalImpliedVol = math.log(forward / strike) / x
    return lognormalImpliedVol
  }

  def normalVolatilityAndreasenHuge(spec: SABRModelSpec, forward: Double, strike: Double, tte: Double): Double = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val onebeta = 1 - beta
    val z = 1.0
    val y = 1.0 / (z * alpha * onebeta) * (math.pow(forward, onebeta) - math.pow(strike, onebeta))
    val jy = math.sqrt((1 + nu * nu * y * y - 2 * rho * nu * y))
    val x = 1.0 / nu * math.log((jy - rho + nu * y) / (1 - rho))
    val normalImpliedVol = (forward - strike) / x
    return normalImpliedVol
  }

  def priceAndreasenHuge(isCall: Boolean, strike: Double, spec: SABRModelSpec, forward: Double, tte: Double): Double = {
    val spline = priceAndreasenHuge(isCall, spec, forward, tte, true)
    return spline.value(strike)
  }

  def priceAndreasenHugeZABR(isCall: Boolean, strike: Double, spec: ZABRModelSpec, forward: Double, tte: Double): Double = {
    val spline = priceAndreasenHugeZABR(isCall, spec, forward, tte, true)
    return spline.value(strike)
  }

  def computeFmax(spec: SABRModelSpec, forward: Double, tte: Double): Double = {
    val dev = 4
    val theta = 0.5 * spec.nu * dev * math.sqrt(tte)
    val Fmax = forward * math.exp(2 * theta)
    val z = 2 / spec.nu * math.sinh(theta) * (math.cosh(theta) + spec.rho * math.sinh(theta))
    val fpow = z * (spec.alpha) * (1 - spec.beta) + math.pow(forward, 1 - spec.beta)
    val Fmax2 = math.pow(fpow, 1 / (1 - spec.beta))
    return Fmax
  }

  def priceAndreasenHuge(isCall: Boolean, spec: SABRModelSpec, forward: Double, tte: Double, withLaplaceCorrection: Boolean = true): CubicPP = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val onebeta = 1 - beta
    val z = 1.0
    //later the forward vol is used to price with 1 step euler FD scheme.
    //we could invert back this price to find real implied vol for comparison

    val fsteps = AndreasenHugeSize
    val fmin = 0
    var fmax = 10 * forward //math.min(10*forward,computeFmax(spec, forward, tte))
    var deltaf = (fmax - fmin) / fsteps
    val findex = ((forward - fmin) / deltaf).toInt
    deltaf = (forward - fmin) / findex
    fmax = (fsteps - 1) * deltaf
    val mesh = new UniformMesh1D(fsteps, new Mesh1DBoundaries(fmin, fmax))
    val rhs = Array.ofDim[Double](mesh.size)
    val lhs = new TridiagonalMatrix(mesh.size)
    val varianceVector = Array.ofDim[Double](mesh.size)
    val sign = if (isCall) 1 else -1
    val forwardonebeta = math.pow(forward, onebeta)
    for (i <- 0 until mesh.size) {
      val k = mesh.x(i)
      rhs(i) = math.max(sign * (forward - k), 0)
      val y = 1.0 / (z * alpha * onebeta) * (forwardonebeta - math.pow(k, onebeta))
      val jy = math.sqrt((1 + nu * nu * y * y - 2 * rho * nu * y))
      val sigmak = alpha * math.pow(k, beta)
      val forwardvol = jy * z * sigmak
      var thetak2 = forwardvol * forwardvol
      if (withLaplaceCorrection) {
        val x = 1.0 / nu * math.log((jy - rho + nu * y) / (1 - rho))
        val xsi = math.abs(x) / math.sqrt(tte)
        val factor = 2.0 * (1.0 - xsi * CumulativeNormalDistribution.value(-xsi) / NormalDistribution.value(xsi))
        //        val normalVol = (forward - k) / x
        //        val eps = 1e-7
        //        val factorter = (BachelierVanillaEuropean.price(isCall, k, forward, normalVol, tte) - math.max(sign * (forward - k), 0)) / (BachelierVanillaEuropean.price(isCall, k, forward, normalVol, tte + eps) - BachelierVanillaEuropean.price(isCall, k, forward, normalVol, tte)) * eps / tte
        //        if (math.abs(factorter - factor) > 1e-4) {
        //          println("factor was " + factor + " but factorter was " + factorter)
        //        }
        thetak2 *= factor
      }
      varianceVector(i) = thetak2
      lhs.middle(i) = 1.0
    }
    val diffCache = new DifferentialCache(mesh.x)
    lhs.plusD2(1, lhs.size - 1, diffCache, varianceVector, -0.5 * tte)
    //zero gamma boundary
    //    lhs.upper(0) = 0.0
    //    lhs.middle(0) = 1.0
    //    lhs.lower(lhs.size - 1) = 0.0
    //    lhs.middle(lhs.size - 1) = 1.0

    //andreasen huge boundary
    val normalImpliedVol0 = normalVolatilityAndreasenHuge(spec, forward, mesh.x(0), tte)
    rhs(0) = BachelierVanillaEuropean.price(isCall, mesh.x(0), forward, normalImpliedVol0, tte)
    val normalImpliedVoln = normalVolatilityAndreasenHuge(spec, forward, mesh.x(lhs.size - 1), tte)
    rhs(lhs.size - 1) = BachelierVanillaEuropean.price(isCall, mesh.x(lhs.size - 1), forward, normalImpliedVoln, tte)

    //order 1 boundary
    //    val    firstLine = new OperatorLine(0, 3)
    //    ForwardPartialOrder2Parabolic1DBoundaryFactory.makeODELine(0, mesh.x, -0.5*varianceVector(0)*tte, 0, 1.0, firstLine)
    //    val lastLine = new OperatorLine(lhs.size - 3,lhs.size) 
    //    BackwardPartialOrder2Parabolic1DBoundaryFactory.makeODELine(lhs.size - 1, mesh.x, -0.5*varianceVector(lhs.size-1)*tte, 0, 1.0, lastLine)
    //    lhs.setBoundaries(firstLine, lastLine)

    val price = Array.ofDim[Double](mesh.size)
    val solver = new ThomasTridiagonalSolver()
    solver.init(mesh.size)
    solver.solve(lhs, rhs, price)
    return CubicSpline.makeBesselSpline(mesh.x, price)
  }


  private def solveZABRODE(spec: ZABRModelSpec, y: Double): (Double, Double) = {
    val gamma = spec.gamma
    val nu = spec.nu
    val rho = spec.rho

    val odeFunction = new Function2D {
      def value(y: Double, f: Double): Double = {
        val A = 1 + (gamma - 2) * (gamma - 2) * nu * nu * y * y + 2 * rho * (gamma - 2) * nu * y
        val B = 2 * rho * (1 - gamma) * nu + 2 * (1 - gamma) * (gamma - 2) * nu * nu * y
        val C = (1 - gamma) * (1 - gamma) * nu * nu
        val delta = math.max(0, B * B * f * f - 4 * A * (C * f * f - 1))
        val oderhs = (-B * f + math.sqrt(delta)) / (2 * A)
        return oderhs
      }
    }
    val size = 50
    val ya = Array.ofDim[Double](size)
    var i = 0
    while (i < size) {
      ya(i) = i * y / (size - 1)
      i += 1
    }
    val f = Array.ofDim[Double](size)
    new RungeKuttaODESolver().solve(odeFunction, ya, 0, f)
    return (f(size - 1), odeFunction.value(y, f(size - 1)))
  }

  def priceAndreasenHugeZABR(isCall: Boolean, spec: ZABRModelSpec, forward: Double, tte: Double, withLaplaceCorrection: Boolean = true): CubicPP = {
    val strikePriceVol = priceAndreasenHugeZABRStrikePriceVol(isCall, spec, forward, tte, withLaplaceCorrection)
    return CubicSpline.makeBesselSpline(strikePriceVol._1, strikePriceVol._2)
  }

  def priceAndreasenHugeZABRStrikePriceVol(isCall: Boolean, spec: ZABRModelSpec, forward: Double, tte: Double, withLaplaceCorrection: Boolean = true): (Array[Double], Array[Double], Array[Double]) = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val onebeta = 1 - beta
    val gamma = spec.gamma
    val z = 1.0
    //later the forward vol is used to price with 1 step euler FD scheme.
    //we could invert back this price to find real implied vol for comparison

    val fsteps = AndreasenHugeSize
    var fmin = 0.005
    var fmax = 6 * (forward) + fmin;
    //math.min(10*forward,computeFmax(spec, forward, tte))
    var deltaf = (fmax - fmin) / fsteps
    //val findex = ((forward-fmin)/deltaf).toInt
    //deltaf = (forward-fmin)/findex


    fmax = (fsteps) * deltaf
    val mesh = new Mesh1D() {
      private val x_ = Array.ofDim[Double](fsteps)
      for (i <- 0 until fsteps) {
        x_(i) = fmin + (i) * deltaf
      }

      def x = x_

      def size = fsteps
    } //new UniformMesh1D(fsteps, new Mesh1DBoundaries(fmin, fmax))
    val rhs = Array.ofDim[Double](mesh.size)
    val lhs = new TridiagonalMatrix(mesh.size)
    val varianceVector = Array.ofDim[Double](mesh.size)
    val sign = if (isCall) 1 else -1
    val forwardonebeta = math.pow(forward, onebeta)

    var x = Array.ofDim[Double](mesh.size)
    for (i <- 0 until mesh.size) {
      val k = mesh.x(i)
      rhs(i) = math.max(sign * (forward - k), 0)
      val y = Math.pow(z, gamma - 2) / (alpha * onebeta) * (forwardonebeta - math.pow(k, onebeta))
      val fAndF = solveZABRODE(spec, y)
      val fy = fAndF._1
      val Fy = 1.0 / fAndF._2
      //val jy = Math.sqrt((1 + nu * nu * y * y - 2 * rho * nu * y))
      //println("k="+k+" Fy="+Fy+" jy="+jy)
      val sigmak = alpha * math.pow(k, beta)
      val forwardvol = Fy * z * sigmak
      var thetak2 = forwardvol * forwardvol
      x(i) = math.pow(z, 1 - gamma) * fy
      //val directBlackVol = if (math.abs(forward-k)<1e-6) alpha*math.pow(forward, beta-1) else (Math.log(forward/k))/x(i)
      //println("Direct "+gamma+" "+k+" "+directBlackVol)
      if (withLaplaceCorrection) {
        //1.0 / nu * math.log((jy - rho + nu * y) / (1 - rho))
        val xsi = math.abs(x(i)) / math.sqrt(tte)
        val factor = 2.0 * (1.0 - xsi * CumulativeNormalDistribution.value(-xsi) / NormalDistribution.value(xsi))
        //        val normalVol = (forward - k) / x
        //        val eps = 1e-7
        //        val factorter = (BachelierVanillaEuropean.price(isCall, k, forward, normalVol, tte) - math.max(sign * (forward - k), 0)) / (BachelierVanillaEuropean.price(isCall, k, forward, normalVol, tte + eps) - BachelierVanillaEuropean.price(isCall, k, forward, normalVol, tte)) * eps / tte
        //        if (math.abs(factorter - factor) > 1e-4) {
        //          println("factor was " + factor + " but factorter was " + factorter)
        //        }
        thetak2 *= factor
      }
      varianceVector(i) = thetak2
      lhs.middle(i) = 1.0
    }
    val diffCache = new DifferentialCache(mesh.x)
    lhs.plusD2(1, lhs.size - 1, diffCache, varianceVector, -0.5 * tte)
    //zero gamma boundary
    //    lhs.upper(0) = 0.0
    //    lhs.middle(0) = 1.0
    //    lhs.lower(lhs.size - 1) = 0.0
    //    lhs.middle(lhs.size - 1) = 1.0

    //andreasen huge boundary
    val normalImpliedVol0 = (forward - mesh.x(0)) / x(0);
    rhs(0) = BachelierVanillaEuropean.price(isCall, mesh.x(0), forward, normalImpliedVol0, tte)
    val normalImpliedVoln = (forward - mesh.x(lhs.size - 1)) / x(lhs.size - 1)
    rhs(lhs.size - 1) = BachelierVanillaEuropean.price(isCall, mesh.x(lhs.size - 1), forward, normalImpliedVoln, tte)

    //order 1 boundary
    //    val    firstLine = new OperatorLine(0, 3)
    //    ForwardPartialOrder2Parabolic1DBoundaryFactory.makeODELine(0, mesh.x, -0.5*varianceVector(0)*tte, 0, 1.0, firstLine)
    //    val lastLine = new OperatorLine(lhs.size - 3,lhs.size)
    //    BackwardPartialOrder2Parabolic1DBoundaryFactory.makeODELine(lhs.size - 1, mesh.x, -0.5*varianceVector(lhs.size-1)*tte, 0, 1.0, lastLine)
    //    lhs.setBoundaries(firstLine, lastLine)

    val price = Array.ofDim[Double](mesh.size)
    val solver = new ThomasTridiagonalSolver()
    solver.init(mesh.size)
    solver.solve(lhs, rhs, price)

    return (mesh.x, price, x)
  }

  def priceAndreasenHugeMulti(isCall: Boolean, spec: SABRModelSpec, forward: Double, tte: Double, timeSteps: Int): CubicPP = {
    val nu = spec.nu
    val alpha = spec.alpha
    val beta = spec.beta
    val rho = spec.rho
    val onebeta = 1 - beta
    val z = 1.0
    //later the forward vol is used to price with 1 step euler FD scheme.
    //we could invert back this price to find real implied vol for comparison
    //TODO make sure strike is on the grid
    val mesh = new UniformMesh1D(1600 * 1 + 1, new Mesh1DBoundaries(0, 4 * forward))
    var rhs = Array.ofDim[Double](mesh.size)
    val lhs = new TridiagonalMatrix(mesh.size)
    val varianceVector = Array.ofDim[Double](mesh.size)
    val sign = if (isCall) 1 else -1
    val forwardonebeta = math.pow(forward, onebeta)
    for (i <- 0 until mesh.size) {
      val k = mesh.x(i)
      rhs(i) = math.max(sign * (forward - k), 0)
      val y = 1.0 / (z * alpha * onebeta) * (forwardonebeta - math.pow(k, onebeta))
      val jy = math.sqrt((1 + nu * nu * y * y - 2 * rho * nu * y))
      val sigmak = alpha * math.pow(k, beta)
      val forwardvol = jy * z * sigmak
      var thetak2 = forwardvol * forwardvol
      varianceVector(i) = thetak2
      lhs.middle(i) = 1.0
    }
    val diffCache = new DifferentialCache(mesh.x)
    var t = tte
    val dt = tte / timeSteps
    lhs.plusD2(1, lhs.size - 1, diffCache, varianceVector, -0.5 * tte / timeSteps)
    //zero gamma boundary
    lhs.upper(0) = 0.0
    lhs.middle(0) = 1.0
    lhs.lower(lhs.size - 1) = 0.0
    lhs.middle(lhs.size - 1) = 1.0
    var price = Array.ofDim[Double](mesh.size)
    val solver = new ThomasTridiagonalSolver()
    solver.init(mesh.size)
    var tmp = price
    while (t > Epsilon.MACHINE_EPSILON_SQRT) {
      t -= dt
      price = tmp
      solver.solve(lhs, rhs, price)
      tmp = rhs
      rhs = price
    }

    //andreasen huge boundary
    //    val normalImpliedVol0 = normalVolatilityAndreasenHuge(spec, forward, mesh.x(0), tte)
    //    rhs(0) = BachelierVanillaEuropean.price(isCall, mesh.x(0), forward, normalImpliedVol0, tte)
    //    val normalImpliedVoln = normalVolatilityAndreasenHuge(spec, forward, mesh.x(lhs.size - 1), tte)
    //    rhs(lhs.size - 1) = BachelierVanillaEuropean.price(isCall, mesh.x(lhs.size - 1), forward, normalImpliedVoln, tte)

    //order 1 boundary
    //    val    firstLine = new OperatorLine(0, 3)
    //    ForwardPartialOrder2Parabolic1DBoundaryFactory.makeODELine(0, mesh.x, -0.5*varianceVector(0)*tte, 0, 1.0, firstLine)
    //    val lastLine = new OperatorLine(lhs.size - 3,lhs.size) 
    //    BackwardPartialOrder2Parabolic1DBoundaryFactory.makeODELine(lhs.size - 1, mesh.x, -0.5*varianceVector(lhs.size-1)*tte, 0, 1.0, lastLine)
    //    lhs.setBoundaries(firstLine, lastLine)

    return CubicSpline.makeBesselSpline(mesh.x, price)
  }
}

class DensityGaussianApproxPricer(spec: SABRModelSpec, forward: Double, tte: Double, steps: Int, Fmax: Double) {
  private var Q: Array[Double] = null
  private var F: Array[Double] = null
  private var dx = 0.0
  init()

  def init() {
    dx = (Fmax - spec.b) / steps
    var sum = 0.0
    var j = 0
    val forwardonebeta = math.pow(forward, 1 - spec.beta)
    F = Array.ofDim(steps)
    Q = Array.ofDim(steps)
    val C0 = math.pow(forward + spec.b, spec.beta)
    while (j < steps) {
      F(j) = spec.b + dx * (j + 0.5)
      val f = math.abs(F(j) + spec.b) - spec.b
      val C = math.pow(f + spec.b, spec.beta)
      val fonebeta = math.pow(f, 1 - spec.beta)
      val z = (fonebeta - forwardonebeta) / (spec.alpha * (1 - spec.beta))
      val y = 1 / spec.nu * math.log((math.sqrt(1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z) + spec.rho + spec.nu * z) / (1 + spec.rho))
      // val gamma = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) spec.beta / forwardonebeta else (C - C0) / (f - forward)
      // val expTerm = math.exp(spec.rho * spec.nu * spec.alpha * gamma * tte)
      val expTerm = 1
      val D = (1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z) * expTerm * C * C
      val dy = 1 / (spec.alpha * math.sqrt(D))
      Q(j) = 1 / math.sqrt(2 * math.Pi * tte) * math.exp(-y * y / (2 * tte)) * dy
      j += 1
    }
  }

  def price(isCall: Boolean, strike: Double): Double = {
    var j = 0
    var price = 0.0
    var Qsum = 0.0
    var j0 = -1
    if (isCall) {
      while (j < F.length) {
        if (F(j) > strike) {
          price += (F(j) - strike) * Q(j) * dx
          if (j0 < 0) j0 = j
        }
        Qsum += Q(j) * dx
        j += 1
      }
      val term = 0.5 * dx
      price += 0.5 * term * term * Q(j0)

    } else {
      while (j < F.length) {
        if (F(j) < strike) {
          price += (strike - F(j)) * Q(j) * dx
        } else {
          if (j0 < 0) j0 = j
        }
        Qsum += Q(j) * dx
        j += 1
      }
      val term = 0.5 * dx
      price += 0.5 * term * term * Q(j0)

    }
    return price
  }

}
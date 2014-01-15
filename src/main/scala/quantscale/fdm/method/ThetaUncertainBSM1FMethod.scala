package quantscale.fdm.method

import org.slf4j.LoggerFactory
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.UncertainBSM1FFDSpec
import quantscale.fdm.State

class ThetaUncertainBSM1FMethod(private var thetaV: Double = ThetaParabolic1DMethod.THETA_CRANK_NICOLSON) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private var tridiagonal: TridiagonalMatrix = null;
  private var rhs: Array[Double] = null;
  private var _spec: UncertainBSM1FFDSpec = null;
  private var x: Array[Double] = null;
  private var ex: Array[Double] = null;

  val tolerance = 1e-6
  val maxNewtonIteration = 3

  private var dtFracV = 1.0;

  override def spec = _spec

  def copy(): Parabolic1DMethod = {
    return new ThetaUncertainBSM1FMethod(thetaV)
  }

  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV.asInstanceOf[UncertainBSM1FFDSpec]
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    rhs = new Array[Double](spec.grid.spaceVector.size);
    x = spec.grid.spaceVector;
    ex = spec.grid.spaceTransform.transform(x);
  }

  def theta = thetaV;

  def theta_=(v: Double) = {
    thetaV = v;
  }

  def dtFrac = dtFracV;

  def dtFrac_=(v: Double) = {
    dtFracV = v;
  }

  def populateSystem(t: Double, dt: Double, f: Array[Double], fNew: Array[Double]) {
    var j: Int = x.length - 2;

    val isDriftSpaceDependent = spec.bIsStateDependent()
    var mu = if (isDriftSpaceDependent) 0.0 else spec.b(t, dt, 0);
    val r = -spec.c(t, dt, 0);
    while (j > 0) {
      //TODO optimize
      val xj = x(j);
      val dxj = x(j) - x(j - 1);
      val dxjplus = x(j + 1) - x(j);
      if (isDriftSpaceDependent) mu = spec.b(t, dt, ex(j))

      val gamma =
      //        (fNew(j + 1) - 2 * fNew(j) + fNew(j - 1)) / (dxj * dxjplus)
        2 * (fNew(j + 1) - fNew(j)) / ((x(j + 1) - x(j - 1)) * Math.abs(dxjplus)) +
          2 * (fNew(j - 1) - fNew(j)) / ((x(j + 1) - x(j - 1)) * Math.abs(dxj));

      var sigma = spec.vol(gamma)

      val variance = sigma * sigma * xj * xj;

      val sign = spec.volMin * spec.volMin * xj * xj / (dxj) - mu * xj

      var a = 0.0
      var c = 0.0
      var b = 0.0

      if (sign >= 0) {
        a = (mu * xj - variance / dxj) / (dxj + dxjplus)
        c = (mu * xj + variance / dxjplus) / (dxj + dxjplus)
        b = (variance / (dxj * dxjplus) + r)
        tridiagonal.middle(j) = 1 + (1 - theta) * dt * dtFrac * b;
        val betaTemp = (1 - theta) * dt * dtFrac;
        tridiagonal.lower(j) = betaTemp * a;
        tridiagonal.upper(j) = -betaTemp * c;
      } else {
        a = -variance / dxj / (dxj + dxjplus)
        c = (mu * xj / dxjplus + variance / (dxjplus * (dxj + dxjplus)))
        b = variance / (dxj * dxjplus) + r + mu * xj / dxjplus
        val betaTemp = (1 - theta) * dt * dtFrac
        tridiagonal.lower(j) = betaTemp * a;
        tridiagonal.upper(j) = -betaTemp * c;
        tridiagonal.middle(j) = 1 + (1 - theta) * dt * dtFrac * b;
      }
      if (tridiagonal.lower(j) > 0) {
        logger.error("found positive lower diag");
      }
      if (tridiagonal.upper(j) > 0) {
        logger.error("found positive upper diag");
      }
      if (tridiagonal.middle(j) < 0) {
        logger.error("diagonal negative, timestep is too large");
      }

      val bStar = 1 - theta * dt * dtFrac * b;
      if (theta == ThetaParabolic1DMethod.THETA_CRANK_NICOLSON) {
        rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
      } else {
        val betaTemp = theta * dt * dtFrac;
        val aStar = -betaTemp * a;
        val cStar = -betaTemp * c;
        rhs(j) = aStar * f(j - 1) + bStar * f(j) + cStar * f(j + 1);
      }
      j -= 1;
    }
    val mu0 = if (spec.isInstanceOf[ConstantLogBSM1FFDSpec]) spec.b(t, dt, ex(0)) + spec.vol(t, dt, ex(0)) * spec.vol(t, dt, ex(0)) / 2 else spec.b(t, dt, ex(0))
    initLowerBoundary(t, dt, f, x, r, mu0)
    val muM = if (spec.isInstanceOf[ConstantLogBSM1FFDSpec]) spec.b(t, dt, ex(x.length - 1)) + spec.vol(t, dt, ex(x.length - 1)) * spec.vol(t, dt, ex(x.length - 1)) / 2 else spec.b(t, dt, ex(x.length - 1))
    initUpperBoundary(t, dt, f, x, r, muM)
  }

  def initLowerBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val dx = x(1) - x(0);
    val c = dt * dtFrac * mu * x(0) / dx;
    val b = dt * dtFrac * (r + mu * x(0) / dx);
    tridiagonal.middle(0) = 1 + (1 - theta) * b;
    tridiagonal.upper(0) = -(1 - theta) * c;
    val rMiddle = 1 - theta * b;
    rhs(0) = rMiddle * f(0) + theta * c * f(1);
  }

  def initUpperBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val m = x.length - 1;
    val dx = x(m) - x(m - 1);
    val a = dt * dtFrac * mu * x(m) / dx;
    val b = dt * dtFrac * (r - mu * x(m) / dx);
    tridiagonal.middle(m) = 1 + (1 - theta) * b;
    tridiagonal.lower(m) = (1 - theta) * a;
    val rMiddle = 1 - theta * b;
    rhs(m) = -theta * a * f(m - 1) + rMiddle * f(m);
    if (logger.isDebugEnabled()) {
      logger.debug("dx=" + dx + " dt=" + dt + " mu=" + mu + " theta=" + theta + " xm=" + x(m));
    }
  }

  private def printGamma(f: Array[Double]) {
    val x = spec.grid.spaceVector

    var i = 1
    var gamma = new Array[Double](f.length)
    println("Spot\tGamma");
    while (i < f.length - 1) {
      gamma(i) =
        //2 * (f(i + 1) - f(i)) / ((x(i + 1) - x(i - 1)) * Math.abs(x(i+1)-x(i))) +
        //2 * (f(i - 1) - f(i)) / ((x(i + 1) - x(i - 1)) * Math.abs(x(i)-x(i-1)));

        (f(i + 1) - 2 * f(i) + f(i - 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)))
      println(x(i) + "\t" + gamma(i))
      i += 1
    }
    throw new RuntimeException("gamma after first step");

  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    val fInitial = f.price.clone()
    var error = 0.0
    var countNewton = 0
    do {
      error = 0.0
      val fNewtonIteration = f.price.clone();
      populateSystem(currentTime, dt, fInitial, f.price); //f used for gamma
      if (theta != ThetaParabolic1DMethod.THETA_EXPLICIT) {
        solver.solve(tridiagonal, rhs, f.price);
      }
      //            printGamma(f);

      var i = 1
      while (i < f.price.length - 1) {
        error = Math.max(
          Math.abs(f.price(i) - fNewtonIteration(i)) / Math.max(Math.abs(f.price(i)), 1.0),
          error)
        i += 1
      }
      countNewton += 1
      if (logger.isDebugEnabled()) {
        logger.debug("error=" + error);
      }
    } while (error > tolerance && countNewton < maxNewtonIteration);
    if (logger.isDebugEnabled()) {
      logger.debug("countNewton=" + countNewton);
    }
    //      if (countNewton >= maxNewtonIteration) {
    //        throw new RuntimeException("Newton method did not converge, error=" + error);
    //      }
  }
}
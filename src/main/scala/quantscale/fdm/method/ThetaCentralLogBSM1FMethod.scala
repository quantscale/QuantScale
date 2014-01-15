package quantscale.fdm.method {

import java.util.Arrays

import org.slf4j.LoggerFactory

import quantscale.fdm.BSM1FFDSpec
import quantscale.fdm.TridiagonalMatrix

object ThetaCentralLogBSM1FMethod {
  val THETA_CRANK_NICOLSON = 0.5
  val THETA_IMPLICIT = 0.0
  val THETA_EXPLICIT = 1.0
}

class ThetaCentralLogBSM1FMethod(var theta: Double = ThetaCentralLogBSM1FMethod.THETA_CRANK_NICOLSON) extends BSM1FMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private var tridiagonal: TridiagonalMatrix = null;
  private var rhs: Array[Double] = null;
  private var spec: BSM1FFDSpec = null;

  private var dtFracV = 1.0;
  private var x: Array[Double] = null;
  private var ex: Array[Double] = null;

  override def initSystem(specV: BSM1FFDSpec) {
    spec = specV
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    rhs = new Array[Double](spec.grid.spaceVector.size);
    x = spec.grid.spaceVector;
    ex = spec.grid.spaceTransform.transform(x);

  }

  def dtFrac = dtFracV;

  def dtFrac_=(v: Double) = {
    dtFracV = v;
  }

  def populateSystem(t: Double, dt: Double, f: Array[Double]) {
    val x = spec.grid.spaceVector;
    var j: Int = x.length - 2;

    val mu = spec.drift(t, dt);
    val r = spec.discountRate(t, dt);
    while (j > 0) {
      //TODO optimize
      val xj = x(j);
      val exj = ex(j);

      val sigma = spec.vol(t, dt, exj);
      val variance = sigma * sigma;
      val muMod = mu - 0.5 * variance;

      val dxj = x(j) - x(j - 1);
      val dxjplus = x(j + 1) - x(j);

      val b = (variance / (dxj * dxjplus) + r);
      val a = (muMod - variance / dxj)
      val c = (muMod + variance / dxjplus)
      tridiagonal.middle(j) = 1 + (1 - theta) * dt * dtFrac * b;
      val betaTemp = (1 - theta) * dt * dtFrac / (dxj + dxjplus);
      tridiagonal.lower(j) = betaTemp * a;
      tridiagonal.upper(j) = -betaTemp * c;

      val bStar = 1 - theta * dt * dtFrac * b;
      if (theta == ThetaParabolic1DMethod.THETA_CRANK_NICOLSON) {
        rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
      } else {
        val betaTemp = theta * dt * dtFrac / (dxj + dxjplus);
        val aStar = -betaTemp * a;
        val cStar = -betaTemp * c;
        rhs(j) = aStar * f(j - 1) + bStar * f(j) + cStar * f(j + 1);
      }
      j -= 1;
    }
    val sigma0 = spec.vol(t, dt, ex(0));
    initLowerBoundary(t, dt, f, x, r, mu - sigma0 * sigma0 * 0.5);
    val sigmaM = spec.vol(t, dt, ex(x.length - 1));
    initUpperBoundary(t, dt, f, x, r, mu - sigmaM * sigmaM * 0.5);
  }

  def initLowerBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val dx = x(1) - x(0);
    val c = dt * dtFrac * mu / dx;
    val b = dt * dtFrac * (r + mu / dx);
    tridiagonal.middle(0) = 1 + (1 - theta) * b;
    tridiagonal.upper(0) = -(1 - theta) * c;
    val rMiddle = 1 - theta * b;
    rhs(0) = rMiddle * f(0) + theta * c * f(1);
  }

  def initUpperBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val m = x.length - 1;
    val dx = x(m) - x(m - 1);
    val a = dt * dtFrac * mu / dx;
    val b = dt * dtFrac * (r - mu / dx);
    tridiagonal.middle(m) = 1 + (1 - theta) * b;
    tridiagonal.lower(m) = (1 - theta) * a;
    val rMiddle = 1 - theta * b;
    rhs(m) = -theta * a * f(m - 1) + rMiddle * f(m);
    if (logger.isDebugEnabled()) {
      logger.debug("dx=" + dx + " dt=" + dt + " mu=" + mu + " theta=" + theta + " xm=" + x(m));
    }
  }

  override def solve(currentTime: Double, dt: Double, f: Array[Double]) {
    populateSystem(currentTime, dt, f);
    if (logger.isDebugEnabled()) {
      logger.debug("a=" + Arrays.toString(tridiagonal.lower));
      logger.debug("b=" + Arrays.toString(tridiagonal.middle));
      logger.debug("c=" + Arrays.toString(tridiagonal.upper));
      logger.debug("d=" + Arrays.toString(rhs));
    }
    if (theta != ThetaParabolic1DMethod.THETA_EXPLICIT) {
      solver.solve(tridiagonal, rhs, f);
    }
  }
}

}
package quantscale.fdm.method
import quantscale.fdm.UncertainBSM1FFDSpec
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.UncertainBSM1FFDSpec
import org.slf4j.LoggerFactory
import quantscale.fdm.Epsilon
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State

class TRBDF2UncertainBSM1FMethod(payoff: FDPayoff) extends Parabolic1DMethod {
  private final val logger = LoggerFactory.getLogger(getClass());

  final val theta = 0.5;
  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = (1 - theta) * alpha;

  private var tridiagonal: TridiagonalMatrix = null;
  private var rhs: Array[Double] = null;
  private var _spec: UncertainBSM1FFDSpec = null;
  private var fTr: Array[Double] = null;
  private var x: Array[Double] = null;
  private var ex: Array[Double] = null;

  val tolerance = 1e-6
  val maxNewtonIteration = 3

  override def spec = _spec

  def copy(): Parabolic1DMethod = {
    return new TRBDF2UncertainBSM1FMethod(payoff)
  }
  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV.asInstanceOf[UncertainBSM1FFDSpec];
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    rhs = new Array[Double](spec.grid.spaceVector.size);
    fTr = new Array[Double](spec.grid.spaceVector.size);
    x = spec.grid.spaceVector;
    ex = spec.grid.spaceTransform.transform(x);
  }

  def populateSystem(t: Double, dt: Double, f: Array[Double], fNew: Array[Double]) {
    val beta = backcoeff * dt;
    var j: Int = x.length - 2;

    val isDriftSpaceDependent = spec.bIsStateDependent
    val mu = if (isDriftSpaceDependent) 0.0 else spec.b(t, dt, 0);
    val r = -spec.c(t, dt, 0);
    while (j > 0) {
      //TODO optimize
      val xj = x(j);

      val dxj = x(j) - x(j - 1);
      val dxjplus = x(j + 1) - x(j);
      val muMod = if (isDriftSpaceDependent) spec.b(t, dt, ex(j)) else mu;

      val gamma =
        //        (fNew(j + 1) - 2 * fNew(j) + fNew(j - 1)) / (dxj * dxjplus)
        2 * (fNew(j + 1) - fNew(j)) / ((x(j + 1) - x(j - 1)) * Math.abs(dxjplus)) +
          2 * (fNew(j - 1) - fNew(j)) / ((x(j + 1) - x(j - 1)) * Math.abs(dxj));

      var sigma = spec.vol(gamma)

      val variance = sigma * sigma * xj * xj;

      val sign = spec.volMin * spec.volMin * xj * xj / (dxj) - muMod * xj
      //      val dtMax = 1.0/(alpha)/(spec.volMax*spec.volMax*xj*xj/(dxj*dxjplus)+mu*xj/(dxjplus+dxj)+r);
      //      if (dt > dtMax ) {
      //        throw new RuntimeException("schema will not be monotone, max dt="+dtMax+" but was "+dt);
      //      }
      val betaTemp = beta / (dxj + dxjplus);
      if (sign >= 0) {
        tridiagonal.lower(j) = betaTemp * (muMod * xj - variance / dxj);
        tridiagonal.upper(j) = -betaTemp * (muMod * xj + variance / dxjplus);
        tridiagonal.middle(j) = 1 + beta * (variance / (dxj * dxjplus) + r);
        val bStar = 1 - theta * alpha * dt * (variance / (dxj * dxjplus) + r);
        rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
      } else {
        tridiagonal.lower(j) = -betaTemp * (variance / dxj);
        tridiagonal.upper(j) = -beta * (muMod * xj / dxjplus + variance / (dxjplus * (dxj + dxjplus)));
        tridiagonal.middle(j) = 1 + beta * (variance / (dxj * dxjplus) + r + muMod * xj / dxjplus);
        val bStar = 1 - theta * alpha * dt * (variance / (dxj * dxjplus) + r + muMod * xj / dxjplus);
        rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
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

      //      tridiagonal.middle(j) = 1 + beta * (variance / (dxj * dxjplus) + r);
      //      val betaTemp = beta / (dxj + dxjplus);
      //      tridiagonal.lower(j) = betaTemp * (mu * xj - variance / dxj);
      //      tridiagonal.upper(j) = -betaTemp * (mu * xj + variance / dxjplus);
      //
      //      val bStar = 1 - theta * alpha * dt * (variance / (dxj * dxjplus) + r);
      //      rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
      j -= 1;
    }
    val mu0 = if (spec.isInstanceOf[ConstantLogBSM1FFDSpec]) spec.b(t, dt, ex(0)) + spec.vol(t, dt, ex(0)) * spec.vol(t, dt, ex(0)) / 2 else spec.b(t, dt, ex(0))
    initLowerBoundary(t, dt, f, x, r, mu0)
    val muM = if (spec.isInstanceOf[ConstantLogBSM1FFDSpec]) spec.b(t, dt, ex(x.length - 1)) + spec.vol(t, dt, ex(x.length - 1)) * spec.vol(t, dt, ex(x.length - 1)) / 2 else spec.b(t, dt, ex(x.length - 1))
    initUpperBoundary(t, dt, f, x, r, muM)
  }

  def initLowerBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val dx = x(1) - x(0);
    tridiagonal.middle(0) = 1 + backcoeff * dt * (r + mu * x(0) / dx);
    tridiagonal.upper(0) = -backcoeff * dt * mu * x(0) / dx;
    val rMiddle = 1 - backcoeff * dt * (r + mu * x(0) / dx);
    rhs(0) = rMiddle * f(0) - tridiagonal.upper(0) * f(1);
  }

  def initUpperBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val m = x.length - 1;
    val dx = x(m) - x(m - 1);
    tridiagonal.middle(m) = 1 + backcoeff * dt * (r - mu * x(m) / dx);
    tridiagonal.lower(m) = backcoeff * dt * mu * x(m) / dx;
    val rMiddle = 1 - backcoeff * dt * (r - mu * x(m) / dx);
    rhs(m) = -tridiagonal.lower(m) * f(m - 1) + rMiddle * f(m);
  }

  def computeBDF2Rhs(price: Array[Double], prevPrice: Array[Double], rhs: Array[Double]) {
    val frac = 1 / (alpha * (2 - alpha));
    var i: Int = spec.grid.spaceSize - 1;
    while (i >= 0) {
      rhs(i) = frac * (price(i) - (1 - alpha) * (1 - alpha) * prevPrice(i));
      i -= 1;
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    solve1(currentTime, dt, f.price)
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

        (f(i + 1) - 2 * f(i) + f(i - 1)) / ((ex(i + 1) - ex(i)) * (ex(i) - ex(i - 1)))
      println(x(i) + "\t" + gamma(i))
      i += 1
    }
    throw new RuntimeException("gamma after first step");

  }

  def solve1(currentTime: Double, dt: Double, f: Array[Double]) {
    val fInitial = f.clone()
    var error = 0.0
    var countNewton = 0
    do {
      error = 0.0
      val fNewtonIteration = f.clone();
      populateSystem(currentTime, dt, fInitial, f); //f used for gamma
      if (payoff != null) payoff.setTime(currentTime + dt * alpha);
      solver.solve(tridiagonal, rhs, fTr);
      computeBDF2Rhs(fTr, fInitial, rhs);
      if (payoff != null) payoff.setTime(currentTime);
      solver.solve(tridiagonal, rhs, f);

      var i = 1
      while (i < f.length - 1) {
        error = Math.max(
          Math.abs(f(i) - fNewtonIteration(i)) / Math.max(Math.abs(f(i)), 1.0),
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
    //    if (countNewton >= maxNewtonIteration) {
    //      throw new RuntimeException("Newton method did not converge, error=" + error);
    //    }
    //TODO compute next dt. The f & q in bank are probably L & f. TODO verify
    //    val C = (-3 * alpha * alpha + 4 * alpha - 2) / (12 * (2 - alpha))
    //    val tau = new Array[Double](f.length)
    //    var i = 0
    //    while (i < tau.length) {
    //      tau(i) = 2 * C * dt * (fInitial(i) / alpha - fTr(i) / (alpha * (1 - alpha)) +
    //        f(i) / (1 - alpha))
    //      i += 1
    //    }

  }

  def solve2(currentTime: Double, dt: Double, f: Array[Double]) {
    val fInitial = f.clone()
    var error = 0.0
    var countNewton = 0

    populateSystem(currentTime, dt, fInitial, fTr); //f used for gamma
    if (payoff != null) payoff.setTime(currentTime + dt * alpha);
    solver.solve(tridiagonal, rhs, fTr);

    //    do {
    //      error = 0.0
    //      val fNewtonIteration = if (countNewton ==0) f else fTr.clone();
    //      populateSystem(currentTime, dt, fInitial, fTr); //f used for gamma
    //      if (payoff != null) payoff.setTime(currentTime + dt * alpha);
    //      solver.solve(tridiagonal, rhs, fTr);
    //
    //      var i = 1
    //      while (i < f.length - 1) {
    //        error = Math.max(
    //          Math.abs(fTr(i) - fNewtonIteration(i)) / Math.max(Math.abs(fTr(i)), 1.0),
    //          error)
    //        i += 1
    //      }
    //      countNewton += 1
    //      if (logger.isDebugEnabled()) {
    //        logger.debug("TR error=" + error);
    //      }
    //    } while (error > tolerance && countNewton < maxNewtonIteration);
    //    if (logger.isDebugEnabled()) {
    //      logger.debug("TR countNewton=" + countNewton);
    //    }
    //    if (countNewton >= maxNewtonIteration) {
    //      throw new RuntimeException("TR Newton method did not converge, error=" + error);
    //    }

    countNewton = 0
    do {
      error = 0.0
      val fNewtonIteration = if (countNewton == 0) fTr else f.clone()
      populateSystem(currentTime, dt, fInitial, f); //f used for gamma
      computeBDF2Rhs(fTr, fInitial, rhs);
      if (payoff != null) payoff.setTime(currentTime);
      solver.solve(tridiagonal, rhs, f);
      var i = 1
      while (i < f.length - 1) {
        error = Math.max(
          Math.abs(f(i) - fNewtonIteration(i)) / Math.max(Math.abs(f(i)), 1.0),
          error)
        i += 1
      }
      countNewton += 1
      if (logger.isDebugEnabled()) {
        logger.debug("BDF2 error=" + error);
      }
    } while (error > tolerance && countNewton < maxNewtonIteration);
    if (logger.isDebugEnabled()) {
      logger.debug("BDF2 countNewton=" + countNewton);
    }
    if (countNewton >= maxNewtonIteration) {
      throw new RuntimeException("BDF2 Newton method did not converge, error=" + error);
    }

  }

}
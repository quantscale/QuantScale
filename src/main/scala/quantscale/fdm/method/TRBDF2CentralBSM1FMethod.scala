package quantscale.fdm.method

import quantscale.fdm.BSM1FFDSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff

class TRBDF2CentralBSM1FMethod(payoff: FDPayoff) extends BSM1FMethod {

  final val theta = 0.5;
  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = (1 - theta) * alpha;

  private var tridiagonal: TridiagonalMatrix = null;
  private var rhs: Array[Double] = null;
  private var spec: BSM1FFDSpec = null;
  private var fTr: Array[Double] = null;

  override def initSystem(specV: BSM1FFDSpec) {
    spec = specV
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    rhs = new Array[Double](spec.grid.spaceVector.size);
    fTr = new Array[Double](spec.grid.spaceVector.size);
  }

  def populateSystem(t: Double, dt: Double, f: Array[Double]) {
    val beta = backcoeff * dt;
    val x = spec.grid.spaceVector;
    var j: Int = x.length - 2;

    val mu = spec.drift(t, dt);
    val r = spec.discountRate(t, dt);
    while (j > 0) {
      //TODO optimize
      val xj = x(j);
      val sigma = spec.vol(t, dt, xj);
      val variance = sigma * sigma * xj * xj;
      val dxj = x(j) - x(j - 1);
      val dxjplus = x(j + 1) - x(j);

      tridiagonal.middle(j) = 1 + beta * (variance / (dxj * dxjplus) + r);
      val betaTemp = beta / (dxj + dxjplus);
      tridiagonal.lower(j) = betaTemp * (mu * xj - variance / dxj);
      tridiagonal.upper(j) = -betaTemp * (mu * xj + variance / dxjplus);

      val bStar = 1 - theta * alpha * dt * (variance / (dxj * dxjplus) + r);
      rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
      j -= 1;
    }
    initLowerBoundary(t, dt, f, x, r, mu);
    initUpperBoundary(t, dt, f, x, r, mu);
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

  override def solve(currentTime: Double, dt: Double, f: Array[Double]) {
    val fInitial = f.clone();
    populateSystem(currentTime, dt, f);
    if (payoff != null) payoff.setTime(currentTime + dt * alpha);
    solver.solve(tridiagonal, rhs, fTr);
    computeBDF2Rhs(fTr, f, rhs);
    if (payoff != null) payoff.setTime(currentTime);
    solver.solve(tridiagonal, rhs, f);

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

}
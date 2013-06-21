package quantscale.fdm.method
import quantscale.fdm.BSM1FFDSpec
import quantscale.fdm.TridiagonalMatrix

/**
 * Solves the log transformed problem: the space contains log(assetPrice)
 */
class TRBDF2CentralLogBSM1FMethod() extends BSM1FMethod {

  final val theta = 0.5;
  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = (1 - theta) * alpha;

  private var tridiagonal: TridiagonalMatrix = null;
  private var rhs: Array[Double] = null;
  private var spec: BSM1FFDSpec = null;
  private var fTr: Array[Double] = null;
  private var x : Array[Double] = null;
  private var ex: Array[Double] = null;
  
  override def initSystem(specV: BSM1FFDSpec) {
    spec = specV
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    rhs = new Array[Double](spec.grid.spaceVector.size);
    fTr = new Array[Double](spec.grid.spaceVector.size);
    x = spec.grid.spaceVector;
    ex = spec.grid.spaceTransform.transform(x);
  }

  def populateSystem(t: Double, dt: Double, f: Array[Double]) {
    val beta = backcoeff * dt;
    var j: Int = x.length - 2;

    val mu = spec.drift(t, dt);
    val r = spec.discountRate(t, dt);
    while (j > 0) {
      //TODO optimize
      val xj = x(j);
      val exj = ex(j);
      val sigma = spec.vol(t, dt, exj);
      val variance = sigma * sigma;
      val dxj = x(j) - x(j - 1);
      val dxjplus = x(j + 1) - x(j);
      val muMod = mu - 0.5 * variance;

      tridiagonal.middle(j) = 1 + beta * (variance / (dxj * dxjplus) + r);
      val betaTemp = beta / (dxj + dxjplus);
      tridiagonal.lower(j) = betaTemp * (muMod - variance / dxj);
      tridiagonal.upper(j) = -betaTemp * (muMod + variance / dxjplus);

      val bStar = 1 - theta * alpha * dt * (variance / (dxj * dxjplus) + r);
      rhs(j) = -tridiagonal.lower(j) * f(j - 1) + bStar * f(j) - tridiagonal.upper(j) * f(j + 1);
      j -= 1;
    }
    val sigma0 = spec.vol(t, dt, ex(0));
    initLowerBoundary(t, dt, f, x, r, mu - sigma0 * sigma0 * 0.5);
    val sigmaM = spec.vol(t, dt, ex(x.length - 1));
    initUpperBoundary(t, dt, f, x, r, mu - sigmaM * sigmaM * 0.5);
  }

  def initLowerBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val dx = x(1) - x(0);
    tridiagonal.middle(0) = 1 + backcoeff * dt * (r + mu / dx);
    tridiagonal.upper(0) = -backcoeff * dt * mu / dx;
    val rMiddle = 1 - backcoeff * dt * (r + mu / dx);
    rhs(0) = rMiddle * f(0) - tridiagonal.upper(0) * f(1);
  }

  def initUpperBoundary(t: Double, dt: Double, f: Array[Double], x: Array[Double], r: Double, mu: Double) {
    val m = x.length - 1;
    val dx = x(m) - x(m - 1);
    tridiagonal.middle(m) = 1 + backcoeff * dt * (r - mu / dx);
    tridiagonal.lower(m) = backcoeff * dt * mu / dx;
    val rMiddle = 1 - backcoeff * dt * (r - mu / dx);
    rhs(m) = -tridiagonal.lower(m) * f(m - 1) + rMiddle * f(m);
  }

  def computeBDF2Rhs(price: Array[Double], prevPrice: Array[Double]) {
    val frac = 1 / (alpha * (2 - alpha));
    var i: Int = spec.grid.spaceSize - 1;
    while (i >= 0) {
      rhs(i) = frac * (price(i) - (1 - alpha) * (1 - alpha) * prevPrice(i));
      i -= 1;
    }
  }

  override def solve(currentTime: Double, dt: Double, f: Array[Double]) {
    populateSystem(currentTime, dt, f);
    solver.solve(tridiagonal, rhs, fTr);
    computeBDF2Rhs(fTr, f);
    solver.solve(tridiagonal, rhs, f);
  }

}
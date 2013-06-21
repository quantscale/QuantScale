package quantscale.fdm.method
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State

class TRBDF2SingleCentralBSM1FMethod() extends Parabolic1DMethod {

  final val theta = 0.5;
  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = (1 - theta) * alpha;

  private var tridiagonal: TridiagonalMatrix = null;
  private var rhs: State = null;
  private var _spec: Parabolic1DFDSpec = null;
  private var fTr: State = null;
  private var thomasSolver = new ThomasTridiagonalSolver();
  private val frac = 1 / (alpha * (2 - alpha));
  private var x: Array[Double] = null;
  private var ex: Array[Double] = null;

  override def spec = _spec
      def copy() : Parabolic1DMethod = {
    return new TRBDF2SingleCentralBSM1FMethod
  }
  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    thomasSolver.init(spec.grid.spaceVector.size);
    x = spec.grid.spaceVector;
    ex = spec.grid.spaceTransform.transform(x);

  }

  def populateSystem(t: Double, dt: Double, f: State) {
    val beta = backcoeff * dt;
    var j: Int = x.length - 2;

    val isDriftSpaceDependent = spec.bIsStateDependent()
    val mu = if (isDriftSpaceDependent) 0.0 else spec.b(t, dt, 0);
    val rC = if (spec.cIsStateDependent) 0.0 else -spec.c(t, dt, 0);
    while (j > 0) {
      //TODO optimize
      val variance = spec.a(t, dt, ex(j)) * 2
      val dxj = x(j) - x(j - 1);
      val dxjplus = x(j + 1) - x(j);
      val muMod = if (isDriftSpaceDependent) spec.b(t, dt, ex(j)) else mu
      val r = if (spec.cIsStateDependent) -spec.c(t, dt, ex(j)) else rC
      if (smearingReducer != null) {
        smearingReducer.filter(variance, muMod, dxj, dxjplus)
      }
      tridiagonal.middle(j) = 1 + beta * (variance / (dxj * dxjplus) + r);
      val betaTemp = beta / (dxj + dxjplus);
      tridiagonal.lower(j) = betaTemp * (muMod - variance / dxj);
      tridiagonal.upper(j) = -betaTemp * (muMod + variance / dxjplus);

      val bStar = 1 - theta * alpha * dt * (variance / (dxj * dxjplus) + r);
      for (d <- 0 until f.stateDimensions) {
        rhs.values(d)(j) = -tridiagonal.lower(j) * f.values(d)(j - 1) + bStar * f.values(d)(j) - tridiagonal.upper(j) * f.values(d)(j + 1);
      }
      j -= 1;
    }
    val mu0 = if (spec.isInstanceOf[ConstantLogBSM1FFDSpec]) spec.b(t, dt, ex(0)) + spec.a(t, dt, ex(0)) else spec.b(t, dt, ex(0))
    initLowerBoundary(t, dt, f, x, -spec.c(t, dt, ex(0)), mu0)
    val muM = if (spec.isInstanceOf[ConstantLogBSM1FFDSpec]) spec.b(t, dt, ex(x.length - 1)) + spec.a(t, dt, ex(x.length - 1)) else spec.b(t, dt, ex(x.length - 1))
    initUpperBoundary(t, dt, f, x, -spec.c(t, dt, ex(ex.length - 1)), muM)
  }

  def initLowerBoundary(t: Double, dt: Double, f: State, x: Array[Double], r: Double, mu: Double) {
    val dx = x(1) - x(0);
    val sign = -1 //important for Brennan Schwartz, upwinding would be if (mu > 0) 1 else -1
    tridiagonal.middle(0) = 1 - sign * backcoeff * dt * (r + mu / dx);
    tridiagonal.upper(0) = sign * backcoeff * dt * mu / dx;
    val rMiddle = 1 + sign * backcoeff * dt * (r + mu / dx);
    for (d <- 0 until f.stateDimensions) {
      rhs.values(d)(0) = rMiddle * f.values(d)(0) - tridiagonal.upper(0) * f.values(d)(1);
    }
  }

  def initUpperBoundary(t: Double, dt: Double, f: State, x: Array[Double], r: Double, mu: Double) {
    val m = x.length - 1;
    val dx = x(m) - x(m - 1);
    val sign = -1 //if (mu > 0) -1 else 1
    tridiagonal.middle(m) = 1 + sign * backcoeff * dt * (-r + mu / dx);
    tridiagonal.lower(m) = -sign * backcoeff * dt * mu / dx;
    val rMiddle = 1 - sign * backcoeff * dt * (-r + mu / dx);
    for (d <- 0 until f.stateDimensions) {
      rhs.values(d)(m) = -tridiagonal.lower(m) * f.values(d)(m - 1) + rMiddle * f.values(d)(m);
    }
  }

  def computeBDF2Rhs(price: Array[Double], prevPrice: Array[Double], rhs: Array[Double]) {
    var i: Int = spec.grid.spaceSize - 1;
    while (i >= 0) {
      rhs(i) = frac * (price(i) - (1 - alpha) * (1 - alpha) * prevPrice(i));
      i -= 1;
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    if (rhs == null) {
      rhs = new State(f.stateDimensions, f.size)
      fTr = new State(f.stateDimensions, f.size)
    }
    populateSystem(currentTime, dt, f);
    for (d <- 0 until f.stateDimensions) {
      thomasSolver.solve(tridiagonal, rhs.values(d), fTr.values(d));
      computeBDF2Rhs(fTr.values(d), f.values(d), rhs.values(d));
      solver.solve(tridiagonal, rhs.values(d), f.values(d));
    }
  }

}
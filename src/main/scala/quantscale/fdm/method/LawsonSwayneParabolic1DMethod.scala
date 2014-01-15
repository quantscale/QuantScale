package quantscale.fdm.method

import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm._

class LawsonSwayneParabolic1DMethod(payoff: FDPayoff = null) extends Parabolic1DMethod {
  private val b = 1 - math.sqrt(2) / 2
  private val sqrt2 = math.sqrt(2)

  var tridiagonal: TridiagonalMatrix = null;
  private var rhs: State = null;
  private var _spec: Parabolic1DFDSpec = null;
  private var fTr: State = null;
  private var x: Array[Double] = null;
  var ex: Array[Double] = null;

  var driftVector: Array[Double] = null
  var varianceVector: Array[Double] = null
  var discountVector: Array[Double] = null
  var diffCache: DifferentialCache = null
  private var firstLine: OperatorLine = null
  private var lastLine: OperatorLine = null
  private var solverND: TridiagonalSolverND = null

  override def spec = _spec

  def copy(): Parabolic1DMethod = {
    return new LawsonSwayneParabolic1DMethod(payoff)
  }

  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV
    tridiagonal = new TridiagonalMatrix(specV.grid.spaceVector.size);
    x = specV.grid.spaceVector;
    driftVector = new Array[Double](x.length)
    discountVector = new Array[Double](x.length)
    varianceVector = new Array[Double](x.length)
    ex = specV.grid.spaceTransform.transform(x);
    diffCache = new DifferentialCache(x)
    firstLine = new OperatorLine(0, 3)
    lastLine = new OperatorLine(x.length - 3, x.length)
    fTr = null
    rhs = null
    solverND = null
  }

  def initLeftHandSide(t: Double, dt: Double) {
    computeDiscount(ex, t, dt, discountVector)
    computeDrift(ex, t, dt, driftVector)
    computeVariance(ex, t, dt, driftVector, varianceVector)
    val multiplier = -dt
    val A = tridiagonal
    A.fill(0.0)
    A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier)
      .plusD1(1, x.length - 1, diffCache, driftVector, multiplier)
      .plusD0(1, x.length - 1, discountVector, -multiplier, 1.0)
    //    A.parabolicOperator(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier, driftVector, multiplier, 1.0 - multiplier * r)
  }

  def initBoundaries(t: Double, dt: Double) {
    val multiplier = -dt
    lowerBoundary.makeLine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), multiplier, firstLine)
    val m = x.length - 1
    upperBoundary.makeLine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), multiplier, lastLine)
    tridiagonal.setBoundaries(firstLine, lastLine)
  }

  //supposes that left hand side is initialized
  def initRightHandSide(f: State) {
    if (rhs == null) rhs = new State(f.stateDimensions, f.size)
    for (d <- 0 until f.stateDimensions) {
      Array.copy(f.values(d), 0, rhs.values(d), 0, f.size)
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    initLeftHandSide(currentTime + dt * (1 - b), dt * b)
    initRightHandSide(f)
    initBoundaries(currentTime + dt * (1 - b), dt * b)
    if (payoff != null) payoff.setTime(currentTime + dt * (1 - b))
    if (fTr == null) fTr = new State(f.stateDimensions, f.size)
    if (solverND == null) {
      solverND = new TridiagonalSolverND(solver, f.stateDimensions)
    }
    solverND.solve(tridiagonal, rhs.values, fTr.values)
    initLeftHandSide(currentTime + dt * (1 - 2 * b), dt * b)
    //initRightHandSide(f)
    initBoundaries(currentTime + dt * (1 - 2 * b), dt * b)
    if (payoff != null) payoff.setTime(currentTime); //or +dt*2b?
    solverND.solve(tridiagonal, fTr.values, f.values);
    for (d <- 0 until rhs.stateDimensions) {
      var i = 0
      val fValues = f.values(d)
      val fTrValues = fTr.values(d)
      val size = f.size
      while (i < size) {
        fValues(i) = (sqrt2 + 1) * fValues(i) - sqrt2 * fTrValues(i)
        i += 1
      }
    }
  }
}

package quantscale.fdm.method
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.OperatorLine
import quantscale.fdm.DifferentialCache
import quantscale.fdm.State
import quantscale.fdm.TridiagonalSolverND
import quantscale.fdm.Epsilon

/**
 * Solves the log transformed problem: the space contains log(assetPrice)
 */
class TGASerialParabolic1DMethod(payoff: FDPayoff = null, private var _a: Double = 1.25 + math.sqrt(2) / 2) extends Parabolic1DMethod {

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

  private var r1, r2: Double = 0.0

  override def spec = _spec

  def copy(): Parabolic1DMethod = {
    return new TGASerialParabolic1DMethod(payoff)
  }
  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV
    tridiagonal = new TridiagonalMatrix(spec.grid.spaceVector.size);
    x = spec.grid.spaceVector;
    driftVector = new Array[Double](x.length)
    discountVector = new Array[Double](x.length)
    varianceVector = new Array[Double](x.length)
    ex = spec.grid.spaceTransform.transform(x);
    diffCache = new DifferentialCache(x)
    firstLine = new OperatorLine(0, 3)
    lastLine = new OperatorLine(x.length - 3, x.length)
    initFactors()
  }

  private def initFactors() {
    val sqrt = math.sqrt(_a * _a - 4 * _a + 2)
    r1 = (2 * _a - 1) / (_a + sqrt)
    r2 = (2 * _a - 1) / (_a - sqrt)
    //    println(_a+" "+sqrt+" "+r1+" "+r2+" "+s1+" "+s2)
  }

  def initLeftHandSide(t: Double, dt: Double, ratio: Double) {
    computeDiscount(ex, t, dt, discountVector)
    computeDrift(ex, t, dt, driftVector)
    computeVariance(ex, t, dt, driftVector, varianceVector)
    val multiplier = -ratio * dt
    val A = tridiagonal
    A.fill(0.0)
    A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier)
      .plusD1(1, x.length - 1, diffCache, driftVector, multiplier)
      .plusD0(1, x.length - 1, discountVector, -multiplier, 1.0)
    //    A.parabolicOperator(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier, driftVector, multiplier, 1.0 - multiplier * r)
  }

  private def initRightHandSideBoundaries(f: State) {
    //init rhs from lhs
    for (d <- 0 until rhs.stateDimensions) {
      val rhsValues = rhs.values(d)
      val fValues = f.values(d)
      val ratio = (1 - _a) / (r1)
      rhsValues(0) = (1 - (firstLine.value(0) - 1) * ratio) * fValues(0)
      for (i <- 1 until firstLine.iEnd) {
        rhsValues(0) -= firstLine.value(i) * ratio * fValues(i)
      }
      val m = x.length - 1
      rhsValues(m) = (1 - (lastLine.value(m) - 1) * ratio) * fValues(m)
      for (i <- lastLine.iStart until m) {
        rhsValues(m) -= ratio * lastLine.value(i) * fValues(i)
      }
    }
  }

  def initBoundaries(ratio: Double, t: Double, dt: Double, f: State) {
    val multiplier = -ratio * dt
    lowerBoundary.makeLine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), multiplier, firstLine)
    val m = x.length - 1
    upperBoundary.makeLine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), multiplier, lastLine)

    tridiagonal.setBoundaries(firstLine, lastLine)
  }

  //supposes that left hand side is initialized
  def initRightHandSide(f: State) {
    if (rhs == null) rhs = new State(f.stateDimensions, f.size)
    val ratio = (1 - _a) / r1
    for (d <- 0 until rhs.stateDimensions) {
      val rhsValues = rhs.values(d)
      val fValues = f.values(d)
      var j: Int = x.length - 2;
      while (j > 0) {
        val bStar = 1 - (tridiagonal.middle(j) - 1) * ratio;
        rhsValues(j) = -tridiagonal.lower(j) * ratio * fValues(j - 1) + bStar * fValues(j) - tridiagonal.upper(j) * ratio * fValues(j + 1);
        j -= 1
      }
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    initLeftHandSide(currentTime, dt, r1)
    initRightHandSide(f)
    initBoundaries(r1, currentTime, dt, f)
    initRightHandSideBoundaries(f)

    if (fTr == null) fTr = new State(f.stateDimensions, f.size)
    if (solverND == null) {
      solverND = new TridiagonalSolverND(solver, f.stateDimensions)
    }
    solverND.solve(tridiagonal, rhs.values, fTr.values)

    if (math.abs(r2 - r1) > Epsilon.MACHINE_EPSILON_SQRT) {
      initLeftHandSide(currentTime, dt, r2)
      initBoundaries(r2, currentTime, dt, f)
    }
    if (payoff != null) payoff.setTime(currentTime);
    solverND.solve(tridiagonal, fTr.values, f.values);
  }

}
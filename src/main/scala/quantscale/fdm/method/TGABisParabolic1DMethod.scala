package quantscale.fdm.method

;

import org.slf4j.LoggerFactory
import quantscale.fdm.TridiagonalMatrix
import java.util.Arrays
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State
import quantscale.fdm.OperatorLine
import quantscale.fdm.DifferentialCache
import quantscale.fdm.TridiagonalSolverND

/**
 * a must be between 1/2 and 2-sqrt(2) or over 2+sqrt(2)
 */
class TGABisParabolic1DMethod(payoff: FDPayoff, private var _a: Double = 0.54) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private var tridiagonal1, tridiagonal2: TridiagonalMatrix = null
  private var rhs1, rhs2: State = null

  private var specialIndex: Int = 0;
  private var r1, r2, s1, s2: Double = 0.0;
  private var _spec: Parabolic1DFDSpec = null;
  private var x: Array[Double] = null;
  var ex: Array[Double] = null;
  var driftVector: Array[Double] = null
  var varianceVector: Array[Double] = null
  var discountVector: Array[Double] = null
  var diffCache: DifferentialCache = null
  private var firstLine1, firstLine2: OperatorLine = null
  private var lastLine1, lastLine2: OperatorLine = null
  private var solverND: TridiagonalSolverND = null
  private var A: TridiagonalMatrix = null

  def copy(): Parabolic1DMethod = {
    return new TGAParabolic1DMethod(payoff, _a)
  }

  override def spec = _spec

  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV
    A = new TridiagonalMatrix(spec.grid.spaceVector.size)

    x = spec.grid.spaceVector;
    driftVector = new Array[Double](x.length)
    discountVector = new Array[Double](x.length)
    varianceVector = new Array[Double](x.length)
    ex = spec.grid.spaceTransform.transform(x);
    diffCache = new DifferentialCache(x)
    firstLine1 = new OperatorLine(0, 3)
    lastLine1 = new OperatorLine(x.length - 3, x.length)
    firstLine2 = new OperatorLine(0, 3)
    lastLine2 = new OperatorLine(x.length - 3, x.length)
    initFactors()
  }

  private def initFactors() {
    val sqrt = math.sqrt(_a * _a - 4 * _a + 2)
    r1 = (2 * _a - 1) / (_a + sqrt)
    r2 = (2 * _a - 1) / (_a - sqrt)
    s1 = (1 - _a + r1) / (r1 - r2)
    s2 = (1 - _a + r2) / (r2 - r1)
    //    println(_a+" "+sqrt+" "+r1+" "+r2+" "+s1+" "+s2)
  }

  def a_=(Value: Double): Unit = {
    _a = Value
    initFactors()
  }

  def initLeftHandSide(t: Double, dt: Double) {
    computeDiscount(ex, t, dt, discountVector)
    computeDrift(ex, t, dt, driftVector)
    computeVariance(ex, t, dt, driftVector, varianceVector)

    A.fill(0.0)
    A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5)
      .plusD1(1, x.length - 1, diffCache, driftVector, 1.0)
      .plusD0(1, x.length - 1, discountVector, -1.0)

  }

  def initBoundaries(t: Double, dt: Double, f: State) {
    //other ode factory without +1 and mult
    lowerBoundary.makeODELine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), firstLine1)
    val m = x.length - 1
    upperBoundary.makeODELine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), lastLine1)

    A.setBoundaries(firstLine1, lastLine1)

  }

  def initRightHandSide(f: State) {
    if (rhs1 == null) {

      rhs1 = new State(f.stateDimensions, f.size)
      rhs2 = new State(f.stateDimensions, f.size)
    }
    for (d <- 0 until rhs1.stateDimensions) {
      val rhs1Values = rhs1.values(d)
      val rhs2Values = rhs2.values(d)
      val fValues = f.values(d)
      System.arraycopy(fValues, 0, rhs1Values, 0, x.length)
      System.arraycopy(fValues, 0, rhs2Values, 0, x.length)

    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    initLeftHandSide(currentTime, dt)
    initRightHandSide(f)
    initBoundaries(currentTime, dt, f)
    //boundaries will be updated here as well
    tridiagonal1 = TridiagonalMatrix.identity(A.size).plus(A, -r1 * dt)
    tridiagonal2 = TridiagonalMatrix.identity(A.size).plus(A, -r2 * dt)

    if (solverND == null) {
      solverND = new TridiagonalSolverND(solver, f.stateDimensions)
    }
    solverND.solve(tridiagonal1, rhs1.values, f.values)
    solverND.solve(tridiagonal2, rhs2.values, rhs2.values);

    for (d <- 0 until f.stateDimensions) {
      val fValues = f.values(d)
      val rhs2Values = rhs2.values(d)
      var i = 0
      while (i < fValues.length) {
        fValues(i) = s1 * fValues(i) + s2 * rhs2Values(i)
        i += 1
      }
    }
  }

}
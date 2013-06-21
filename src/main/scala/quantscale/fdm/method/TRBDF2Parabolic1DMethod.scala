package quantscale.fdm.method
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.OperatorLine
import quantscale.fdm.DifferentialCache
import quantscale.fdm.State
import quantscale.fdm.TridiagonalSolverND

/**
 * Solves the log transformed problem: the space contains log(assetPrice)
 */
class TRBDF2Parabolic1DMethod(payoff: FDPayoff = null) extends Parabolic1DMethod {

  final val theta = 0.5;
  final val alpha = 2 - Math.sqrt(2);
  final val backcoeff = (1 - theta) * alpha;

  var tridiagonal: TridiagonalMatrix = null;
  private var rhs: State = null;
  private var _spec: Parabolic1DFDSpec = null;
  private var fTr: State = null;
  private var x: Array[Double] = null;
  var ex: Array[Double] = null;
  private val frac = 1 / (alpha * (2 - alpha));

  var driftVector: Array[Double] = null
  var varianceVector: Array[Double] = null
  var discountVector: Array[Double] = null
  var diffCache: DifferentialCache = null
  private var firstLine: OperatorLine = null
  private var lastLine: OperatorLine = null
  private var solverND: TridiagonalSolverND = null
  
  override def spec = _spec

  def copy(): Parabolic1DMethod = {
    return new TRBDF2Parabolic1DMethod(payoff)
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
  }

  def initLeftHandSide(t: Double, dt: Double) {
    computeDiscount(ex, t, dt, discountVector)
    computeDrift(ex, t, dt, driftVector)
    computeVariance(ex, t, dt, driftVector, varianceVector)
    val multiplier = -backcoeff * dt
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
      val ratio = theta / (1 - theta)
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

  def initBoundaries(t: Double, dt: Double, f: State) {
    val multiplier = -backcoeff * dt
    lowerBoundary.makeLine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), multiplier, firstLine)
    val m = x.length - 1
    upperBoundary.makeLine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), multiplier, lastLine)

    initRightHandSideBoundaries(f)
    tridiagonal.setBoundaries(firstLine, lastLine)
  }

  //supposes that left hand side is initialized
  def initRightHandSide(f: State) {
    if (rhs == null) rhs = new State(f.stateDimensions, f.size)
    val ratio = theta / (1 - theta)
    for (d <- 0 until rhs.stateDimensions) {
      val rhsValues = rhs.values(d)
      val fValues = f.values(d)
      var j: Int = x.length - 2;
      while (j > 0) {
        val bStar = 1 - (tridiagonal.middle(j) - 1) * ratio;
        rhsValues(j) = -tridiagonal.lower(j) * fValues(j - 1) + bStar * fValues(j) - tridiagonal.upper(j) * fValues(j + 1);
        j -= 1
      }
    }
  }

  def computeBDF2Rhs(price: State, prevPrice: State) {
    for (d <- 0 until rhs.stateDimensions) {
      val rhsValues = rhs.values(d)
      val priceValues = price.values(d)
      val prevPriceValues = prevPrice.values(d)
      var i: Int = spec.grid.spaceSize - 1;
      while (i >= 0) {
        rhsValues(i) = frac * (priceValues(i) - (1 - alpha) * (1 - alpha) * prevPriceValues(i));
        i -= 1;
      }
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    initLeftHandSide(currentTime, dt)
    initRightHandSide(f)
    initBoundaries(currentTime, dt, f)
    if (payoff != null) payoff.setTime(currentTime + dt * alpha)
    if (fTr == null) fTr = new State(f.stateDimensions, f.size)
    if (solverND == null) {
      solverND = new TridiagonalSolverND(solver, f.stateDimensions)
    }
    solverND.solve(tridiagonal, rhs.values, fTr.values)
    computeBDF2Rhs(fTr, f);
    if (payoff != null) payoff.setTime(currentTime);
    solverND.solve(tridiagonal, rhs.values, f.values);
  }

}
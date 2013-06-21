package quantscale.fdm.method

import java.util.Arrays

import org.slf4j.LoggerFactory

import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.DifferentialCache
import quantscale.fdm.OperatorLine
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.TridiagonalSolverND

object ThetaParabolic1DMethod {
  val THETA_CRANK_NICOLSON = 0.5
  val THETA_EXPLICIT = 1.0
  val THETA_IMPLICIT = 0.0 
}
class ThetaParabolic1DMethod(var theta: Double = ThetaParabolic1DMethod.THETA_CRANK_NICOLSON) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  var tridiagonal: TridiagonalMatrix = null;
  var rhs: State = null;
  private var _spec: Parabolic1DFDSpec = null;

  private var x: Array[Double] = null
  private var ex: Array[Double] = null

  private var driftVector: Array[Double] = null
  private var varianceVector: Array[Double] = null
  private var discountVector: Array[Double] = null
  private var diffCache: DifferentialCache = null
  private var firstLine: OperatorLine = null
  private var lastLine: OperatorLine = null
  private var solverND : TridiagonalSolverND = null

  override def spec = _spec
  
  override def copy() : Parabolic1DMethod = {
    val c = new ThetaParabolic1DMethod(theta)
    c.lowerBoundary = lowerBoundary
    c.upperBoundary = upperBoundary
    c.smearingReducer = smearingReducer
    return c
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
    if (theta == ThetaParabolic1DMethod.THETA_EXPLICIT) {
       val multiplier = dt
        val A = tridiagonal
        A.fill(0.0)
        A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier)
          .plusD1(1, x.length - 1, diffCache, driftVector, multiplier)
          .plusD0(1, x.length - 1, discountVector, -multiplier, 1.0)
    } else {
      val multiplier = -(1 - theta) * dt 
      val A = tridiagonal
      A.fill(0.0)
      A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier)
        .plusD1(1, x.length - 1, diffCache, driftVector, multiplier)
        .plusD0(1, x.length - 1, discountVector, -multiplier, 1.0)
    }
    //    A.parabolicOperator(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier, driftVector, multiplier, 1.0 - multiplier * r)

  }

  private def initRightHandSideBoundaries(f: State) {
    //init rhs from lhs
    if (theta != ThetaParabolic1DMethod.THETA_IMPLICIT && theta != ThetaParabolic1DMethod.THETA_EXPLICIT) {
      //otherwise already handled in copy from f to rhs
      val ratio = theta / (1 - theta)
      for (d <- 0 until rhs.stateDimensions) {
        val rhsValues = rhs.values(d)
        val fValues = f.values(d)
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
  }

  def initBoundaries(t: Double, dt: Double, f: State) {
    val multiplier = if (theta == ThetaParabolic1DMethod.THETA_EXPLICIT) dt else -(1 - theta) * dt 
    lowerBoundary.makeLine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), multiplier, firstLine)
    val m = x.length - 1
    upperBoundary.makeLine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), multiplier, lastLine)

    initRightHandSideBoundaries(f)
    tridiagonal.setBoundaries(firstLine, lastLine)
  }

  //supposes that left hand side is initialized
  def initRightHandSide(f: State, alpha: Double=1.0) {
    if (rhs == null) rhs = new State(f.stateDimensions, f.size)
    for (d <- 0 until rhs.stateDimensions) {
      val rhsValues = rhs.values(d)
      val fValues = f.values(d)
      if (theta == ThetaParabolic1DMethod.THETA_IMPLICIT) {
        if (alpha == 1.0) {
        System.arraycopy(fValues, 0, rhsValues, 0, x.length)
        } else {
          var j = x.length-1
          while (j>=0) {
            rhsValues(j) = alpha*fValues(j)
            j-=1
          }
        }
      } else if (theta == ThetaParabolic1DMethod.THETA_EXPLICIT) {
        System.arraycopy(fValues, 0, rhsValues, 0, x.length)
      } else {
        val ratio = theta / (1 - theta)
        var j: Int = x.length - 2;
        while (j > 0) {
          val bStar = 1 - (tridiagonal.middle(j) - 1) * ratio;
          rhsValues(j) = -tridiagonal.lower(j) * fValues(j - 1) + bStar * fValues(j) - tridiagonal.upper(j) * fValues(j + 1);
          j -= 1
        }
      }
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    initLeftHandSide(currentTime, dt)
    initRightHandSide(f)
    initBoundaries(currentTime, dt, f)
    if (logger.isDebugEnabled()) {
      logger.debug("a=" + Arrays.toString(tridiagonal.lower));
      logger.debug("b=" + Arrays.toString(tridiagonal.middle));
      logger.debug("c=" + Arrays.toString(tridiagonal.upper));
      logger.debug("d=" + rhs);
    }
    if (theta == ThetaParabolic1DMethod.THETA_EXPLICIT) {
      for (d <- 0 until f.stateDimensions) {
        tridiagonal.multiply(rhs.values(d),f.values(d))
      }
    } else {
      if (solverND == null) {
        solverND = new TridiagonalSolverND(solver, f.stateDimensions)
      }
      solverND.solve(tridiagonal, rhs.values, f.values)
//      for (d<-0 until f.stateDimensions) {
//        solver.solve(tridiagonal, rhs.values(d), f.values(d))
//      }
    }
  }
}
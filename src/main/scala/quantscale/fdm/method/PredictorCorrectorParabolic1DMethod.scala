package quantscale.fdm.method

import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Parabolic1DFDSpec
import org.slf4j.LoggerFactory
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.State

class PredictorCorrectorParabolic1DMethod(val payoff: FDPayoff) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private val thetaMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_EXPLICIT);
  private var tridiagonalHalf: TridiagonalMatrix = null
  private var fFull, fTemp, fCorrector: State = null

  private var specialIndex: Int = 0;

  def copy(): Parabolic1DMethod = {
    return new PredictorCorrectorParabolic1DMethod(payoff)
  }

  override def spec = thetaMethod.spec

  override def initSystem(specV: Parabolic1DFDSpec) {
    thetaMethod.initSystem(specV);
    thetaMethod.solver = solver;
    thetaMethod.smearingReducer = smearingReducer
    thetaMethod.lowerBoundary = lowerBoundary
    thetaMethod.upperBoundary = upperBoundary
    tridiagonalHalf = new TridiagonalMatrix(thetaMethod.tridiagonal.size)
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    if (fFull == null) {
      fFull = new State(f.stateDimensions, f.size)
      fTemp = new State(f.stateDimensions, f.size)
      fCorrector = new State(f.stateDimensions, f.size)
    }
    thetaMethod.initLeftHandSide(currentTime, dt)
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime, dt, f)
    for (d <- 0 until f.stateDimensions) {
      thetaMethod.tridiagonal.multiply(thetaMethod.rhs.values(d), fTemp.values(d))
    }
    val m = thetaMethod.tridiagonal.size
    var i = 0
    while (i < m) {
      tridiagonalHalf.lower(i) = thetaMethod.tridiagonal.lower(i) * 0.5
      tridiagonalHalf.upper(i) = thetaMethod.tridiagonal.upper(i) * 0.5
      tridiagonalHalf.middle(i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) * 0.5
      i += 1
    }
    thetaMethod.tridiagonal = tridiagonalHalf
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime, dt * 0.5, f)
    for (d <- 0 until f.stateDimensions) {
      tridiagonalHalf.multiply(thetaMethod.rhs.values(d), fFull.values(d))
    }

    thetaMethod.initRightHandSide(fTemp)
    thetaMethod.initBoundaries(currentTime, dt * 0.5, f)
    for (d <- 0 until f.stateDimensions) {
      tridiagonalHalf.multiply(thetaMethod.rhs.values(d), fCorrector.values(d))
      val fValues = f.values(d)
      val fTempValues = fTemp.values(d)
      val fFullValues = fFull.values(d)
      val fCorrectorValues = fCorrector.values(d)
      i = 0
      while (i < m) {
        fValues(i) = fFullValues(i) + (fCorrectorValues(i) - fTempValues(i))
        i += 1
      }
    }
  }

}
package quantscale.fdm.method

import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.Parabolic1DFDSpec
import org.slf4j.LoggerFactory
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.State

class RK2Parabolic1DMethod(val payoff: FDPayoff) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private val thetaMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_EXPLICIT);
  private var tridiagonalHalf: TridiagonalMatrix = null
  private var fFull, fTemp: State = null

  private var specialIndex: Int = 0;

  override def spec = thetaMethod.spec

  def copy(): Parabolic1DMethod = {
    return new RK2Parabolic1DMethod(payoff)
  }

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
    }
    thetaMethod.initLeftHandSide(currentTime, dt)
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime, dt, f)

    val m = thetaMethod.tridiagonal.size
    var i = 0
    while (i < m) {
      tridiagonalHalf.lower(i) = thetaMethod.tridiagonal.lower(i) * 0.5
      tridiagonalHalf.upper(i) = thetaMethod.tridiagonal.upper(i) * 0.5
      tridiagonalHalf.middle(i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) * 0.5
      i += 1
    }

    val tridiagonal = thetaMethod.tridiagonal
    thetaMethod.tridiagonal = tridiagonalHalf
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f)
    if (payoff != null) payoff.setTime(currentTime + dt * 0.5);
    for (d <- 0 until f.stateDimensions) {
      tridiagonalHalf.multiply(thetaMethod.rhs.values(d), fTemp.values(d))
    }

    //    while (i < m) {
    //      tridiagonal.middle(i) = tridiagonal.middle(i) - 1
    //      i += 1
    //    }
    thetaMethod.tridiagonal = tridiagonal
    thetaMethod.initBoundaries(currentTime, dt, fTemp)
    if (payoff != null) payoff.setTime(currentTime);
    for (d <- 0 until f.stateDimensions) {
      val fValues = f.values(d)
      val fTempValues = fTemp.values(d)
      val fFullValues = fFull.values(d)
      tridiagonal.multiply(fTempValues, fFullValues)
      i = 0
      while (i < m) {
        fValues(i) = fValues(i) - fTempValues(i) + fFullValues(i)
        i += 1
      }
    }
  }

}
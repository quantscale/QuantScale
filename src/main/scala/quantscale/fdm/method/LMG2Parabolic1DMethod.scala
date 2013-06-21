package quantscale.fdm.method;

import org.slf4j.LoggerFactory
import quantscale.fdm.TridiagonalMatrix
import java.util.Arrays
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State
import quantscale.fdm.TridiagonalSolverND

class LMG2Parabolic1DMethod(payoff: FDPayoff) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private val thetaMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
  private var tridiagonalHalf: TridiagonalMatrix = null
  private var fFull, fTemp: State = null
  private var solverND: TridiagonalSolverND = null

  private var specialIndex: Int = 0;

  def copy(): Parabolic1DMethod = {
    return new LMG2Parabolic1DMethod(payoff)
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
    }
    thetaMethod.initLeftHandSide(currentTime, dt)
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime, dt, f)
    if (solverND == null) {
      solverND = new TridiagonalSolverND(solver, f.stateDimensions)
    }
    
    
    //elliot ockendon does not preserve rhs
    solverND.solve(thetaMethod.tridiagonal, thetaMethod.rhs.values, fFull.values)
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
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f)
    if (payoff != null) payoff.setTime(currentTime + dt * 0.5);
    solverND.solve(tridiagonalHalf, thetaMethod.rhs.values, fTemp.values)
    //    thetaMethod.initRightHandSide(fTemp)
    thetaMethod.initBoundaries(currentTime, dt * 0.5, fTemp)
    if (payoff != null) payoff.setTime(currentTime);
    solverND.solve(tridiagonalHalf, fTemp.values, f.values)

    for (d <- 0 until f.stateDimensions) {
      val fValues = f.values(d)
      val fFullValues = fFull.values(d)
      i = 0
      while (i < m) {
        fValues(i) = 2 * fValues(i) - fFullValues(i)
        i += 1
      }
    }
  }

}
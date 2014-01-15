package quantscale.fdm.method

;

import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import org.slf4j.LoggerFactory
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.TridiagonalSolverND

class LRE3Parabolic1DMethod(payoff: FDPayoff) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private val thetaMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
  private var tridiagonalHalf, tridiagonalThird: TridiagonalMatrix = null
  private var f1, f2, f2tmp, f3tmp: State = null

  private var specialIndex: Int = 0;

  override def spec = thetaMethod.spec

  private var solverND: TridiagonalSolverND = null

  def copy(): Parabolic1DMethod = {
    return new LRE3Parabolic1DMethod(payoff)
  }

  override def initSystem(specV: Parabolic1DFDSpec) {
    thetaMethod.initSystem(specV)
    thetaMethod.solver = solver
    thetaMethod.smearingReducer = smearingReducer
    thetaMethod.lowerBoundary = lowerBoundary
    thetaMethod.upperBoundary = upperBoundary
    tridiagonalHalf = new TridiagonalMatrix(thetaMethod.tridiagonal.size)
    tridiagonalThird = new TridiagonalMatrix(thetaMethod.tridiagonal.size)
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    if (f1 == null) {
      f1 = new State(f.stateDimensions, f.size)
      f2 = new State(f.stateDimensions, f.size)
      f2tmp = new State(f.stateDimensions, f.size)
      f3tmp = new State(f.stateDimensions, f.size)
    }
    thetaMethod.initLeftHandSide(currentTime, dt)
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime, dt, f)
    //TODO create an ND solver => work in parallel
    if (solverND == null) {
      solverND = new TridiagonalSolverND(solver, f.stateDimensions)
    }

    //clone f?
    solverND.solve(thetaMethod.tridiagonal, thetaMethod.rhs.values, f1.values)

    val m = thetaMethod.tridiagonal.size
    var i = 0
    while (i < m) {
      tridiagonalHalf.lower(i) = thetaMethod.tridiagonal.lower(i) / 2
      tridiagonalHalf.upper(i) = thetaMethod.tridiagonal.upper(i) / 2
      tridiagonalHalf middle (i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) / 2
      tridiagonalThird.lower(i) = thetaMethod.tridiagonal.lower(i) / 3
      tridiagonalThird.upper(i) = thetaMethod.tridiagonal.upper(i) / 3
      tridiagonalThird middle (i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) / 3
      i += 1
    }
    thetaMethod.tridiagonal = tridiagonalHalf
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f)
    if (payoff != null) payoff.setTime(currentTime + dt / 2);
    //clone f?
    solverND.solve(tridiagonalHalf, thetaMethod.rhs.values, f2tmp.values)
    thetaMethod.tridiagonal = tridiagonalHalf
    thetaMethod.initRightHandSide(f2tmp)
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f2tmp)
    if (payoff != null) payoff.setTime(currentTime);
    solverND.solve(tridiagonalHalf, thetaMethod.rhs.values, f2.values)

    thetaMethod.tridiagonal = tridiagonalThird
    thetaMethod.initRightHandSide(f)
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f)
    if (payoff != null) payoff.setTime(currentTime + 2 * dt / 3);

    solverND.solve(tridiagonalThird, thetaMethod.rhs.values, f2tmp.values)

    thetaMethod.initRightHandSide(f2tmp)
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f2tmp)
    if (payoff != null) payoff.setTime(currentTime + dt / 3);
    solverND.solve(tridiagonalThird, thetaMethod.rhs.values, f3tmp.values)
    if (payoff != null) payoff.setTime(currentTime);

    thetaMethod.initRightHandSide(f3tmp)
    thetaMethod.initBoundaries(currentTime + dt * 0.5, dt * 0.5, f3tmp)
    solverND.solve(tridiagonalThird, thetaMethod.rhs.values, f.values)
    for (d <- 0 until f.stateDimensions) yield {
      val fv = f.values(d)
      val f2v = f2.values(d)
      val f1v = f1.values(d)
      i = 0
      while (i < m) {
        fv(i) = 4.5 * fv(i) - 4 * f2v(i) + 0.5 * f1v(i);
        i += 1
      }
    }
  }

}
package quantscale.fdm.method

import quantscale.fdm.Epsilon
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.State

/**
 * specialPoints must be sorted
 */
class RannacherCentralBSM1FMethod(specialPoints: Array[Double], payoff: FDPayoff = null) extends Parabolic1DMethod {

  val thetaMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
  private var tridiagonalRan: TridiagonalMatrix = null
  private var fTmp: State = null
  private var specialIndex: Int = 0;

  override def spec = thetaMethod.spec

  def copy(): Parabolic1DMethod = {
    val c = new RannacherCentralBSM1FMethod(specialPoints, payoff)
    c.solver = solver.copy()
    return c
  }
  
  override def initSystem(specV: Parabolic1DFDSpec) {
    thetaMethod.initSystem(specV);
    thetaMethod.solver = solver;
    thetaMethod.smearingReducer = smearingReducer
    thetaMethod.lowerBoundary = lowerBoundary
    thetaMethod.upperBoundary = upperBoundary
    tridiagonalRan = new TridiagonalMatrix(thetaMethod.tridiagonal.size)
    specialIndex = specialPoints.length - 1;
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    if (isSpecialPoint(currentTime)) {

      thetaMethod.theta = ThetaParabolic1DMethod.THETA_IMPLICIT;
      thetaMethod.initLeftHandSide(currentTime, dt)
      //      thetaMethod.initRightHandSide(f)
      thetaMethod.initBoundaries(currentTime, dt, f)
      val m = thetaMethod.tridiagonal.size
      var i = 0
      while (i < m) {
        tridiagonalRan.lower(i) = thetaMethod.tridiagonal.lower(i) / 4
        tridiagonalRan.upper(i) = thetaMethod.tridiagonal.upper(i) / 4
        tridiagonalRan.middle(i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) / 4
        i += 1
      }
      //4 1/4 steps implicit euler (M. Giles)
      thetaMethod.tridiagonal = tridiagonalRan
      //      thetaMethod.initRightHandSide(f)
      thetaMethod.initBoundaries(currentTime + 3 * dt / 4, dt * 0.25, f)
      if (payoff != null) payoff.setTime(currentTime + 3 * dt / 4)
      if (fTmp == null) fTmp = new State(f.stateDimensions, f.size)
      for (d <- 0 until f.stateDimensions) {
        solver.solve(tridiagonalRan, f.values(d), fTmp.values(d))
      }

      //      thetaMethod.initRightHandSide(fTmp)
      thetaMethod.initBoundaries(currentTime + 2 * dt / 4, dt * 0.25, fTmp)
      if (payoff != null) payoff.setTime(currentTime + 2 * dt / 4)
      for (d <- 0 until f.stateDimensions) {
        solver.solve(tridiagonalRan, fTmp.values(d), f.values(d))
      }

      //      thetaMethod.initRightHandSide(f)
      thetaMethod.initBoundaries(currentTime + dt / 4, dt * 0.25, f)
      if (payoff != null) payoff.setTime(currentTime + dt / 4)
      for (d <- 0 until f.stateDimensions) {
        solver.solve(tridiagonalRan, f.values(d), fTmp.values(d))
      }
      //      thetaMethod.initRightHandSide(fTmp)
      thetaMethod.initBoundaries(currentTime, dt * 0.25, fTmp)
      if (payoff != null) payoff.setTime(currentTime)
      for (d <- 0 until f.stateDimensions) {
        solver.solve(tridiagonalRan, fTmp.values(d), f.values(d))
      }
      thetaMethod.theta = ThetaParabolic1DMethod.THETA_CRANK_NICOLSON;
    } else {
      thetaMethod.solve(currentTime, dt, f);
    }
  }

  def isSpecialPoint(t: Double): Boolean = {
    if (specialIndex >= 0 && t <= specialPoints(specialIndex) + Epsilon.MACHINE_EPSILON_SQRT) {
      specialIndex -= 1;
      return true;
    }
    if (payoff != null) {
      return payoff.isDiscontinuous
    }
    return false;
  }

}
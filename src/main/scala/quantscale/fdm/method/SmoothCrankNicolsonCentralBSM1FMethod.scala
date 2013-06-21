package quantscale.fdm.method
import quantscale.fdm.Epsilon
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State
import quantscale.fdm.TridiagonalSolverND

/**
 * A better fix for Crank-Nicolson - LMG smoothing
 *
 * Crank-Nicolson is known to be unreliable around discontinuities. Rannacher proposed a fix
 * using 2 1/2 implicit Euler smoothing steps after each discontinuity. Giles showed that in practice,
 * 4 1/4 implicit Euler steps is optimal (is it?).
 * We consider here a smoothing based 2 steps of Lawson-Morris-Gourlay scheme in order to maintain
 * the second order accuracy, even when the number of discontinuities is high.
 *
 * If discontinuities are continuous, it will revert to a classic LMG scheme with 1/2 time step. This is
 * typically the case for an American option.
 *
 * This is very good for the daily barrier kind of problems: it saves time over LMG when the number of time
 * steps is high enough, and it provides higher accuracy as well.
 *
 *
 * specialPoints must be sorted
 */

class SmoothCrankNicolsonCentralBSM1FMethod(specialPoints: Array[Double], var payoff: FDPayoff = null) extends Parabolic1DMethod {

  val thetaMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON);

  private var f1, f2tmp: State = null
  private var tridiagonalHalf, tridiagonalQuarter: TridiagonalMatrix = null

  private var specialIndex: Int = 0;
  private var solverND: TridiagonalSolverND = null
  override def spec = thetaMethod.spec
    def copy() : Parabolic1DMethod = {
    return new SmoothCrankNicolsonCentralBSM1FMethod(specialPoints, payoff)
  }
  override def initSystem(specV: Parabolic1DFDSpec) {
    thetaMethod.initSystem(specV);
    thetaMethod.solver = solver;
    thetaMethod.smearingReducer = smearingReducer
    thetaMethod.lowerBoundary = lowerBoundary
    thetaMethod.upperBoundary = upperBoundary

    tridiagonalHalf = new TridiagonalMatrix(thetaMethod.tridiagonal.size)
    tridiagonalQuarter = new TridiagonalMatrix(thetaMethod.tridiagonal.size)
    specialIndex = specialPoints.length - 1;
    solverND = null
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    if (isSpecialPoint(currentTime)) {
      if (f1 == null) {
        f1 = new State(f.stateDimensions, f.size)
        f2tmp = new State(f.stateDimensions, f.size)
      }
      if (solverND == null) {
        solverND = new TridiagonalSolverND(solver, f.stateDimensions)
      }

      thetaMethod.theta = ThetaParabolic1DMethod.THETA_IMPLICIT;
      thetaMethod.initLeftHandSide(currentTime, dt) //FIXME use dt here and divide later by 2, don't /4 everytime
      thetaMethod.initRightHandSide(f)
      thetaMethod.initBoundaries(currentTime, dt, f)
      /* 1 full step becomes 2* 2 extrapolated 1/4 steps */
      val m = thetaMethod.tridiagonal.size
      var i = 0
      var frac = 0.5
      while (i < m) {
        tridiagonalHalf.lower(i) = thetaMethod.tridiagonal.lower(i) * frac
        tridiagonalHalf.upper(i) = thetaMethod.tridiagonal.upper(i) * frac
        tridiagonalHalf.middle(i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) * frac
        i += 1
      }
      frac = 0.25
      i = 0
      while (i < m) {
        tridiagonalQuarter.lower(i) = thetaMethod.tridiagonal.lower(i) * frac
        tridiagonalQuarter.upper(i) = thetaMethod.tridiagonal.upper(i) * frac
        tridiagonalQuarter.middle(i) = 1 + (thetaMethod.tridiagonal.middle(i) - 1) * frac
        i += 1
      }

      frac = 0.5
      thetaMethod.tridiagonal = tridiagonalHalf
      var fCopy = f.copy()
      thetaMethod.initBoundaries(currentTime + dt * (1 - frac), dt * frac, fCopy)
      if (payoff != null) payoff.setTime(currentTime + dt * (1 - frac))
      solverND.solve(tridiagonalHalf, fCopy.values, f1.values)

      frac = 0.25
      thetaMethod.tridiagonal = tridiagonalQuarter
      fCopy = f.copy()
      thetaMethod.initBoundaries(currentTime + dt * (1 - frac), dt * frac, fCopy)
      if (payoff != null) payoff.setTime(currentTime + dt / 2 + dt / 4);
      solverND.solve(tridiagonalQuarter, fCopy.values, f2tmp.values)

      thetaMethod.initBoundaries(currentTime + dt * (1 - 2 * frac), dt * frac, f2tmp)
      if (payoff != null) payoff.setTime(currentTime + dt * (1 - 2 * frac));
      solverND.solve(tridiagonalQuarter, f2tmp.values, f.values)
      for (d <- 0 until f.stateDimensions) {
        val fValues = f.values(d)
        val f1Values = f1.values(d)
        i = 0
        while (i < m) {
          fValues(i) = 2 * fValues(i) - f1Values(i);
          i += 1
        }
      }

      frac = 0.5
      fCopy = f.copy()
      thetaMethod.tridiagonal = tridiagonalHalf
      thetaMethod.initBoundaries(currentTime, dt * frac, fCopy)
      if (payoff != null) payoff.setTime(currentTime);
      solverND.solve(tridiagonalHalf, fCopy.values, f1.values)

      frac = 0.25
      fCopy = f.copy
      thetaMethod.tridiagonal = tridiagonalQuarter
      thetaMethod.initBoundaries(currentTime + dt * frac, dt * frac, fCopy)
      if (payoff != null) payoff.setTime(currentTime + dt / 4);
      solverND.solve(tridiagonalQuarter, fCopy.values, f2tmp.values)

      thetaMethod.initBoundaries(currentTime, dt * frac, f2tmp)
      if (payoff != null) payoff.setTime(currentTime);
      solverND.solve(tridiagonalQuarter, f2tmp.values, f.values)
      for (d <- 0 until f.stateDimensions) {
        val fValues = f.values(d)
        val f1Values = f1.values(d)
        i = 0
        while (i < m) {
          fValues(i) = 2 * fValues(i) - f1Values(i);
          i += 1
        }

      }

      thetaMethod.theta = ThetaParabolic1DMethod.THETA_CRANK_NICOLSON;
    } else {
      thetaMethod.solve(currentTime, dt, f);
    }
  }

  def isSpecialPoint(t: Double): Boolean = {

    if (specialIndex >= 0 && t <= specialPoints(specialIndex) + Epsilon.MACHINE_EPSILON_SQRT) {
      specialIndex -= 1
      return true;
    }
    if (payoff != null) {
      return payoff.isDiscontinuous
    }
    return false
  }

}
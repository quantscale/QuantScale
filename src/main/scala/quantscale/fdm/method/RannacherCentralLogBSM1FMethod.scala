package quantscale.fdm.method

import quantscale.fdm.BSM1FFDSpec
import quantscale.fdm.Epsilon

/**
 * specialPoints must be sorted
 */
class RannacherCentralLogBSM1FMethod(specialPoints: Array[Double]) extends BSM1FMethod {

  val thetaMethod = new ThetaCentralLogBSM1FMethod();

  private var specialIndex: Int = 0;

  override def initSystem(specV: BSM1FFDSpec) {
    thetaMethod.initSystem(specV);
    thetaMethod.solver = solver;
    specialIndex = specialPoints.length - 1;
  }

  override def solve(currentTime: Double, dt: Double, f: Array[Double]) {
    if (isSpecialPoint(currentTime)) {
      thetaMethod.theta = ThetaParabolic1DMethod.THETA_IMPLICIT;
      thetaMethod.solve(currentTime, dt / 2.0, f);
      thetaMethod.solve(currentTime + dt / 2.0, dt / 2.0, f);
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
    return false;
  }

}
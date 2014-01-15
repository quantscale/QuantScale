package quantscale.fdm.method

import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State

class RannacherParabolic1DMethod extends ThetaParabolic1DMethod() {
  private var step = 0

  override def initSystem(specV: Parabolic1DFDSpec) {
    super.initSystem(specV)
    step = 0
  }

  override def solve(currentTime: Double, dt: Double, f: State) {
    if (step == 0) {
      theta = ThetaParabolic1DMethod.THETA_IMPLICIT
      super.solve(currentTime + dt / 2, dt / 2, f)
      super.solve(currentTime, dt / 2, f)
    } else {
      theta = ThetaParabolic1DMethod.THETA_CRANK_NICOLSON
      super.solve(currentTime, dt, f)
    }
  }
}

package quantscale.fdm.method
import quantscale.fdm.BSM1FFDSpec
import quantscale.fdm.TridiagonalSolver

abstract class BSM1FMethod {
  def initSystem(specV: BSM1FFDSpec);

  def solve(currentTime: Double, dt: Double, f: Array[Double]);

  var solver : TridiagonalSolver = null;
}
package quantscale.math
import quantscale.fdm.Epsilon

object InverseHyperbolics {
  def asinh(x: Double): Double = {
    if (math.abs(x) < Epsilon.MACHINE_EPSILON_FOURTH_ROOT) {
      return x - x * x * x / 6;
    }
    return math.log(x + math.sqrt(1 + x * x));
  }
}
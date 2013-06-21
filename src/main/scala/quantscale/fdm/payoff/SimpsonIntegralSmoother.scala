package quantscale.fdm.payoff
import quantscale.fdm.Epsilon
import quantscale.fdm.mesh.MeshUtil

class SimpsonIntegralSmoother(discontinuities: Array[Double]) extends FDPayoffSmoother {

  var h = 0;
  final val INTEGRATION_STEPS = 10;

  private var v1 = new Array[Double](INTEGRATION_STEPS);
  private var v2 = new Array[Double](INTEGRATION_STEPS);
  private var v3 = new Array[Double](INTEGRATION_STEPS);

  def findNearestIndex(v: Array[Double], z: Double): Int = {
    var i = MeshUtil.findIndex(v, z);
    return i;
  }

  def makeSmooth(payoff: FDPayoff) {
    val originalStateFull = payoff.state
    val originalState = payoff.state.price;
    val originalSpace = payoff.space;
    val ooSpace = payoff.originalSpace;

    for (discontinuity <- discontinuities) {
      val indexDiscontinuity = findNearestIndex(ooSpace, discontinuity);

      if (indexDiscontinuity > 0) {
        val low = (ooSpace(indexDiscontinuity) + ooSpace(indexDiscontinuity - 1)) / 2;
        val high = (ooSpace(indexDiscontinuity) + ooSpace(indexDiscontinuity + 1)) / 2;
        val averageLow = simpsonIntegration(payoff, low, ooSpace(indexDiscontinuity) - Epsilon.MACHINE_EPSILON_SQRT)
        val averageHigh = simpsonIntegration(payoff, ooSpace(indexDiscontinuity) + Epsilon.MACHINE_EPSILON_SQRT, high)
        originalState(indexDiscontinuity) = (averageLow + averageHigh) / (high - low);
      }
    }
    payoff.space = originalSpace;
    payoff.originalSpace = ooSpace;
    payoff.state = originalStateFull
    payoff.state.price = originalState;
  }

  private def simpsonIntegration(payoff: quantscale.fdm.payoff.FDPayoff, low: Double, high: Double): Double = {
    val dy = (high - low) / INTEGRATION_STEPS;
    for (k <- 0 until INTEGRATION_STEPS) {
      val y = low + k * dy;
      v1(k) = y;
      v2(k) = y + dy / 2;
      v3(k) = y + dy;
    }
    payoff.initState(v1);
    payoff.eval();
    v1 = payoff.state.price;
    payoff.initState(v2);
    payoff.eval();
    v2 = payoff.state.price;
    payoff.initState(v3);
    payoff.eval();
    v3 = payoff.state.price;
    // average
    var average = 0.0;
    //FIXME integrate to discont-eps and then from discont+eps to h
    for (k <- 0 until INTEGRATION_STEPS) {
      // simpson formula for integral approximation
      average += (v1(k) + 4 * v2(k) + v3(k)) * dy / 6;
    }
    average
  }
}
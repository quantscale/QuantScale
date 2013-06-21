package quantscale.fdm.payoff
import quantscale.fdm.Epsilon

class VanillaFDPayoff(
  isCall: Boolean,
  strike: Double,
  expiryTime: Double) extends FDPayoff {
  override def copy(): FDPayoff = {
    return new VanillaFDPayoff(isCall, strike, expiryTime)
  }

  val putCallSign: Int = if (isCall) 1 else -1;

  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous

  override def eval(): Unit = {
    if (math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      var i: Int = 0;
      while (i < state.size) {
        state.price(i) = math.max(putCallSign * (space(i) - strike), 0);
        i += 1;
      }
      _isDiscontinuous = true
    } else {
      _isDiscontinuous = false
    }
  }

}
package quantscale.fdm.payoff

import quantscale.fdm.Epsilon

class VanillaBondOptionFDPayoff(val isCall: Boolean, val strike: Double, val optionExpiryTime: Double, val bondMaturityTime: Double, val coupon: Double = 1.0) extends FDPayoff {
  override def copy(): FDPayoff = {
    return new VanillaBondOptionFDPayoff(isCall, strike, optionExpiryTime, bondMaturityTime, coupon)
  }
  val putCallSign: Int = if (isCall) 1 else -1;

  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous

  override def eval(): Unit = {
    if (math.abs(bondMaturityTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      for (i <- 0 until state.size) {
        state.price(i) = coupon
      }
    }
    if (math.abs(optionExpiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      for (i <- 0 until state.size) {
        state.price(i) = math.max(putCallSign * (state.price(i) - strike), 0);
      }
      _isDiscontinuous = true
    } else {
      _isDiscontinuous = false
    }
  }
}


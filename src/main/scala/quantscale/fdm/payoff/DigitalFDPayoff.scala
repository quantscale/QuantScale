package quantscale.fdm.payoff
import quantscale.fdm.Epsilon

class DigitalFDPayoff (
    isCall: Boolean,
    strike: Double, 
    payout : Double,
    expiryTime: Double) extends FDPayoff {
  
      override def copy(): FDPayoff = {
        return new DigitalFDPayoff(isCall, strike, payout, expiryTime)
      }

  val putCallSign: Int = if (isCall) 1 else -1;

  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous

  override def eval(): Unit = {
    if (Math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      var i: Int = 0;
      while (i < state.size) {
        state.price(i) = if (putCallSign * (space(i) - strike) > 0) payout else 0
        i += 1;
      }
      _isDiscontinuous = true
    } else {
      _isDiscontinuous = false
    }
  }


}
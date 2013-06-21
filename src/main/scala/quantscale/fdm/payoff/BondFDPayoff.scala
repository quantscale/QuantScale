package quantscale.fdm.payoff

import quantscale.fdm.Epsilon

class BondFDPayoff(val expiryTime: Double, val coupon: Double=1.0) extends FDPayoff {

  override def copy(): FDPayoff = {
    return new BondFDPayoff(expiryTime, coupon)
  }
    
    override def eval(): Unit = {
    if (math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      var i: Int = 0;
      while (i < state.size) {
        state.price(i) = coupon
        i += 1;
      }
    }
  }
}
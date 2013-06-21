package quantscale.fdm.payoff
import quantscale.fdm.Epsilon

class ButterflyFDPayoff(
  isCall: Boolean,
  K1: Double, 
  K2: Double,
  expiryTime: Double) extends FDPayoff {
  
      override def copy(): FDPayoff = {
        return new ButterflyFDPayoff(isCall, K1, K2, expiryTime)
      }
        
  val putCallSign: Int = if (isCall) 1 else -1
  val Kmid = (K1 + K2) * 0.5

  private var _isDiscontinuous = false
  
  override def isDiscontinuous = _isDiscontinuous
  
  override def eval(): Unit = {
    if (math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      var i = 0
      while (i < state.size) {
        state.price(i) = math.max(putCallSign * (space(i) - K1), 0) -
          2 * math.max(putCallSign * (space(i) - Kmid), 0) +
          math.max(putCallSign * (space(i) - K2), 0)
        i += 1
      }
      _isDiscontinuous = true
    } else {
      _isDiscontinuous = false
    }
  }

}
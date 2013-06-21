package quantscale.fdm.payoff
import quantscale.fdm.Epsilon

class KOBarrierFDPayoff(
  val isCall: Boolean,
  val strike: Double,
  val upLevel: Double = Double.MaxValue,
  val downLevel: Double = Double.MinValue,
  val expiryTime: Double,
  val koTimes: Array[Double],
  val isContinuous: Boolean = false) extends FDPayoff {

      override def copy(): FDPayoff = {
        return new KOBarrierFDPayoff(isCall, strike, upLevel, downLevel, expiryTime, koTimes, isContinuous)
      }
        
  val putCallSign: Int = if (isCall) 1 else -1

  private var koTimeIndex = 0
  
  private var _isDiscontinuous = false
  
  override def isDiscontinuous = _isDiscontinuous
  
  override def initState(underlyingLevel : Array[Double]) {
        super.initState(underlyingLevel)
        koTimeIndex = koTimes.length - 1
    }
  
  override def eval(): Unit = {
      _isDiscontinuous = false
    if (math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      var i: Int = 0
      while (i < state.size) {
        state.price(i) = math.max(putCallSign * (space(i) - strike), 0)
        i += 1
      }
      _isDiscontinuous = true
    }

    while (koTimeIndex >= 0 && (koTimes(koTimeIndex) > time + Epsilon.MACHINE_EPSILON_SQRT)) {
      koTimeIndex -= 1
    }
    if (isContinuous || (koTimeIndex >= 0 && math.abs(koTimes(koTimeIndex) - time) < Epsilon.MACHINE_EPSILON_SQRT)) {
      //note : mesh should be ideally truncated if barrier is fully continuous
      var i: Int = 0
      _isDiscontinuous = true
      while (i < state.size) {
        if (space(i) > upLevel) {
          state.price(i) = 0
        } else if (space(i) < downLevel) {
          state.price(i) = 0
        }
        i += 1
      }
    }
  }
}
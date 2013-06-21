package quantscale.fdm.payoff
import quantscale.fdm.Epsilon
import quantscale.math.CubicSpline
import quantscale.fdm.State
import quantscale.math.Extrapolation1D
import quantscale.math.LinearExtrapolation

class AccumulatorKODAFDPayoff(
  val buySell: Int,
  val periodAccrual: BackwardSchedule,
  val knockouts: BackwardSchedule,
  val payments: BackwardSchedule,
  val forwardPrice: Double,
  val accrualAbove: Double,
  val accrualBelow: Double,
  val knockOutBarrier: Double,
  val isKOAfterAccrual: Boolean = false) extends FDPayoff {

  override def copy(): FDPayoff = {
    return new AccumulatorKODAFDPayoff(buySell, periodAccrual.copy(), knockouts.copy(), payments.copy(), forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKOAfterAccrual)
  }
  
  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous
  private var tp = 0.0

  override def stateDimensions(): Int = { return 2 }

  override def setTime(t: Double) {
    super.setTime(t)
    periodAccrual.advance(t)
    knockouts.advance(t)
    payments.advance(t)
  }

  override def eval(): Unit = {
    _isDiscontinuous = false
    if (payments.isMeshTime) {
      val pv = state.values(1)
      for (j <- 0 until pv.length) {
        pv(j) = 1
      }
      tp = time
      //new coupon +  curr value
      //      val discountPV = state.values(AccumulatorFDPayoff.discretizationSize)
      //      for (j <- 0 until discountPV.length) {
      //          discountPV(j) = 1.0
      //        }
    }

    if (isKOAfterAccrual && knockouts.isMeshTime) {
      val pv = state.values(0)
      for (j <- 0 until pv.length) {
        if (space(j) >= knockOutBarrier) {
          pv(j) = 0
        }
      }
      _isDiscontinuous = true
    }

    if (periodAccrual.isMeshTime) {
      val pv = state.values(0)
      val discountPV = state.values(1)
      //      val df = math.exp((time-tp)*0.03)
      for (j <- 0 until pv.length) {
        val quote = space(j)
        val amount = if (quote > forwardPrice) accrualAbove else accrualBelow
        val df = discountPV(j)
        pv(j) += buySell * amount * (space(j) - forwardPrice * df)
      }
    }

    if (!isKOAfterAccrual && knockouts.isMeshTime) {
      val pv = state.values(0)
      for (j <- 0 until pv.length) {
        if (space(j) >= knockOutBarrier) {
          pv(j) = 0 //buySell * shares * (space(j) - forwardPrice);
        }
      }
      _isDiscontinuous = true
    }
  }

  override def initState(underlyingLevel: Array[Double]) {
    super.initState(underlyingLevel)
    periodAccrual.reset()
    knockouts.reset()
    payments.reset()
  }
}

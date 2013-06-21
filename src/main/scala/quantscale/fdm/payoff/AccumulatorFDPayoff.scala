package quantscale.fdm.payoff
import quantscale.fdm.Epsilon
import quantscale.math.CubicSpline
import quantscale.fdm.State
import quantscale.math.Extrapolation1D
import quantscale.math.LinearExtrapolation

object AccumulatorFDPayoff {
  val discretizationSize = 3
}

class AccumulatorFDPayoff(
  val buySell: Int,
  val periodAccrual: BackwardSchedule,
  val knockouts: BackwardSchedule,
  val payments: BackwardSchedule,
  val forwardPrice: Double,
  val accrualAbove: Double,
  val accrualBelow: Double,
  val knockOutBarrier: Double,
  val isKOAfterAccrual: Boolean = true) extends FDPayoff {

  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous

  private var totalPeriodAccruingShares: Array[Double] = null
  private var callValueAtAccruingShares: Array[Double] = null
  private var stateClone: State = null
  private var extrapolation: Extrapolation1D = new LinearExtrapolation()

  override def copy() : FDPayoff = {
    return new AccumulatorFDPayoff(buySell,  periodAccrual.copy(), knockouts.copy(), payments.copy(), forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKOAfterAccrual)
  }
  override def stateDimensions() = AccumulatorFDPayoff.discretizationSize
  
  override def setTime(t: Double) {
    super.setTime(t)
    periodAccrual.advance(t)
    knockouts.advance(t)
    payments.advance(t)
  }

  override def eval(): Unit = {
    if (payments.isMeshTime) {
      var pv0 = state.values(0)
      var i = AccumulatorFDPayoff.discretizationSize-1
      while (i>=0) {
        val shares = totalPeriodAccruingShares(i)
        val pv = state.values(i)
        for (j <- 0 until pv.length) {
          pv(j) = buySell * shares * (space(j) - forwardPrice) + pv0(j)
        }
        //new coupon +  curr value
        i-=1
      }
    }

    if (isKOAfterAccrual && knockouts.isMeshTime) {
      for (i <- 0 until AccumulatorFDPayoff.discretizationSize) {
        var shares = totalPeriodAccruingShares(i)
        val pv = state.values(i)
        for (j <- 0 until pv.length) {
          if (space(j) >= knockOutBarrier) {
            pv(j) = buySell * shares * (space(j) - forwardPrice);
          }
        }
      }
      _isDiscontinuous = true
    } else {
      _isDiscontinuous = false
    }
    
    if (periodAccrual.isMeshTime) {
      for (i <- 0 until AccumulatorFDPayoff.discretizationSize) {
        //call ti-(S,I) = call ti+ (S,I+S)
        val pv = stateClone.values(i)
        for (j <- 0 until pv.length) {
          val quote = space(j)
          val amount = if (quote > forwardPrice) accrualAbove else accrualBelow
          pv(j) = interpolateValue(state.values, j, totalPeriodAccruingShares(i) + amount)
        }
      }
      for (i <- 0 until AccumulatorFDPayoff.discretizationSize) {
        System.arraycopy(stateClone.values(i), 0, state.values(i), 0, state.values(i).length);
      }
    }

     if (!isKOAfterAccrual && knockouts.isMeshTime) {
      for (i <- 0 until AccumulatorFDPayoff.discretizationSize) {
        var shares = totalPeriodAccruingShares(i)
        val pv = state.values(i)
        for (j <- 0 until pv.length) {
          if (space(j) >= knockOutBarrier) {
            pv(j) = buySell * shares * (space(j) - forwardPrice);
          }
        }
      }
      _isDiscontinuous = true
    } else {
      _isDiscontinuous = false
    }
  }

  def interpolateValue(pv: Array[Array[Double]], underlyingIndex: Int, average: Double) : Double = {
    //call[i] corresponds to average[i].

    for (i <- 0 until AccumulatorFDPayoff.discretizationSize) {
      callValueAtAccruingShares(i) = pv(i)(underlyingIndex)
    }
    //        interp.setY(callValueAtAverage);
    extrapolation.init(totalPeriodAccruingShares, callValueAtAccruingShares)
    if (extrapolation.isInExtrapolationRange(average)) {
      return extrapolation.value(average)
    } else {
      val interp = CubicSpline.makeBesselSpline(totalPeriodAccruingShares, callValueAtAccruingShares)
      return interp.value(average)
    }
  }

  override def initState(underlyingLevel: Array[Double]) {
    super.initState(underlyingLevel)
    periodAccrual.reset()
    knockouts.reset()
    payments.reset()
    totalPeriodAccruingShares = new Array[Double](AccumulatorFDPayoff.discretizationSize)
    val periodAccrualSize = periodAccrual.length / payments.length;

    for (i <- 0 until AccumulatorFDPayoff.discretizationSize) {
      totalPeriodAccruingShares(i) = periodAccrualSize * ((accrualAbove) * i / (AccumulatorFDPayoff.discretizationSize - 1));
    }
    stateClone = state.copy()
    callValueAtAccruingShares = new Array[Double](AccumulatorFDPayoff.discretizationSize)
  }

}
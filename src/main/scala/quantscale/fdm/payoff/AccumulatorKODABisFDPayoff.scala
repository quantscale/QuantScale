package quantscale.fdm.payoff

import quantscale.fdm.State
import quantscale.math.CubicSpline
import quantscale.math.Extrapolation1D
import quantscale.math.LinearExtrapolation
import java.util.Arrays
import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting
import quantscale.fdm.Epsilon

class AccumulatorKODABisFDPayoff(
  val buySell: Int,
  val periodAccrual: BackwardSchedule,
  val knockouts: BackwardSchedule,
  val payments: BackwardSchedule,
  val forwardPrice: Double,
  val accrualAbove: Double,
  val accrualBelow: Double,
  val knockOutBarrier: Double,
  val isKOAfterAccrual: Boolean = false,
  val tarnLevel: Double = Double.MaxValue) extends FDPayoff() {

  override def copy(): FDPayoff = {
    return new AccumulatorKODABisFDPayoff(buySell,  periodAccrual.copy(), knockouts.copy(), payments.copy(), forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKOAfterAccrual, tarnLevel)
  }
    
  private var _isDiscontinuous = false

  override def isDiscontinuous(): Boolean = {
    return _isDiscontinuous
  }

  private var totalPeriodAccruingShares: Array[Double] = null
  private var totalShares: Array[Double] = null
  private var callValueAtAccruingShares: Array[Double] = null
  private var stateClone: State = null
  private var extrapolation: Extrapolation1D = new LinearExtrapolation()

  private var _discretizationSize = 0

  initTotalSharesArray()

  def initTotalSharesArray() {
    var totalSharesList = new ArrayBuffer[Double]()
    val maxAccrual = math.max(accrualAbove, accrualBelow)
    val tarnLevelMax = tarnLevel + maxAccrual

    val na = tarnLevelMax / accrualAbove
    val nb = tarnLevelMax / accrualBelow
    if (na > 200 || nb > 200) {
      val periodAccrualSize = periodAccrual.length // / payments.length;
      if (tarnLevelMax > periodAccrualSize * maxAccrual) {
        //tarn not effective
        _discretizationSize = 2
        totalPeriodAccruingShares = new Array[Double](_discretizationSize)
        totalPeriodAccruingShares(0) = 0
        totalPeriodAccruingShares(1) = math.min(accrualAbove, accrualBelow)
        //        totalPeriodAccruingShares(2) = maxAccrual
      } else {
        _discretizationSize = 200
        totalPeriodAccruingShares = new Array[Double](_discretizationSize)
        for (i <- 0 until _discretizationSize) {
          totalPeriodAccruingShares(i) = tarnLevelMax * i / (_discretizationSize - 1);
        }
      }
    } else {
      var i = 0
      var j = 0
      while (i < na) {
        j = 0
        while (j < nb) {
          val shares = i * accrualAbove + j * accrualBelow
          if (shares <= tarnLevelMax) {
            totalSharesList += shares
          }
          j += 1
        }
        i += 1
      }
      totalPeriodAccruingShares = totalSharesList.toArray
      Sorting.quickSort(totalPeriodAccruingShares)
      totalSharesList = new ArrayBuffer[Double]()
      var previous = totalPeriodAccruingShares(0)
      totalSharesList += previous
      i = 1
      while (i < totalPeriodAccruingShares.length) {
        if (totalPeriodAccruingShares(i) - previous > Epsilon.MACHINE_EPSILON_SQRT) {
          previous = totalPeriodAccruingShares(i)
          totalSharesList += previous
        }
        i += 1
      }
      totalPeriodAccruingShares = totalSharesList.toArray
      _discretizationSize = totalPeriodAccruingShares.length
      println("size=" + _discretizationSize + " " + totalSharesList)
    }
  }

  override def stateDimensions(): Int = {
    return _discretizationSize
  }

  override def setTime(t: Double) {
    super.setTime(t)
    periodAccrual.advance(t)
    knockouts.advance(t)
    payments.advance(t)
  }

  override def eval(): Unit = {
    _isDiscontinuous = false
    if (payments.isMeshTime) {
      var pv0 = state.values(0)
      var i = _discretizationSize - 1
      while (i >= 0) {
        val shares = totalPeriodAccruingShares(i)
        
        val pv = state.values(i)
          for (j <- 0 until pv.length) {
            pv(j) = buySell * shares * (space(j) - forwardPrice) + pv0(j)
          }
        
        //new coupon +  curr value
        i -= 1
      }
    }

    if (isKOAfterAccrual && knockouts.isMeshTime) {
      for (i <- 0 until _discretizationSize) {
        var shares = totalPeriodAccruingShares(i)
        val pv = state.values(i)
        for (j <- 0 until pv.length) {
          if (space(j) >= knockOutBarrier) {
            pv(j) = buySell * shares * (space(j) - forwardPrice);
          }
        }
      }
      _isDiscontinuous = true
    }

    if (periodAccrual.isMeshTime) {
      for (i <- 0 until _discretizationSize) {
        //call ti-(S,I) = call ti+ (S,I+S)
        val pv = stateClone.values(i)
        for (j <- 0 until pv.length) {
          val quote = space(j)
          val amount = if (quote > forwardPrice) accrualAbove else accrualBelow
          val total = totalPeriodAccruingShares(i) + amount
          totalShares(i) += amount
          if (totalShares(i) >= tarnLevel) {
            pv(j) = buySell * total * (space(j) - forwardPrice)
          } else {
            pv(j) = interpolateValue(state.values, j, total)
          }

        }
      }
      for (i <- 0 until _discretizationSize) {
        System.arraycopy(stateClone.values(i), 0, state.values(i), 0, state.values(i).length);
      }
    }

    if (!isKOAfterAccrual && knockouts.isMeshTime) {
      for (i <- 0 until _discretizationSize) {
        var shares = totalPeriodAccruingShares(i)
        val pv = state.values(i)
        for (j <- 0 until pv.length) {
          if (space(j) >= knockOutBarrier) {
            pv(j) = buySell * shares * (space(j) - forwardPrice);
          }
        }
      }
      _isDiscontinuous = true
    }
  }

  def interpolateValue(pv: Array[Array[Double]], underlyingIndex: Int, average: Double): Double = {
    var i = 0
    if (average < totalPeriodAccruingShares(0)) {
      i = 0
    } else if (average >= totalPeriodAccruingShares(totalPeriodAccruingShares.length - 1)) {
      i = totalPeriodAccruingShares.length - 2
    } else {
      i = Arrays.binarySearch(totalPeriodAccruingShares, average);
      i = if (i < 0) -i - 2 else i;
    }
    val value = (pv(i + 1)(underlyingIndex) - pv(i)(underlyingIndex)) / (totalPeriodAccruingShares(i + 1) - totalPeriodAccruingShares(i)) * (average - totalPeriodAccruingShares(i)) + pv(i)(underlyingIndex)
    return value
    //    //call[i] corresponds to average[i].
    //
    //    for (i <- 0 until AccumulatorKODABisFDPayoff.discretizationSize) {
    //      callValueAtAccruingShares(i) = pv(i)(underlyingIndex)
    //    }
    //    //        interp.setY(callValueAtAverage);
    //    extrapolation.init(totalPeriodAccruingShares, callValueAtAccruingShares)
    //    if (extrapolation.isInExtrapolationRange(average)) {
    //      return extrapolation.value(average)
    //    } else {
    //      val interp = CubicSpline.makeBesselSpline(totalPeriodAccruingShares, callValueAtAccruingShares)
    //      return interp.value(average)
    //    }
  }

  override def initState(underlyingLevel: Array[Double]) {
    super.initState(underlyingLevel)
    periodAccrual.reset()
    knockouts.reset()
    payments.reset()

    stateClone = state.copy()
    callValueAtAccruingShares = new Array[Double](_discretizationSize)
    totalShares = new Array[Double](_discretizationSize)
    for (i <- 0 until _discretizationSize) {
      totalShares(i) = totalPeriodAccruingShares(i)
    }
  }

}
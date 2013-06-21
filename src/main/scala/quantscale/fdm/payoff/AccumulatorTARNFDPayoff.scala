package quantscale.fdm.payoff
import quantscale.fdm.Epsilon

import quantscale.math.CubicSpline
import quantscale.fdm.State
import quantscale.math.Extrapolation1D
import quantscale.math.LinearExtrapolation
import scala.util.Sorting
import scala.collection.mutable.ArrayBuffer
import java.util.Arrays
import quantscale.math.CubicPP

object AccumulatorTARNInterpolation extends Enumeration {
  type AccumulatorTARNInterpolation = Value
  val LINEAR, CUBIC = Value
}

import AccumulatorTARNInterpolation._

class AccumulatorTARNFDPayoff(
  val buySell: Int,
  val periodAccrual: BackwardSchedule,
  val knockouts: BackwardSchedule,
  val payments: BackwardSchedule,
  val forwardPrice: Double,
  val accrualAbove: Double,
  val accrualBelow: Double,
  val knockOutBarrier: Double,
  val isKOAfterAccrual: Boolean = false,
  val tarnLevel: Double = Double.MaxValue) extends FDPayoff {

  override def copy(): FDPayoff = {
    return new AccumulatorTARNFDPayoff(buySell, periodAccrual.copy(), knockouts.copy(), payments.copy(), forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKOAfterAccrual, tarnLevel)
  }

  private var _isDiscontinuous = false

  private var _discretizationSize = 1
  private var totalQuantity: Array[Double] = null
  private var stateClone: State = null
  private var tp = 0.0
  private var y, yPrime: Array[Double] = null
  private var quantityAboveIndex, quantityBelowIndex: Array[Int] = null

  private var interp: CubicPP = null
  private var interpolatorType = CUBIC
  private var currentAccrualIndex = 0
  initTotalSharesArray()

  override def isDiscontinuous = _isDiscontinuous
  def initTotalSharesArray() {
    var totalSharesList = new ArrayBuffer[Double]()
    val maxAccrual = math.max(accrualAbove, accrualBelow)
    val tarnLevelMax = tarnLevel + maxAccrual

    val na = tarnLevelMax / accrualAbove
    val nb = tarnLevelMax / accrualBelow
    if (na > 200 || nb > 200) {
      val periodAccrualSize = periodAccrual.length;
      if (tarnLevelMax > periodAccrualSize * maxAccrual) {
        //tarn not effective
        _discretizationSize = 1
        totalQuantity = new Array[Double](_discretizationSize)
        totalQuantity(0) = 0
      } else {
        _discretizationSize = 200
        totalQuantity = new Array[Double](_discretizationSize)
        for (i <- 0 until _discretizationSize) {
          totalQuantity(i) = tarnLevelMax * i / (_discretizationSize - 1);
        }
      }
    } else {
      var i = 0
      //      while (i <= na) {
      //        totalSharesList += i * accrualAbove
      //        i += 1
      //      }
      //      
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
      totalQuantity = totalSharesList.toArray
      Sorting.quickSort(totalQuantity)
      totalSharesList = new ArrayBuffer[Double]()
      var previous = totalQuantity(0)
      totalSharesList += previous
      i = 1
      while (i < totalQuantity.length) {
        if (totalQuantity(i) - previous > Epsilon.MACHINE_EPSILON_SQRT) {
          previous = totalQuantity(i)
          totalSharesList += previous
        }
        i += 1
      }
      totalQuantity = totalSharesList.toArray
      _discretizationSize = totalQuantity.length

    }
    println("size=" + _discretizationSize + " " + totalSharesList)
  }
  override def stateDimensions(): Int = { return _discretizationSize + 1 }

  override def setTime(t: Double) {
    super.setTime(t)
    periodAccrual.advance(t)
    knockouts.advance(t)
    payments.advance(t)
  }

  override def eval(): Unit = {
    _isDiscontinuous = false
    if (payments.isMeshTime) {
      val pv = state.values(_discretizationSize)
      for (j <- 0 until pv.length) {
        pv(j) = 1
      }
    }

    if (isKOAfterAccrual && knockouts.isMeshTime) {

      val discountPV = state.values(_discretizationSize)

      for (i <- 0 until _discretizationSize) {
        val pv = state.values(i)
        for (j <- 0 until pv.length) {
          val quote = space(j)
          if (space(j) >= knockOutBarrier) {
            pv(j) = 0 //  buySell*accrualAbove * (forwardPrice*discountPV(j)-forwardPrice)
          }
        }
      }
      _isDiscontinuous = true
    }

    if (periodAccrual.isMeshTime) {

      //      for (i <- 0 until _discretizationSize) {
      //        for (j <- 0 until state.size) {
      //          val quote = space(j)
      //          val amount = if (quote > forwardPrice) accrualAbove else accrualBelow
      //          initInterpolator(totalQuantity, state.values, j)
      //          val pv = stateClone.values(i)
      //          pv(j) = interpolateValue(state.values, j, totalQuantity(i) + amount)
      //        }
      //      }

      //for a given time, graph of state.values(i,j) to see where it is sensitive
      for (j <- 0 until state.size) {
        val quote = space(j)
        val amount = if (quote > forwardPrice) accrualAbove else accrualBelow
        initInterpolator(totalQuantity, state.values, j)
        for (i <- 0 until _discretizationSize) {
          val pv = stateClone.values(i)
          pv(j) = interpolateValue(state.values, j, totalQuantity(i) + amount)
        }
      }
      for (i <- 0 until _discretizationSize) {
        System.arraycopy(stateClone.values(i), 0, state.values(i), 0, state.values(i).length);
      }
      val discountPV = state.values(_discretizationSize)

      for (i <- 0 until _discretizationSize) {
        val pv = state.values(i)
        if (totalQuantity(i) >= tarnLevel) {
          for (j <- 0 until pv.length) {
            pv(j) = 0 //buySell * shares * (space(j) - forwardPrice);
          }
        } else {
          for (j <- 0 until pv.length) {
            val quote = space(j)
            val amount = if (quote > forwardPrice) accrualAbove else accrualBelow
            val df = discountPV(j)
            pv(j) += buySell * amount * (space(j) - forwardPrice * df)
          }
        }
      }
      currentAccrualIndex += 1

    }

    if (!isKOAfterAccrual && knockouts.isMeshTime) {
      val discountPV = state.values(_discretizationSize)
      val pv = state.values(0)
      for (j <- 0 until pv.length) {
        if (space(j) >= knockOutBarrier) {
          pv(j) = 0

        }
      }
      _isDiscontinuous = true
    }
  }

  def initInterpolator(x: Array[Double], pv: Array[Array[Double]], underlyingIndex: Int) {
    if (interpolatorType != LINEAR && totalQuantity.length > 1) {
      var i = 0
      while (i < y.length) {
        y(i) = pv(i)(underlyingIndex)
        i += 1
      }
      interp = CubicSpline.makeBesselSpline(x, y)
    }
  }

  def interpolateValue(pv: Array[Array[Double]], underlyingIndex: Int, average: Double): Double = {
    interpolatorType match {
      case LINEAR =>
        return linearInterpolateValue(pv, underlyingIndex, average)
      case CUBIC =>
        return cubicInterpolateValue(pv, underlyingIndex, average)
    }
  }

  private def cubicInterpolateValue(pv: Array[Array[Double]], underlyingIndex: Int, average: Double): Double = {
    if (totalQuantity.length == 1) {
      return pv(0)(underlyingIndex)
    }
    return interp.value(average)
  }

  private def linearInterpolateValue(pv: Array[Array[Double]], underlyingIndex: Int, average: Double): Double = {
    var i = 0
    if (average < totalQuantity(0)) {
      i = 0
    } else if (average >= totalQuantity(totalQuantity.length - 1)) {
      i = totalQuantity.length - 2
    } else {
      i = Arrays.binarySearch(totalQuantity, average);
      i = if (i < 0) -i - 2 else i;
    }
    val value = (pv(i + 1)(underlyingIndex) - pv(i)(underlyingIndex)) / (totalQuantity(i + 1) - totalQuantity(i)) * (average - totalQuantity(i)) + pv(i)(underlyingIndex)
    return value
  }

  override def initState(underlyingLevel: Array[Double]) {
    super.initState(underlyingLevel)
    periodAccrual.reset()
    knockouts.reset()
    payments.reset()
    stateClone = state.copy()
    y = new Array[Double](totalQuantity.length)
    yPrime = new Array[Double](totalQuantity.length)
    currentAccrualIndex = 0
  }
}

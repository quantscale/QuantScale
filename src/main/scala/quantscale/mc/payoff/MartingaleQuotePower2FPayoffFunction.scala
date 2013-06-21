package quantscale.mc.payoff

import quantscale.mc.BSM2FMartingaleQuotePower

class MartingaleQuotePower2FPayoffFunction(val martingalePower: BSM2FMartingaleQuotePower) {
  private var _timeIndex = 0
  private var _times : Array[Double]= null
  def reset() {
    _timeIndex = 0
    _times = martingalePower.martingalePowerPath(0).time
  }
  def setTime(time: Double) {
    while (_timeIndex < _times.length-1 && time > _times(_timeIndex)) {
      _timeIndex +=1
    }
  }

  def eval(power: Int, asset: Int): Double = {
    if (power == 0) return 1.0
    val logS = martingalePower.martingalePowerPath(power-1).values(_timeIndex)(asset)
    return math.exp(logS)
  }
  
  def evalCross(): Double = {
    val logS = martingalePower.martingalePowerPathCross(0).values(_timeIndex)(0)
    return math.exp(logS)
  }
}
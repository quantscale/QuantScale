package quantscale.mc.payoff

import quantscale.mc.BSM1FMartingaleQuotePower

class MartingaleQuotePowerPayoffFunction(val martingalePower: BSM1FMartingaleQuotePower) {
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
}
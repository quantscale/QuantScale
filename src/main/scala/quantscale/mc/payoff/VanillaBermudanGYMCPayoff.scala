package quantscale.mc.payoff

import quantscale.mc.MCPath
import quantscale.mc.BSM1FMartingaleQuotePower

class VanillaBermudanGYMCPayoff(isCall: Boolean, strike: Double, _exerciseTimes: Array[Double], _paymentTimes: Array[Double]) extends MCPayoff {

  private var _currentTime: Double = 0.0
  private var _exerciseCashflow = new Cashflow(0, 0, 0)
  private var _cashflows = Array[Cashflow](_exerciseCashflow)
  private var sign = if (isCall) 1 else -1
  private var _indexNextExercise = 0
  private var _isExercise = false
  private var _x = 0.0

  def regressionVariablesDimension = 5
  var martingalePower: MartingaleQuotePowerPayoffFunction = new MartingaleQuotePowerPayoffFunction(new BSM1FMartingaleQuotePower(4, 0))

  def reset() {
    _exerciseCashflow.amount = 0.0
    _indexNextExercise = 0
    _isExercise = false
    martingalePower.reset()
  }

  def setTime(time: Double) {
    _currentTime = time

    while (_indexNextExercise < _exerciseTimes.length - 1 && time > _exerciseTimes(_indexNextExercise)) {
      _indexNextExercise += 1
    }

    if (time == _exerciseTimes(_indexNextExercise)) {
      _isExercise = true
    }
    martingalePower.setTime(time)
  }

  def exerciseTimes = evaluationTimes

  def payments = _cashflows


  def evalRegressionVariables(quotes: Array[Double], independentVariables: Array[Double]) {
    val quote = quotes(0)
    independentVariables(0) = 1
    independentVariables(1) = martingalePower.eval(1, 0)
    independentVariables(2) = martingalePower.eval(2, 0)
    independentVariables(3) = martingalePower.eval(3, 0)
    independentVariables(4) = martingalePower.eval(4, 0)
  }

  def eval(quotes: Array[Double]) {
    if (_isExercise) {
      _x = math.exp(quotes(0))
      //FIXME move out exp of here
      _exerciseCashflow.amount = math.max(0.0, sign * (_x - strike))
      _exerciseCashflow.fixingTime = _currentTime
      _exerciseCashflow.paymentTime = _paymentTimes(_indexNextExercise)
      _exerciseCashflow.isExercise = true
    } else {
      _exerciseCashflow.amount = 0
      _exerciseCashflow.isExercise = false
    }
  }

  def paymentTimes = _paymentTimes

  def evaluationTimes = _exerciseTimes
}
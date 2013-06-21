package quantscale.mc.payoff

import quantscale.mc.MCPath

//Basket of 2
class VanillaBasketBermudanMCPayoff(isCall: Boolean, strike: Double, weights: Array[Double], fixings: Array[Double], _exerciseTimes: Array[Double], _paymentTimes: Array[Double]) extends MCPayoff {

  private var _currentTime : Double = 0.0
  private var _exerciseCashflow = new Cashflow(0, 0, 0)
  private var _cashflows = Array[Cashflow](_exerciseCashflow)
  private var sign = if (isCall) 1 else -1
  private var _indexNextExercise = 0
  private var _isExercise = false
  private var _x = 0.0
  private var _y = 0.0
  
  def reset() {
    _exerciseCashflow.amount = 0.0
    _indexNextExercise = 0
    _isExercise = false
  }
  
  def setTime(time : Double) {
    _currentTime = time
    
    while (_indexNextExercise < _exerciseTimes.length-1 && time > _exerciseTimes(_indexNextExercise)) {
      _indexNextExercise +=1
    }
    
    if (time == _exerciseTimes(_indexNextExercise)) {
      _isExercise = true
    }
  }
  
  def exerciseTimes = evaluationTimes
  
  def payments = _cashflows
  
  def regressionVariablesDimension = 6
  
  def evalRegressionVariables(quotes: Array[Double], independentVariables: Array[Double]) {
    val quote = quotes(0)
    independentVariables(0) = 1
    independentVariables(1) = _x
    independentVariables(2) = _x*_x
    independentVariables(3) = _y
    independentVariables(4) = _y*_y
    independentVariables(5) = _x*_y
  }
  
  def eval(quotes : Array[Double]) {
    if (_isExercise) {
       _x = math.exp(quotes(0))
       _y = math.exp(quotes(1))
      //FIXME move out exp of here
       var sum = 0.0
       var i = 0
       while (i < weights.length) {
         sum += math.exp(quotes(i))/fixings(i)*weights(i)
         i+=1
       }
      _exerciseCashflow.amount = math.max(0.0,sign*( sum-strike)) 
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
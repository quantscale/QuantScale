package quantscale.mc.payoff

import quantscale.mc.MCPath

class VanillaMCPayoff(isCall: Boolean, strike: Double, maturityTime: Double, paymentTime: Double) extends MCPayoff {

  private var _currentTime : Double = 0.0
  private var _exerciseCashflow = new Cashflow(0, maturityTime, paymentTime)
  private var _cashflows = Array[Cashflow](_exerciseCashflow)
  private var sign = if (isCall) 1 else -1
  
  def reset() {
    _exerciseCashflow.amount = 0.0
  }
  
  def setTime(time : Double) {
    _currentTime = time
  }
  
  def exerciseTimes = evaluationTimes
  
  def payments = _cashflows
  
  def regressionVariablesDimension = 3
  
  def evalRegressionVariables(quotes: Array[Double], independentVariables: Array[Double]) {
    val quote = quotes(0)
    independentVariables(0) = 1
    independentVariables(1) = quote
    independentVariables(2) = quote*quote
  }
  
  def eval(quotes : Array[Double]) {
    if (_currentTime == maturityTime) {
      //FIXME move out exp of here
      _exerciseCashflow.amount = math.max(0.0,sign*( math.exp(quotes(0))-strike)) 
      _exerciseCashflow.isExercise = true
    } else {
      _exerciseCashflow.amount = 0
      _exerciseCashflow.isExercise = false
    }
  }

  def paymentTimes: Array[Double] = {
    Array(paymentTime)
  }

  def evaluationTimes: Array[Double] = {
    Array(maturityTime)
  }
}
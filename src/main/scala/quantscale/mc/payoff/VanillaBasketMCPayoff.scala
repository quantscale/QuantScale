package quantscale.mc.payoff

import quantscale.mc.payoff.MCPayoff
import scala.Array

class VanillaBasketMCPayoff(isCall: Boolean, strike: Double, weights: Array[Double], maturityTime: Double, paymentTime: Double) extends AdjointMCPayoff {

  private var _currentTime : Double = 0.0
  private val _exerciseCashflow = new Cashflow(0, maturityTime, paymentTime)
  private val _cashflows = Array[Cashflow](_exerciseCashflow)
  private val sign = if (isCall) 1 else -1

  def reset() {
    _exerciseCashflow.amount = 0.0
  }

  def setTime(time : Double) {
    _currentTime = time
  }

  def exerciseTimes = evaluationTimes

  def payments = _cashflows

  def regressionVariablesDimension = weights.length+1

  def evalRegressionVariables(quotes: Array[Double], independentVariables: Array[Double]) {
    independentVariables(0) = 1
    var i = 1
    while (i <= weights.length) {
      independentVariables(i) = quotes(i-1)
      i+=1
    }
  }

  def evalAdjoint(quotes: Array[Double], quotes_b: Array[Double]) {
    if (_currentTime == maturityTime) {
      //FIXME move out exp of here
      var i = 0
      var sum = 0.0
      while (i < quotes.length) {
        val xi = math.exp(quotes(i))
        quotes_b(i) = xi
        sum += weights(i)*xi
        i+=1
      }
      val x = sign*(sum-strike)
      val payoff = math.max(0.0,x)
      var payoff_b = 1.0
      var x_b = 0.0
      if (x > 0) {
        x_b = payoff_b
      }
      var sum_b = sign*x_b
      i = quotes.length-1
      while (i>=0) {
        val xi_b = weights(i)*sum_b
        quotes_b(i) = xi_b*quotes_b(i) //xi_b*math.exp(quotes(i))
        i-=1
      }
      _exerciseCashflow.amount = payoff
      _exerciseCashflow.isExercise = true
    } else {
      _exerciseCashflow.amount = 0
      _exerciseCashflow.isExercise = false
    }
  }

  def eval(quotes : Array[Double]) {
    if (_currentTime == maturityTime) {
      //FIXME move out exp of here
      var i = 0
      var sum = 0.0
      while (i < quotes.length) {
        val xi = math.exp(quotes(i))
        sum += weights(i)*xi
        i+=1
      }
      val x = sign*(sum-strike)
      val payoff = math.max(0.0,x)
      _exerciseCashflow.amount = payoff
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

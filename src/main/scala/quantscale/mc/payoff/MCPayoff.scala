package quantscale.mc.payoff

import quantscale.mc.MCPath

class Cashflow(var amount: Double, var fixingTime: Double, var paymentTime: Double, var isExercise: Boolean = false) {

}

trait AdjointMCPayoff extends MCPayoff {
  def evalAdjoint(quotes: Array[Double], quotes_b:Array[Double])

}
trait MCPayoff {
  //implementation should reuse the payments to avoid reconstructing the object in the monte carlo loop
  def eval(quotes: Array[Double])
  def payments: Array[Cashflow] //payments at current time
  def paymentTimes: Array[Double]
  def exerciseTimes: Array[Double]
  def evalRegressionVariables(quotes: Array[Double], independentVariables: Array[Double])
  def regressionVariablesDimension: Int
  def evaluationTimes: Array[Double]
  def setTime(time: Double)
  def reset()
}
package quantscale.mc

import quantscale.mc.payoff.MCPayoff

class LongstaffSchwartzPathPricer(payoff: MCPayoff, discounter: Discounter) {

  private var _evaluationTimes = payoff.evaluationTimes
  private var _exerciseTimes = payoff.exerciseTimes
  private var _independentVariablesDimension = payoff.regressionVariablesDimension //need to be fixed dimension

  private var _dfCache: Map[Double, Double] = Map()

  initDiscountFactors(payoff.paymentTimes)

  def initDiscountFactors(paymentTimes: Array[Double]) {
    for (time <- paymentTimes) {
      val df = discounter.df(time)
      _dfCache += (time -> df)
    }
  }

  def evaluationTimes = _evaluationTimes
  def exerciseTimes = _exerciseTimes
  def independentVariablesDimension = _independentVariablesDimension

  def getDiscountFactor(time: Double): Double = {
    val dfO = _dfCache.get(time)
    if (dfO.isEmpty) {
      val df = discounter.df(time)
      _dfCache += (time -> df)
      return df
    } else {
      return dfO.get
    }
  }

  
  def eval(path: MCPath, exerciseFlows: Array[Double], independentVariables: Array[Array[Double]]) : Double = {
    var sum = 0d
    var ie = 0
    payoff.reset()
    var indexExercise = 0
    var exerciseTime = _exerciseTimes(indexExercise)
    while (ie < evaluationTimes.length) {
      val quotes = path.values(ie)
      val evaluationTime = _evaluationTimes(ie)
      payoff.setTime(evaluationTime)
      payoff.eval(quotes)
      val payments = payoff.payments
      if (evaluationTime == exerciseTime) {
        payoff.evalRegressionVariables(quotes, independentVariables(indexExercise))
        var ip = 0
        exerciseFlows(indexExercise) = 0
        while (ip < payments.length) {
          val p = payments(ip)
          val df = getDiscountFactor(p.paymentTime)
          if (p.isExercise) {
            exerciseFlows(indexExercise) += p.amount*df
          } else {
            sum += (p.amount * df)
          }
          ip += 1
        }
        if (indexExercise < _exerciseTimes.length - 1) {
          indexExercise += 1
          exerciseTime = _exerciseTimes(indexExercise)
        }
      } else {
        var ip = 0
        while (ip < payments.length) {
          val p = payments(ip)
          val df = getDiscountFactor(p.paymentTime)
          sum += (p.amount * df)
          ip += 1
        }
      }
      ie += 1
    }
    return sum
  }
}

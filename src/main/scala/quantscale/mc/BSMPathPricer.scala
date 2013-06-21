package quantscale.mc

import quantscale.mc.payoff.MCPayoff

class BSMPathPricer(payoff: MCPayoff, discounter: Discounter) extends PathPricer {

  private var _evaluationTimes = payoff.evaluationTimes

  private var _dfCache: Map[Double, Double] = Map()

  initDiscountFactors(payoff.paymentTimes)

  def initDiscountFactors(paymentTimes: Array[Double]) {
    for (time <- paymentTimes) {
      val df = discounter.df(time)
      _dfCache += (time -> df)
    }
  }

  def evaluationTimes = _evaluationTimes

  private def getDiscountFactor(time: Double): Double = {
    val dfO = _dfCache.get(time)
    if (dfO.isEmpty) {
      val df = discounter.df(time)
      _dfCache += (time -> df)
      return df
    } else {
      return dfO.get
    }
  }

  def eval(path: MCPath): Double = {
    var sum = 0d
    var ie = 0
    payoff.reset()
    while (ie < evaluationTimes.length) {
      val quotes = path.values(ie)
      val evaluationTime = evaluationTimes(ie)
      payoff.setTime(evaluationTime)
      payoff.eval(quotes)
      val payments = payoff.payments
      var ip = 0
      while (ip < payments.length) {
        val p = payments(ip)
        val df = getDiscountFactor(p.paymentTime)
        sum += (p.amount * df)
        ip += 1
      }
      ie += 1
    }
    return sum
  }
}

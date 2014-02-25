package quantscale.mc

import quantscale.mc.payoff.{AdjointMCPayoff, MCPayoff}

class BSMPathPricer(payoff: MCPayoff, discounter: Discounter) extends PathPricer {

  private val _evaluationTimes = payoff.evaluationTimes

  private var _dfCache: Map[Double, Double] = Map()

  private val _payoffA: AdjointMCPayoff = if (payoff.isInstanceOf[AdjointMCPayoff]) payoff.asInstanceOf[AdjointMCPayoff] else null

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

  def evalAdjoint(path: AdjointMCPath): Double = {
    var sum = 0d
    var ie = 0
    _payoffA.reset()
    while (ie < evaluationTimes.length) {
      val quotes = path.values(ie)
      val quotes_b = path.adjointValues(ie)
      val evaluationTime = evaluationTimes(ie)
      _payoffA.setTime(evaluationTime)
      _payoffA.evalAdjoint(quotes, quotes_b)
      val payments = _payoffA.payments
      var ip = 0
      while (ip < payments.length) {
        val p = payments(ip)
        val df = getDiscountFactor(p.paymentTime)
        sum += (p.amount * df)
        var ib = 0
        //FIXME only works for a single payment
        while (ib < quotes_b.length) {
          quotes_b(ib) *= df
          ib += 1
        }

        ip += 1
      }


      ie += 1
    }
    return sum
  }
}

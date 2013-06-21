package quantscale.mc

trait BSM1FMCSpec {
  def initialPrice: Double
  def variance(time: Double): Double
  def driftCf(time: Double): Double
}

class ConstantBSM1FMCSpec(private val _initialPrice : Double, mu: Double, vol: Double) extends BSM1FMCSpec {
  def initialPrice = _initialPrice
  def variance(time: Double): Double = vol*vol*time
  def driftCf(time: Double): Double = math.exp(mu*time)
}
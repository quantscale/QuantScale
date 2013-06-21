package quantscale.analytic

import quantscale.fdm.Epsilon

object BachelierVanillaEuropean {

  def price(isCall: Boolean, strike: Double, forward: Double, normalVol: Double, tte: Double): Double = {
    val sign = if (isCall) 1 else -1
    val sqrttte = math.sqrt(tte)
    val sqrtvar = normalVol * sqrttte
    val d = sign * (forward - strike) / (sqrtvar)
    if (math.abs(forward - strike) < Epsilon.MACHINE_EPSILON_SQRT) {
      return sqrtvar
    } else {
      val Nd = CumulativeNormalDistribution.value(d)
      val nd = NormalDistribution.value(d)

      val value = sign * (forward - strike) * Nd + sqrtvar * nd
      return value
    }
  }
}
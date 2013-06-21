package quantscale.analytic

import quantscale.fdm.Epsilon

object CoxIngersollRossBond {

  //a*(b-r)+sigma sqrt(r) dW
  def price(a: Double, b: Double, sigma: Double, r0: Double, tte: Double): Double = {
    val gamma = math.sqrt(a * a + 2 * sigma * sigma)
    val denom = 1.0 / ((gamma + a) * (math.exp(gamma * tte) - 1) + 2 * gamma)
    val B = 2 * (math.exp(gamma * tte) - 1) * denom
    val A = math.pow(2 * gamma * math.exp((a + gamma) * tte * 0.5) * denom, 2 * a * b / (sigma * sigma))
    return A * math.exp(-B * r0)
  }

  /**
   * @param T: option expiry
   * @param S: bond maturity (S >= T)
   */
  def priceVanillaOption(isCall: Boolean, strike: Double, a: Double, b: Double, sigma: Double, r0: Double, T: Double, S: Double): Double = {
    val gamma = math.sqrt(a * a + 2 * sigma * sigma)
 
    val tte = S - T
    val denom = 1.0 / ((gamma + a) * (math.exp(gamma * tte) - 1) + 2 * gamma)
    val B = 2 * (math.exp(gamma * tte) - 1) * denom
    val A = math.pow(2 * gamma * math.exp((a + gamma) * tte * 0.5) * denom, 2 * a * b / (sigma * sigma))
    val rho = 2 * gamma / (sigma * sigma * (math.exp(gamma * T) - 1.0))
    val psi = (a + gamma) / (sigma * sigma)
    val rBar = math.log(A / strike) / B
    val chiSquare1 = NonCentralChiSquaredDistribution.cdf(2 * rBar * (rho + psi + B), 4 * a * b / (sigma * sigma), 2 * rho * rho * r0 * math.exp(gamma * T) / (rho + psi + B),errtol=1e-12)
    val chiSquare2 = NonCentralChiSquaredDistribution.cdf(2 * rBar * (rho + psi), 4 * a * b / (sigma * sigma), 2 * rho * rho * r0 * math.exp(gamma * T) / (rho + psi),errtol=1e-12)
    val discountS = price(a, b, sigma, r0, S)
    if (T < Epsilon.MACHINE_EPSILON) {
      return if (isCall) math.max(discountS-strike,0.0) else math.max(strike-discountS, 0)
    }
    val discountT = price(a, b, sigma, r0, T)
    val callPrice = discountS * chiSquare1 - strike * discountT * chiSquare2
    return if (isCall) callPrice else callPrice - discountS + strike * discountT 
  }
}
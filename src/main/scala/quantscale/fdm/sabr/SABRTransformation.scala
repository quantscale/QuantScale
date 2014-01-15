package quantscale.fdm.sabr

import quantscale.math.BijectiveFunction1D
import quantscale.analytic.SABRModelSpec

class SABRTransformation(spec: SABRModelSpec, forward: Double) extends BijectiveFunction1D {
  private val fonebeta = math.pow(forward + spec.b, 1 - spec.beta)

  def value(z: Double): Double = {
    val nu = spec.nu
    val rho = spec.rho
    val alpha = spec.alpha
    val beta = spec.beta
    val e = math.exp(nu * z)
    val sh = 0.5 * (e - 1 / e)
    val ch = 0.5 * (e + 1 / e)
    val y = alpha / nu * (sh + rho * (ch - 1))
    val by = (1 - beta) * y
    val f = if ((fonebeta + by) <= 0) 0 else math.pow(fonebeta + by, 1 / (1 - beta))
    return f
    //k(x)
  }

  def inverseValue(k: Double): Double = {
    val nu = spec.nu
    val rho = spec.rho
    val alpha = spec.alpha
    val beta = spec.beta
    val b = spec.b
    val y = (math.pow(k + b, 1 - beta) - fonebeta) / (1 - beta)
    val temp = (rho + nu * y / alpha)
    val z = -1 / spec.nu * math.log((math.sqrt(1 - rho * rho + temp * temp) - rho - nu * y / alpha) / (1 - spec.rho))
    return z
    //x(k)
  }
}

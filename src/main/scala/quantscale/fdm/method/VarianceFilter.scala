package quantscale.fdm.method

import quantscale.fdm.Epsilon

abstract class VarianceFilter {
  def filter(variance: Double, muMod: Double, dxj: Double, dxjplus: Double): Double
}

class UpWindingVarianceFilter extends VarianceFilter {
  override def filter(variance: Double, muMod: Double, dxj: Double, dxjplus: Double): Double = {
    var correction = 0.
    val mu = math.abs(muMod)
    if (mu > Epsilon.MACHINE_EPSILON_FOURTH_ROOT) {
      val z = mu / variance * (dxjplus + dxj) * 0.5 //peclet number
      correction = z
    }
    return variance * (1.0 + correction) //introduce artificial diffusion of 2nd order (see Fusai)
  }  
}
class ScharfetterGummelVarianceFilter extends VarianceFilter {

  override def filter(variance: Double, muMod: Double, dxj: Double, dxjplus: Double): Double = {
    var correction = 0.
    val mu = math.abs(muMod)
    if (mu > Epsilon.MACHINE_EPSILON_FOURTH_ROOT) {
      val z = mu / variance * (dxjplus + dxj) * 0.5 //peclet number
      correction = z - 1.0 + 2 * z / (math.exp(2 * z) - 1)
    }
    return variance * (1.0 + correction) //introduce artificial diffusion of 2nd order (see Fusai)
  }
}

class ExponentialFittingFilter extends VarianceFilter {
  //same thing as ScharfetterGummel expressed differently - for some reasons it is quite slower
  override def filter(variance: Double, muMod: Double, dxj: Double, dxjplus: Double): Double = {
    if (math.abs(muMod) > Epsilon.MACHINE_EPSILON_FOURTH_ROOT) {
      val z = muMod * (dxjplus + dxj) * 0.5
      return z / math.tanh(z / variance)
    } else {
      return variance
    }
  }
}
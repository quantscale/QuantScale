package quantscale.fdm

import quantscale.fdm.mesh.Mesh2D
import quantscale.analytic.SABRModelSpec
import scala.collection.mutable

class SABRLocalVolFDSpec(
                          gridV: Mesh2D,
                          private val spec: SABRModelSpec,
                          private val forward: Double,
                          private val tte: Double,
                          driftV: Double,
                          discountRateV: Double)
  extends Parabolic1DFDSpec(gridV) {

  private val forwardonebeta = math.pow(forward + spec.b, 1 - spec.beta)

  private val factorC = new mutable.HashMap[Double, Double]()
  private val gammaC = new mutable.HashMap[Double, Double]()
  initCache()

  private def initCache() {
    var i = gridV.spaceSize - 1
    val space = gridV.spaceVector
    val efactor = spec.rho * spec.nu * spec.alpha
    val C0 = math.pow(forward + spec.b, spec.beta)
    while (i >= 0) {
      val f = space(i)
      var gamma = 0.0
      var factor = 0.0
      if (f > Epsilon.MACHINE_EPSILON_SQRT) {
        val C = math.pow(f + spec.b, spec.beta)
        val fonebeta = f / C
        val z = (fonebeta - forwardonebeta) / (spec.alpha * (1 - spec.beta))
        factor = C * C * spec.alpha * spec.alpha * (1 + 2 * spec.rho * spec.nu * z + spec.nu * spec.nu * z * z)
        gamma = if (math.abs(f - forward) < Epsilon.MACHINE_EPSILON_SQRT) spec.beta / forwardonebeta else (C - C0) / (f - forward)
        gamma *= efactor
      }
      gammaC(f) = gamma
      factorC(f) = factor
      i -= 1
    }
  }

  def copy(): Parabolic1DFDSpec = {
    return new SABRLocalVolFDSpec(gridV.copy(), spec, forward, tte, driftV, discountRateV)
  }

  def volSquare(timeEnd: Double, dt: Double, f: Double): Double = {
    //tte-timeEnd for forward, is it ok for backward?
    val Dsquare = factorC(f) * math.exp(gammaC(f) * (tte - timeEnd))
    return Dsquare
  }

  def a(timeEnd: Double, dt: Double, space: Double): Double = {
    val v = volSquare(timeEnd, dt, space)
    return 0.5 * v
  }

  def b(timeEnd: Double, dt: Double, space: Double): Double = driftV * space

  def c(timeEnd: Double, dt: Double, space: Double): Double = -discountRateV

  def aIsStateDependent(): Boolean = true

  def bIsStateDependent(): Boolean = true

  def cIsStateDependent(): Boolean = false

}

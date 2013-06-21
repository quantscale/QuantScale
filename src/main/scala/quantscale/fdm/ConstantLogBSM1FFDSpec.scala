package quantscale.fdm
import quantscale.fdm.mesh.Mesh2D

class ConstantLogBSM1FFDSpec(
  gridV: Mesh2D,
  private val volV: Double,
  driftV: Double,
  discountRateV: Double)
  extends Parabolic1DFDSpec(gridV) {

  def a(timeEnd: Double, dt: Double, space: Double): Double = 0.5 * volV * volV
  def b(timeEnd: Double, dt: Double, space: Double): Double = driftV - 0.5 * volV * volV
  def c(timeEnd: Double, dt: Double, space: Double): Double = -discountRateV

  def aIsStateDependent(): Boolean = false
  def bIsStateDependent(): Boolean = false
  def cIsStateDependent(): Boolean = false

  def copy(): Parabolic1DFDSpec = {
    return new ConstantLogBSM1FFDSpec(gridV.copy(), volV, driftV, discountRateV)
  }
}
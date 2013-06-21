package quantscale.fdm
import quantscale.fdm.mesh.Mesh2D

class ConstantBSM1FFDSpec(
  gridV: Mesh2D,
  private val volV: Double,
  driftV: Double,
  discountRateV: Double)
  extends Parabolic1DFDSpec(gridV) {

  def copy() : Parabolic1DFDSpec = {
    return new ConstantBSM1FFDSpec(gridV.copy(), volV, driftV, discountRateV)
  }
  
  def vol(timeEnd: Double, dt: Double, space: Double) = volV

  def a(timeEnd: Double, dt: Double, space: Double): Double = {
    val v = vol(timeEnd, dt, space)
    return 0.5*v*v*space*space
  } 
  def b(timeEnd: Double, dt: Double, space: Double): Double = driftV*space
  def c(timeEnd: Double, dt: Double, space: Double): Double = -discountRateV

  def aIsStateDependent(): Boolean = true
  def bIsStateDependent(): Boolean = true
  def cIsStateDependent(): Boolean = false

}
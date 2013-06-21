package quantscale.fdm
import quantscale.fdm.mesh.Mesh2D

class UncertainBSM1FFDSpec(
  gridV: Mesh2D,
  val volMin: Double,
  val volMax: Double,
  driftV: Double,
  discountRateV: Double,
  val bestCaseLong: Boolean)
  extends ConstantBSM1FFDSpec(gridV, 0.0, driftV, discountRateV) {

  //initial guess?
  override def vol(timeEnd: Double, dt: Double, space: Double): Double = 0.5 * (volMin + volMax);

  def vol(gamma: Double): Double = {
    if (bestCaseLong) {
      return volBestCaseLong(gamma)
    } else {
      return volWorstCaseLong(gamma)
    }
  }

  def volBestCaseLong(gamma: Double): Double = {
    if (gamma > 0) {
      return volMax;
    } else {
      return volMin;
    }
  }

  def volWorstCaseLong(gamma: Double): Double = {
    if (gamma > 0) {
      return volMin;
    } else {
      return volMax;
    }
  }
}
package quantscale.fdm

import quantscale.fdm.mesh.Mesh2D

class VasicekSpec(gridV: Mesh2D, val k: Double, val theta: Double, val sigma: Double) extends Parabolic1DFDSpec(gridV)  {

  def a(timeEnd: Double, dt: Double, space: Double): Double = {
    return sigma * sigma * 0.5
  }

  def aIsStateDependent(): Boolean = true

  def b(timeEnd: Double, dt: Double, space: Double): Double = {
    return k *(theta- space)
  }
  
  def bIsStateDependent(): Boolean = true

  def c(timeEnd: Double, dt: Double, space: Double): Double = {
    return -space
  }

  def cIsStateDependent(): Boolean = true

    def copy() : Parabolic1DFDSpec = {
    return new VasicekSpec(gridV.copy(), k, theta, sigma)
  }
}
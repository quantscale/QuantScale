package quantscale.fdm

import quantscale.fdm.mesh.Mesh2D

class ConstantCIRSpec(gridV: Mesh2D, val aCIR: Double, val bCIR: Double, val sigmaCIR: Double) extends Parabolic1DFDSpec(gridV) {

  def a(timeEnd: Double, dt: Double, space: Double): Double = {
    return sigmaCIR * sigmaCIR * space * 0.5
  }

  def aIsStateDependent(): Boolean = true

  def b(timeEnd: Double, dt: Double, space: Double): Double = {
    return aCIR - bCIR * space
  }
  def bIsStateDependent(): Boolean = true

  def c(timeEnd: Double, dt: Double, space: Double): Double = {
    return -space
  }

  def cIsStateDependent(): Boolean = true

  def copy() : Parabolic1DFDSpec = {
    return new ConstantCIRSpec(gridV.copy(), aCIR, bCIR, sigmaCIR)
  }
}
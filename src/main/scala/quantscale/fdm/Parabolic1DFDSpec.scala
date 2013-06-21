package quantscale.fdm
import quantscale.fdm.mesh.Mesh2D

/**
 * 
 * dV/dt + a*d2V/dS2 + b*dV/dS + c*V=0  
 * Represents process defined by dS = drift(t,t+dt,S) dt + vol(t,t+dt,S) dW
 * Use mu*S, vol*S for a lognormal process.
 * Use mu-0.5*vol*vol, vol for a log process with normal distribution.
 */
abstract class Parabolic1DFDSpec(gridV: Mesh2D) {
  val grid: Mesh2D = gridV

  def copy() : Parabolic1DFDSpec
   
  //discretized rates, eventually.

  def a(timeEnd: Double, dt: Double, space: Double): Double

  def aIsStateDependent(): Boolean
  def bIsStateDependent(): Boolean
  def cIsStateDependent(): Boolean
  
  def b(timeEnd: Double, dt: Double, space: Double): Double

  def c(timeEnd: Double, dt: Double, space: Double): Double
}
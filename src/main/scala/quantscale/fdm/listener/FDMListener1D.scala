package quantscale.fdm.listener

import quantscale.fdm.State

abstract class FDMListener1D {

  def update(t : Double, f : State, x : Array[Double]) : Unit
}
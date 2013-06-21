package quantscale.fdm.listener
import quantscale.fdm.Epsilon
import quantscale.fdm.State

class GammaListener(var observationTime: Array[Double]) extends FDMListener1D {
  var gamma: Array[Double] = null
  var delta: Array[Double] = null
  var x: Array[Double] = null

  override def update(t: Double, f: State, x: Array[Double]): Unit = {
    if (isObservationTime(t)) {
      var i = 1
      delta = new Array[Double](f.price.length)
      gamma = new Array[Double](f.price.length)
      if (this.x == null) this.x = x
      while (i < f.price.length - 1) {
        val hi = x(i+1)-x(i)
        val him = x(i)-x(i-1)
        val hiOverHim = hi/him
        var denom = 1.0/ (hi*him*(1+hiOverHim))
        gamma(i) = 2*hiOverHim*denom*f.price(i-1) - 2*(1+hiOverHim)*denom*f.price(i)+2*denom*f.price(i+1)
        val hiOverHimSq = hiOverHim*hiOverHim
        denom = denom*him
        delta(i) =  -hiOverHimSq*denom*f.price(i-1)-(1 - hiOverHimSq) * denom*f.price(i)+denom*f.price(i+1)
//          (f.price(i + 1) - 2 * f.price(i) + f.price(i - 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)))
        i += 1
      }
    }
  }

  private def isObservationTime(t: Double): Boolean = {
    var i = observationTime.length - 1
    while (i >= 0) {
//      println("time="+t + " obs "+observationTime(i))
      if (math.abs(t - observationTime(i)) < Epsilon.MACHINE_EPSILON_SQRT) {
        return true
      }
      i -= 1
    }
    return false
  }

  def print() {
    var i = 1
    println("Spot\tGamma");
    if (gamma == null) {
      println("Gamma undefined")
    } else {
      while (i < gamma.length - 1) {
        println(x(i) + "\t" + gamma(i))
        i += 1
      }
    }
  }
}
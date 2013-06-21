package quantscale.math

trait Extrapolation1D {
  def init(x: Array[Double], y: Array[Double])
  def isInExtrapolationRange(z: Double): Boolean
  def value(z: Double): Double
}

class LinearExtrapolation extends Extrapolation1D {
  private var x, y: Array[Double] = null

  def init(x: Array[Double], y: Array[Double]) {
    this.x = x
    this.y = y
  }

  def isInExtrapolationRange(z: Double): Boolean = {
    return z < x(0) || z > x(x.length - 1)
  }

  def value(z: Double): Double = {
    if (z < x(0)) {
      return y(0) + (y(1) - y(0)) / (x(1) - x(0)) * (z - x(0))
    } else if (z > x(x.length - 1)) {
      return y(x.length - 1) + (y(x.length - 2) - y(x.length - 1)) / (x(x.length - 2) - x(x.length - 1)) * (z - x(x.length - 1))
    } else {
      return Double.NaN
    }
  }
}
package quantscale.math
import java.lang.Math

class JacobianConcentrationSingle(var alpha: Double, var yStar: Double, var A:Double) extends Function2D {

  def value(x: Double, y: Double): Double = {
    val z = y - yStar
    return A*Math.sqrt(alpha * alpha + z * z)
  }

  def integralValue(x: Double, Smin: Double, Smax: Double): Double = {
    return yStar + alpha * Math.sinh((1.0-x) * asinh((Smin - yStar) / alpha) +
      x * asinh((Smax - yStar) / alpha))
  }

  def asinh(x: Double): Double = {
    return Math.log(x + Math.pow(x * x + 1, 0.5));
  }

}
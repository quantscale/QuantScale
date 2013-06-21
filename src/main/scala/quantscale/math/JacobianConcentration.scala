package quantscale.math
import java.lang.Math

class JacobianConcentration(var alpha: Double, var yStar: Array[Double], var A:Double) extends Function2D {

  def singleValue(y: Double, yStarSingle: Double) : Double= {
      val z = y - yStarSingle
    return Math.sqrt(alpha * alpha + z * z)

  }
  def value(x: Double, y: Double): Double = {
    var i = yStar.length-1
    var sum = 0.0
    while (i>=0) {
      val z =  1.0/singleValue(y, yStar(i))
      sum += z*z
      i -= 1
    }
    return A/Math.sqrt(sum)
  }

}
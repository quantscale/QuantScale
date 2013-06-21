package quantscale.fdm.payoff
import quantscale.math.CubicSpline
import scala.collection.mutable.ArrayBuffer
import quantscale.math.CubicPP

class BoundaryListener() {
  private var timeSize = 0
  var boundary = new ArrayBuffer[Double](50) // exercise boundary for time t_i
  var t = new ArrayBuffer[Double](50)

  def init(timeSize: Int) {
    boundary.sizeHint(timeSize)
    t.sizeHint(timeSize)
  }
  def setValue(time: Double, space: Double, state: Double) {
    boundary += space
    t += time
  }

  def makeSpline(): CubicPP = {
    val derivative = new Array[Double](t.size)
    val tArray = t.toArray
    val bArray = boundary.toArray
    CubicSpline.computeHarmonicFirstDerivativePCHIM(tArray, bArray, derivative)
    return CubicSpline.makeHermiteSpline(tArray, bArray, derivative)
  }
}
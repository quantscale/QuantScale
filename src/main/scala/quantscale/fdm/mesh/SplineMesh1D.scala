package quantscale.fdm.mesh

import scala.collection.mutable.ArrayBuffer
import quantscale.math.CubicSpline
import scala.util.Sorting
import quantscale.fdm.Epsilon


class SMPoint(override val value: Double, override val isMiddle: Boolean = false, val width: Double = 0, val slope: Double = 0.1)
  extends Point(value: Double, isMiddle: Boolean) {

}

class SplineMesh1D(
                    private var _size: Int,
                    val boundaries: Mesh1DBoundaries,
                    var specialPoints: Array[SMPoint]) extends Mesh1D {

  private var _x: Array[Double] = null

  init()
  build()

  override def x = _x

  override def size = _x.length

  private def init() {
    if (this.specialPoints == null) {
      this.specialPoints = new Array[SMPoint](0);
    } else {
      this.specialPoints = Point.sortAndRemoveIdenticalPoints(this.specialPoints);
      _size = math.max(_size, specialPoints.length + 1)
    }
  }

  private def build() {
    val uPrime = new ArrayBuffer[Double](specialPoints.length);
    val uStar = new ArrayBuffer[Double](specialPoints.length);
    uPrime += (boundaries.min);
    uStar += (boundaries.min);
    val delta = (boundaries.max - boundaries.min) / (_size - 1);
    for (i <- 0 until specialPoints.length) {
      val k = math.round((specialPoints(i).value - boundaries.min) / delta).toInt;
      var pPrime = boundaries.min + delta * k;
      if (specialPoints(i).isMiddle) {
        pPrime += 0.5 * delta;
      }
      if (specialPoints(i).width > 0 && specialPoints(i).slope > 0) {
        uPrime += (pPrime - specialPoints(i).width);
        uStar += (specialPoints(i).value - specialPoints(i).width * specialPoints(i).slope);
        uPrime += (pPrime);
        uStar += (specialPoints(i).value);
        uPrime += (pPrime + specialPoints(i).width);
        uStar += (specialPoints(i).value + specialPoints(i).width * specialPoints(i).slope);
      } else {
        uPrime += (pPrime);
        uStar += (specialPoints(i).value);
      }
    }
    if (math.abs(uPrime(uPrime.length - 1) - boundaries.max) > Epsilon.MACHINE_EPSILON_SQRT) {
      uPrime += (boundaries.max);
      uStar += (boundaries.max);
    }
    //TODO check that xmin and xmax are not already included
    val uPrimeArray = uPrime.toArray
    val uStarArray = uStar.toArray
    val derivatives = new Array[Double](uPrime.length)
    //    CubicSpline.computeHarmonicFirstDerivativePCHIM(uPrimeArray, uStarArray, derivatives)
    CubicSpline.computeC2FirstDerivative(uPrimeArray, uStarArray, derivatives)
    val interpolator = CubicSpline.makeHermiteSpline(uPrimeArray, uStarArray, derivatives)
    //    val interpolator = CubicSpline.makeCubicSpline(uPrimeArray, uStarArray)
    _x = new Array[Double](_size)
    for (i <- 0 until _size) {
      _x(i) = interpolator.value(delta * i + boundaries.min);
    }
    _x(0) = boundaries.min;
    _x(_size - 1) = boundaries.max;
  }

  def insertPoints(s: Array[Double]) {
    var buf = new ArrayBuffer[Double]()
    Sorting.quickSort(s)
    var i = 0
    var h = 0
    while (i < _x.length) {
      while (h < s.length && s(h) < _x(i) + Epsilon.MACHINE_EPSILON_SQRT) {
        if (s(h) < _x(i) - Epsilon.MACHINE_EPSILON_SQRT) {
          buf += s(h)
        }
        h += 1
      }
      buf += _x(i)
      i += 1
    }
    //add points after xmax
    while (h < s.length) {
      buf += s(h)
      h += 1
    }
    _x = buf.toArray
  }
}
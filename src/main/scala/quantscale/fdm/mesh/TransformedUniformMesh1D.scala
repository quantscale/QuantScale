package quantscale.fdm.mesh

import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting
import quantscale.fdm.Epsilon
import quantscale.math.BijectiveFunction1D

/**
 * uniform input between y_min and y_max is going to be transformed
 * to x_min to x_max.
 *
 * @param transformation a bijection
 * @param _size
 * @param boundaries boundaries defined in terms of the uniform input y
 * @param specificPoint point defined in terms of the input y, use transformation.inverse
 */
class TransformedUniformMesh1D(val transformation: BijectiveFunction1D, private val _size: Int, boundaries: Mesh1DBoundaries, val specificPoint: Point, pinMinAndMax: Boolean = true) extends Mesh1D {
  private var _x: Array[Double] = null
  var delta = 0.0
  build()

  override def size = _x.length

  private def build() {
    var xmax = boundaries.max
    var xmin = boundaries.min

    delta = (xmax - xmin) / (_size - 1);
    if (specificPoint == null) {
      buildTransformedUniform(xmin, xmax, delta, _size)
    } else {
      val xstar = specificPoint.value
      xmin = math.min(xstar, xmin)
      xmax = math.max(xstar, xmax)

      val j0: Int = math.floor((xstar - xmin) / delta).toInt;
      val xj0 = xmin + j0 * delta
      if (specificPoint.isMiddle) {
        val mid = 0.5 * (transformation.value(xj0) + transformation.value(xj0 + delta))
        val xmid = transformation.inverseValue(mid);
        val d = xstar - xmid
        xmax = xmax + d
        xmin = xmin + d
        if (d > Epsilon.MACHINE_EPSILON_SQRT) {
          //xmin -= delta
          //add one point
        } else {
          //xmax += delta
          //add one point
        }
      } else {
        val d = xstar - xj0
        if (math.abs(d) > Epsilon.MACHINE_EPSILON_SQRT) {
          xmax = xmax + d
          xmin = xmin + d
        }
      }
      buildTransformedUniform(xmin, xmax, delta, _size)
    }
  }

  private def buildTransformedUniform(xmin: Double, xmax: Double, delta: Double, size: Int) {
    _x = new Array[Double](size)
    if (pinMinAndMax) {
      _x(0) = transformation.value(boundaries.min)
    } else {
      _x(0) = transformation.value(xmin)
    }
    var i = 1
    while (i < size - 1) {
      _x(i) = transformation.value(xmin + i * delta)
      i += 1
    }
    if (pinMinAndMax) {
      _x(size - 1) = transformation.value(boundaries.max)
    } else {
      _x(size - 1) = transformation.value(xmin + (size - 1) * delta)
    }

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

  override def x: Array[Double] = _x
}
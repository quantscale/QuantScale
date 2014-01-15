package quantscale.fdm.mesh

import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting
import quantscale.fdm.Epsilon

class UniformMesh1D(private var _size: Int, boundaries: Mesh1DBoundaries, val dx: Double, specialPoint: Double = Double.NaN, isMiddle: Boolean = false) extends Mesh1D {
  private var _x: Array[Double] = null

  override def size(): Int = {
    return _x.length
  }

  def this(_size: Int, boundaries: Mesh1DBoundaries) {
    this(_size, boundaries, (boundaries.max - boundaries.min) / (_size - 1))
  }


  def this(_size: Int, boundaries: Mesh1DBoundaries, specialPoint: Double, isMiddle: Boolean) {
    this(_size, boundaries, (boundaries.max - boundaries.min) / (_size - 1), specialPoint, isMiddle)
  }

  build()

  private def build() {
    var min = boundaries.min
    var max = boundaries.max
    if (!specialPoint.isNaN) {
      val j0 = math.floor((specialPoint - min) / dx).toInt
      if (isMiddle) {
        val xj0 = min + dx * j0
        val mid = (xj0 + xj0 + dx) / 2
        val d = specialPoint - mid
        if (d > Epsilon.MACHINE_EPSILON_SQRT) {
          max = max + d
          min = min + d - dx
          _size += 1
        } else if (d < -Epsilon.MACHINE_EPSILON_SQRT) {
          max = max + d + dx
          min = min + d
          _size += 1
        }
      } else {
        val d = specialPoint - (min + dx * j0)
        if (math.abs(d) > Epsilon.MACHINE_EPSILON_SQRT) {
          max = max + d
          min = min + d - dx
          _size += 1
        }
      }
    }
    var i = 1
    _x = new Array[Double](_size)
    _x(0) = min
    while (i < _size - 1) {
      _x(i) = _x(i - 1) + dx
      i += 1
    }
    _x(0) = boundaries.min
    _x(_size - 1) = boundaries.max
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
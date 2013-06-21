package quantscale.fdm.mesh
import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting
import quantscale.fdm.Epsilon

class UniformMesh1D(private val _size: Int, boundaries: Mesh1DBoundaries, val dx: Double) extends Mesh1D {
  private var _x = new Array[Double](_size)
  override def size() : Int = { return _x.length }

  def this(_size: Int, boundaries: Mesh1DBoundaries) {
      this(_size, boundaries, (boundaries.max - boundaries.min) / (_size - 1))
  }


  build()

  private def build() {
    var i = 1
    _x(0) = boundaries.min
    while (i < _size - 1) {
      _x(i) = _x(i - 1) + dx
      i += 1
    }
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
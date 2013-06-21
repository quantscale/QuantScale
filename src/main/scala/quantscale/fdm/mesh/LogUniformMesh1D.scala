package quantscale.fdm.mesh
import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting
import quantscale.fdm.Epsilon

class LogUniformMesh1D(private val _size: Int, boundaries: Mesh1DBoundaries, val specificPoint: Double) extends Mesh1D {
  private var _x = new Array[Double](_size)
  var logdx = 0.0

  override def size = _x.length
  
  build()

  private def build() {
    var xmax = math.log(boundaries.max)
    var xmin = math.log(boundaries.min)

    val xstar = math.log(specificPoint)
    
    xmin = math.min(xstar, xmin)
    xmax = math.max(xstar, xmax)
    logdx = (xmax - xmin) / (_size - 1);
    if (xstar-xmin > xmax-xstar) {
      val strikeIndex: Int = math.round((xstar - xmin) / logdx).toInt;
      logdx = (xstar - xmin) / strikeIndex;
      val dx = math.exp(logdx)
      var i = 1
      _x(0) = math.exp(xmin)
      while (i <= _size - 1) {
        _x(i) = _x(i - 1) * dx
        i += 1
      }
      _x(strikeIndex) = specificPoint;
    } else {
      val strikeIndex: Int = math.round((xmax - xstar) / logdx).toInt;
      logdx = (xstar - xmax) / strikeIndex;
      val dx = math.exp(logdx)
      _x(_size-1) = math.exp(xmax)
      var i = _size -2
      while (i >= 0) {
        _x(i) = _x(i + 1) * dx
        i -= 1
      }
      _x(size-1-strikeIndex) = specificPoint;
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
    while (h<s.length) {
      buf+=s(h)
      h+=1
    }
    _x = buf.toArray
  }

  override def x: Array[Double] = _x
}
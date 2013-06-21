package quantscale.fdm.mesh {
  import scala.collection.mutable.ArrayBuffer
import java.util.Arrays
import quantscale.fdm.Epsilon
import scala.util.Sorting

  class UniformMesh2D(spaceSizeV: Int, timeSize: Int, bounds: MeshBoundaries, specialPoint: Double) extends Mesh2D {
    private val x = new Array[Double](spaceSizeV);
    private var t = new Array[Double](timeSize + 1);

    private var dxV = 0.0;

    initTime();
    initSpace();

    override def copy() : Mesh2D = {
      return new UniformMesh2D(spaceSizeV, timeSize, bounds, specialPoint)
    }
    
    private def initTime() {
      var i = t.length - 1;
      val dt = bounds.lastTime / timeSize;
      t(i) = 0;
      while (i > 1) {
        t(i - 1) = t(i) + dt;
        i -= 1;
      }
      t(0) = bounds.lastTime;
    }

    def insertTime(s: Array[Double]) {
      var buf = new ArrayBuffer[Double]()
      Sorting.quickSort(s)
      var i = 0
      var h = s.length - 1
      while (i < t.length) {
        //first t is tte
        while (h >= 0 && s(h) > t(i) - Epsilon.MACHINE_EPSILON_SQRT) {
          if (s(h) > t(i) + Epsilon.MACHINE_EPSILON_SQRT) {
            buf += s(h)
          }
          //else it is same as ti
          h -= 1
        }
        buf += t(i)
        i += 1
      }
      t = buf.toArray
    }

    private def initSpace() {
      val xMax = bounds.topSpace;
      val xMin = bounds.bottomSpace;
      var i = 0;
      dxV = (xMax - xMin) / spaceSize;
      val strikeIndex: Int = math.round((specialPoint - xMin) / dxV).toInt;
      dxV = (specialPoint - xMin) / strikeIndex;
      x(0) = xMin;
      while (i < spaceSize - 1) {
        x(i + 1) = x(i) + dxV;
        i += 1;
      }
      x(strikeIndex) = specialPoint;
    }

    def spaceInterval = dxV;
    override def spaceVector = x;
    override def lastTime = bounds.lastTime;
    override def spaceSize = spaceSizeV;
    override def timeIterator = t.iterator;
  }
}
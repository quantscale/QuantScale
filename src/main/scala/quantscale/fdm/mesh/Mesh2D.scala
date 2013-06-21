package quantscale.fdm.mesh
import java.util.Arrays
import quantscale.fdm.transform.CoordinateTransformation
import quantscale.fdm.transform.IdentityTransformation
import quantscale.fdm.Epsilon
import scala.collection.mutable.ArrayBuffer

abstract class Mesh2D {
  def spaceSize: Int;
  def spaceVector: Array[Double];
  def lastTime: Double;
  def timeIterator: Iterator[Double];

  def copy() : Mesh2D
  
  var spaceTransform : CoordinateTransformation = IdentityTransformation;
}

class DividedMesh2D(originalMesh : Mesh2D, divisor : Int) extends Mesh2D {
  private var timeList = MeshUtil.divideTime(originalMesh.timeIterator, divisor)
   
  override def copy() : Mesh2D = {
    return new DividedMesh2D(originalMesh.copy(), divisor)
  }
  override def spaceSize: Int = originalMesh.spaceSize
  def spaceVector: Array[Double] = originalMesh.spaceVector
  def lastTime: Double = originalMesh.lastTime
  def timeIterator: Iterator[Double] = timeList.iterator 
}

object MeshUtil {
  def divideTime(timeIterator:  Iterator[Double], divisor : Int) : ArrayBuffer[Double] = {
    val it = timeIterator  
    var dividedTime = new ArrayBuffer[Double]
    var t0 = it.next()
    var t1 = t0
    dividedTime += t0
      while (it.hasNext) {
        t1 = it.next()
        //t1 < t0
        var j = 1
        while (j <= divisor) {
            dividedTime += t0 + (t1-t0)*j /divisor
            j+=1
        }
        t0 = t1
      }
    return dividedTime
  }
  
  def locateLowerIndex(x: Array[Double], z: Double): Int = {
    if (z < x(0)) {
      return 0;
    }
    if (z > x(x.length - 1)) {
      return x.length - 1;
    }
    return locateLowerIndexBounded(x, z);
  }

  def locateLowerIndexBounded(x: Array[Double], z: Double): Int = {
    // WARNING the array must be sorted first, we assume that.
    val j = Arrays.binarySearch(x, z);
    return if (j < 0) (-j - 2) else j;
  }

  def findIndex(x: Array[Double], spot: Double): Int = {
    var i: Int = locateLowerIndex(x, spot);
    if (i < 0) i = 0;

    // if there is a near exact place on the grid for spot, it is
    // solutionGridPoint or +1.
    val error = Math.abs(x(i) - spot);
    if (error > Epsilon.MACHINE_EPSILON_SQRT && i < x.length - 1) {
      val upperError = Math.abs(x(i + 1) - spot);
      if (upperError < Epsilon.MACHINE_EPSILON_SQRT) {
        i += 1;
      }
    }

    return i;
  }

  def interpolateLinearly(x: Array[Double], y: Array[Double], z: Double, i: Int): Double = {
    if ((i == x.length - 1) || Math.abs(x(i) - z) < Epsilon.MACHINE_EPSILON_SQRT) {
      return y(i);
    }

    val w = (x(i + 1) - x(i));
    if (Math.abs(w) < Epsilon.MACHINE_EPSILON_SQRT) {
      return y(i);

    }
    // linear interpolation for now, we could do cubic later
    return (y(i) * (x(i + 1) - z) + y(i + 1) * (z - x(i))) / w;
  }
}
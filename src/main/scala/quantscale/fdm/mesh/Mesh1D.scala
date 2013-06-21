package quantscale.fdm.mesh

import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting

import quantscale.fdm.Epsilon

abstract class Mesh1D {
  def x: Array[Double]
  def size: Int
}
//    points.sortWith((a,b) => a.compare(b)>0)

object Point {
  def sortAndRemoveIdenticalPoints[T <: Point: ClassManifest](points: Array[T]): Array[T] = {
    Sorting.quickSort(points.asInstanceOf[Array[Point]])
    if (points.length > 0) {
      val l = new ArrayBuffer[T](points.length)
      var previous = points(0)
      l += points(0)
      for (i <- 1 until points.length) {
        if (points(i).value - previous.value > Epsilon.MACHINE_EPSILON_SQRT) {
          l += points(i)
        }
        previous = points(i)
      }
      return l.toArray
    }
    return points
  }
}

class Point(val value: Double, val isMiddle: Boolean) extends Ordered[Point] {
  def compare(that: Point): Int = {
    return math.signum(this.value - that.value).toInt
  }
}
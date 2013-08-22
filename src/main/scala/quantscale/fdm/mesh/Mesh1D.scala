package quantscale.fdm.mesh

import scala.collection.mutable.ArrayBuffer
import scala.util.Sorting

import quantscale.fdm.Epsilon
import scala.reflect.ClassTag

abstract class Mesh1D {
  def x: Array[Double]

  def size: Int
}

//    points.sortWith((a,b) => a.compare(b)>0)

object Point {

  def sortAndRemoveIdenticalPoints[T <: Point : ClassTag](points: Array[T]): Array[T] = {
    if (points.length > 0) {
      Sorting.quickSort(points)(new PointOrdering[T]())
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

class PointOrdering[T <: Point] extends Ordering[T] {
  def compare(p1: T, p2: T): Int = {
    return math.signum(p1.value - p2.value).toInt
  }
}

class Point(val value: Double, val isMiddle: Boolean) extends Ordered[Point] {
  def compare(that: Point): Int = {
    return math.signum(this.value - that.value).toInt
  }
}
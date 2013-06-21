package quantscale.fdm.payoff

import quantscale.fdm.Epsilon
import scala.util.Sorting
import scala.collection.mutable.ArrayBuffer

/**
 * A payoff schedule that will be evaluated on the time 1D mesh from the last element to the first one.
 * The array will be sorted if not already sorted 
 */

//another possibility is to initialize with the mesh, but this would make difficult adaptive meshes

class BackwardSchedule(private val _times: Array[Double], private val eps: Double = Epsilon.MACHINE_EPSILON_SQRT) {
  
  def copy() : BackwardSchedule = {
    return new BackwardSchedule(_times, eps)
  }
  
  val times = removeIdenticalTimes(_times)
  private var _currentIndex = times.length - 1
  private var _isMeshTime: Boolean = false
  private var _previousMeshTime: Double = 0.0

  private def removeIdenticalTimes(t : Array[Double]) : Array[Double] = {
      Sorting.quickSort(t)
      val newTimes = new ArrayBuffer[Double](t.length)
      var current = t(0)
      newTimes += current
      for (i<- 1 until t.length) {
        if (t(i)-current > eps) {
          newTimes += t(i)
        }
        current = t(i)
      }
      return newTimes.toArray
  }
  
  def length(): Int = times.length
  
  def isMeshTime(): Boolean = {
    return _isMeshTime
  }

  def reset() {
    _currentIndex = times.length - 1
  }
  
  def advance(meshTime: Double) {
//    if (_currentIndex < times.length-1 && meshTime > _previousMeshTime) {
//      _currentIndex += 1 //for scheme that evaluates in between mesh
//    }
    while (_currentIndex >= 0 && times(_currentIndex) > meshTime + eps) {
      _currentIndex -= 1
    }
    if (_currentIndex == -1) {
      _isMeshTime = false
    } else if (times(_currentIndex) > meshTime - eps) {
      _isMeshTime = true
    } else {
      _isMeshTime = false
    }
  }
}
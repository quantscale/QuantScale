package quantscale.mc

import scala.collection.mutable.ArrayBuffer

object BrownianBridgeFactory extends BrownianTransformFactory {
  def makeBrownianTransform(t: Array[Double]): BrownianTransform = {
    return new BrownianBridge(t)
  }
}

class BrownianBridge(t: Array[Double]) extends BrownianTransform {
  private val size = t.length
  private val _sqrtdt = Array.ofDim[Double](size)
  private val _bridgeIndex = Array.ofDim[Int](size)
  private val _stdDev = Array.ofDim[Double](size)
  private val _leftWeight = Array.ofDim[Double](size)
  private val _rightWeight = Array.ofDim[Double](size)
  private val _leftIndex = Array.ofDim[Int](size)
  private val _rightIndex = Array.ofDim[Int](size)
  initialize()

  private def initialize() {
    assert(t(0) > 0, "first step must be > 0 but was " + t(0))
    _sqrtdt(0) = math.sqrt(t(0))
    var i = 1
    while (i < size) {
      _sqrtdt(i) = math.sqrt(t(i) - t(i - 1))
      i += 1
    }

    // map is used to indicate which points are already constructed.
    // If map(i) is zero, path point i is yet unconstructed.
    // map(i)-1 is the index of the variate that constructs
    // the path point # i.
    val map = new Array[Double](size)

    //  The first point in the construction is the global step.
    map(size - 1) = 1
    //  The global step is constructed from the first variate.
    _bridgeIndex(0) = size - 1
    //  The variance of the global step
    _stdDev(0) = math.sqrt(t(size - 1))
    //  The global step to the last point in time is special.
    _leftWeight(0) = 0.0
    _rightWeight(0) = 0.0
    var j = 0
    i = 1
    while (i < size) {
      // Find the next unpopulated entry in the map.
      while (map(j) != 0) j += 1
      var k = j;
      // Find the next populated entry in the map from there.
      while (map(k) == 0) k += 1
      // l-1 is now the index of the point to be constructed next.
      var l = j + ((k - 1 - j) >>> 1)
      map(l) = i;
      // The i-th Gaussian variate will be used to set point l-1.
      _bridgeIndex(i) = l
      _leftIndex(i) = j
      _rightIndex(i) = k
      if (j != 0) {
        _leftWeight(i) = (t(k) - t(l)) / (t(k) - t(j - 1));
        _rightWeight(i) = (t(l) - t(j - 1)) / (t(k) - t(j - 1));
        _stdDev(i) =
          math.sqrt(((t(l) - t(j - 1)) * (t(k) - t(l))) / (t(k) - t(j - 1)));
      } else {
        _leftWeight(i) = (t(k) - t(l)) / t(k);
        _rightWeight(i) = t(l) / t(k);
        _stdDev(i) = math.sqrt(t(l) * (t(k) - t(l)) / t(k));
      }
      j = k + 1
      if (j >= size) j = 0 //  wrap around
      i += 1
    }
  }

  def transform(input: Array[Double], dimension: Int,
                output: Array[Array[Double]]) {
    // We use output to store the path...
    //follow Mike Giles: use first numbers for dimensions, then times
    var d = 0
    while (d < dimension) {
      output(size - 1)(d) = _stdDev(0) * input(d);
      d += 1
    }
    var i = 1
    while (i < size) {
      var j = _leftIndex(i);
      var k = _rightIndex(i);
      var l = _bridgeIndex(i);
      d = 0
      while (d < dimension) {
        if (j != 0) {
          output(l)(d) = _leftWeight(i) * output(j - 1)(d) + _rightWeight(i) * output(k)(d) + _stdDev(i) * input(d + i * dimension);
        } else {
          output(l)(d) = _rightWeight(i) * output(k)(d) + _stdDev(i) * input(d + i * dimension)
        }
        d += 1
      }
      i += 1
    }
    // ...after which, we calculate the variations and
    // normalize to unit times
    i = size - 1
    while (i >= 1) {
      d = 0
      while (d < dimension) {
        output(i)(d) -= output(i - 1)(d)
        output(i)(d) /= _sqrtdt(i)
        d += 1
      }
      i -= 1
    }
    d = 0
    while (d < dimension) {
      output(0)(d) /= _sqrtdt(0)
      d+=1
    }
  }
}
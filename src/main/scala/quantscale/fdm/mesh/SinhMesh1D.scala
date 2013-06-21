package quantscale.fdm.mesh
import quantscale.math.CubicSpline
import quantscale.fdm.Epsilon
import quantscale.math.InverseHyperbolics

object SinhMesh1D {
  val ALPHA_ALMOST_UNIFORM = 10.0
  val ALPHA_HIGHLY_NONUNIFORM = 0.1
}

/**
 * Tavella Randall hyperbolic sine concentration of points
 */
class SinhMesh1D(
  val size: Int,
  val boundaries: Mesh1DBoundaries,
  private val pinMinAndMax: Boolean,
  private var specialPoints: Array[Point],
  private val starPoint: Double,
  private val alpha: Double) extends Mesh1D {

  private var _x: Array[Double] = null
  private val min = boundaries.min
  private val max = boundaries.max

  var placeExactly = false

  init()
  build()

  override def x: Array[Double] = _x

  private def init() {
    if (this.specialPoints == null) {
      this.specialPoints = new Array[Point](0);
    } else {
      this.specialPoints = Point.sortAndRemoveIdenticalPoints(this.specialPoints);
    }
  }

  private def build() {
    val u = new Array[Double](size);
    val alphan = alpha * (max - min) // normalize alpha
    val c1 = InverseHyperbolics.asinh((min - starPoint) / alphan)
    val c2 = InverseHyperbolics.asinh((max - starPoint) / alphan)
    _x = new Array[Double](size)
    var i = 0
    while (i < size) {
      u(i) = (i.toDouble) / (size - 1)
      _x(i) = starPoint + alphan * math.sinh(c2 * u(i) + c1 * (1 - u(i)))
      i += 1
    }
    val uStar = new Array[Double](specialPoints.length + 2)
    val uPrime = new Array[Double](specialPoints.length + 2)
    uStar(0) = 0
    uPrime(0) = 0
    uStar(specialPoints.length + 1) = 1.0
    uPrime(specialPoints.length + 1) = 1.0

    i = 0
    while (i < specialPoints.length) {
      val k = findClosestIndexLessThanOrEqual(specialPoints(i).value, _x);
      uStar(i + 1) = u(k) + (u(k + 1) - u(k)) * (specialPoints(i).value - _x(k)) / (_x(k + 1) - _x(k));
      uPrime(i + 1) = (Math.round(uStar(i + 1) * (size - 1)).toDouble) / (size - 1);
      if (specialPoints(i).isMiddle) {
        uPrime(i + 1) += 0.5 / (size - 1);
      }
      i += 1
    }
    val derivatives = new Array[Double](uPrime.length)
    CubicSpline.computeHarmonicFirstDerivativePCHIM(uPrime, uStar, derivatives)
    val interpolator = CubicSpline.makeHermiteSpline(uPrime, uStar, derivatives)
    _x(0) = min
    _x(size - 1) = max
    i = 1
    while (i < size - 1) {
      val z = interpolator.value(u(i))
      _x(i) = starPoint + alphan * math.sinh(c2 * z + c1 * (1 - z))
      i += 1
    }
    if (placeExactly) {
      i = 0
      while (i < specialPoints.length) {
        val k = findClosestIndexLessThanOrEqual(specialPoints(i).value, _x);
        val forward = _x(k + 1) - specialPoints(i).value;
        val backward = specialPoints(i).value - _x(k);
        if (specialPoints(i).isMiddle) {
          _x(k) = specialPoints(i).value - (_x(k + 1) - _x(k)) / 2;
          _x(k + 1) = specialPoints(i).value + (_x(k + 1) - _x(k)) / 2;
        } else {
          if (forward < backward) {
            _x(k + 1) = specialPoints(i).value;
          } else {
            _x(k) = specialPoints(i).value;
          }
        }
        i += 1
      }
    }
  }

  private def findClosestIndexLessThanOrEqual(target: Double, sortedVect: Array[Double]): Int = {
    var index = -1;
    var i = 1
    while (i < sortedVect.length && (target + Epsilon.MACHINE_EPSILON_SQRT > sortedVect(i))) {
      i += 1
    }
    if (i == sortedVect.length || target + Epsilon.MACHINE_EPSILON_SQRT < sortedVect(i)) {
      index = i - 1
    }
    if (index == -1) {
      throw new RuntimeException(target + " not found in array " + sortedVect);
    }
    return index
  }

  def getSpecificPoints(): Array[Double] = {
    val spec = new Array[Double](specialPoints.length)
    var i = 0
    while (i < spec.length) {
      spec(i) = specialPoints(i).value;
      i += 1
    }
    return spec;
  }
}
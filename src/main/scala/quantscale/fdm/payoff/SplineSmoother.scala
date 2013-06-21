package quantscale.fdm.payoff
import quantscale.fdm.Epsilon
import java.util.Arrays
import quantscale.math.CubicSpline
import quantscale.math.SplineBoundary
import quantscale.math.Natural

/**
 * Only works on uniform grid for piecewise linear payoffs. If not uniform, the 5 points around the discontinuity need to be uniformly spaced.
 */
class SplineSmoother(discontinuities: Array[Double]) extends FDPayoffSmoother {
  def makeSmooth(payoff: FDPayoff) {
    val size = payoff.space.length

    val originalState = payoff.state.price;
    val originalSpace = payoff.space;
    val ooSpace = payoff.originalSpace;

    val splineSize = size - discontinuities.length * 3
//    println(splineSize)
    val x = new Array[Double](splineSize)
    val y = new Array[Double](splineSize)
    x(0) = originalSpace(0)
    y(0) = originalState(0)
    var i = 1
    var j = 1
    var k = 0
    while (i < size - 1) {
      while (k < discontinuities.length && originalSpace(i) > (discontinuities(k) + Epsilon.MACHINE_EPSILON_SQRT)) {
        k += 1
      }
      if (k < discontinuities.length && Math.abs(discontinuities(k) - originalSpace(i + 1)) < Epsilon.MACHINE_EPSILON_SQRT) {
        i += 2
      } else {
        x(j) = originalSpace(i)
        y(j) = originalState(i)
        j += 1
      }
      i += 1
    }
    x(j) = originalSpace(i)
    y(j) = originalState(i)
//    println("x=" + Arrays.toString(x))
//    println("y=" + Arrays.toString(y))
    val spline = CubicSpline.makeBesselSpline(x, y)
    i = 0
    while (i < size) {
      originalState(i) = spline.value(originalSpace(i))
      i += 1
    }
//    println("state=" + Arrays.toString(originalState))

  }
}
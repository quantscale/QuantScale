package quantscale.math

import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.scalatest.FunSuite

@RunWith(classOf[JUnitRunner])
class CubicSplineSuite extends FunSuite {
  test("NotAKnotAgainstScilab") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1.0, 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)

    // CubicInterpolation spline =CalypsoCubicSplineFactory.create(x,y);
    val spline = CubicSpline.makeCubicSpline(x, y, new NotAKnot, new NotAKnot)
    val z = Array.ofDim[Double](21)
    val expected = Array(0.0121951, -0.0436779, -0.0510208, -0.0253661,
      0.0177536, 0.0628058, 0.0942579, 0.1001334, 0.2794847, 0.7838355, 1.0,
      0.5539709, 0.2352764, 0.1163079, 0.0688725, 0.0487366,
      0.0347116, 0.0250999, 0.0188466, 0.0148966, 0.0121951)
    // obtained
    // from
    // Scilab.
    z(0) = x(0)
    var i = 1
    while (i < 21) {
      z(i) = z(i - 1) + (18 / 20.0)
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      assert(math.abs(expected(i) - value) < 1e-6)
      i += 1
    }

  }

  test("PCHIMAgainstScilab") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)

    val spline = CubicSpline.makeHarmonicSpline(x, y)

    val z = Array[Double](-9, -8.1, -4, -2.5, -1, 3, 5, 8, 10)

    val expected = Array[Double](0.0121951, 0.0137555, 0.0596466, 0.1375,
      0.6375, 0.1044618, 0.0398548, 0.0147219, 0.0150123)
    // values obtained from Scilab.
    /*
     * deff('[y]=f(t)','y=ones(t)./(1+t.*t)'); t = [-9 -5 -3 -2 0 1 2 4 7
     * 9]; y=f(t); d = splin(t,y,"monotone"); interp([-9 -8.1 -4 -2.5 -1 3 5 8 10],t,y,d, "natural")
     */
    var i = 0
    while (i < z.length) {
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      assert(math.abs(expected(i) - value) < 1e-5)
      i += 1
    }

  }
  test("NaturalExtrapolationAgainstScilab") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)
    val spline = CubicSpline.makeCubicSpline(x, y, new Natural, new Natural)

    val z = Array[Double](-20, -10, -9, -8.1, 9, 10, 20)

    val expected = Array[Double](-1.1106312, 0.0192726, 0.0121951, 0.0056699,
      0.0121951, 0.0090887, -0.3722188)

    // values obtained from Scilab.
    /*
     * deff('[y]=f(t)','y=ones(t)./(1+t.*t)'); t = [-9 -5 -3 -2 0 1 2 4 7
     * 9]; y=f(t); d = splin(t,y,"natural"); interp([-20 -10 -9 -8.1 9 10
     * 20],t,y,d, "natural")
     */
    var i = 0
    while (i < z.length) {
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      assert(math.abs(expected(i) - value) < 1e-5)
      i += 1
    }
  }

  test("NotAKnotExtrapolationAgainstScilab") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)

    val spline = CubicSpline.makeCubicSpline(x, y, new NotAKnot, new NotAKnot)
    val z = Array[Double](-20, -10, -9, -8.1, 9, 10, 20)
    val expected = Array[Double](10.566406, 0.1500863, 0.0121951, -0.0436779,
      0.0121951, 0.0093745, -0.3240691)

    // values obtained from Scilab.
    /*
     * deff('[y]=f(t)','y=ones(t)./(1+t.*t)'); t = [-9 -5 -3 -2 0 1 2 4 7
     * 9]; y=f(t); d = splin(t,y); interp([-20 -10 -9 -8.1 9 10 20],t,y,d,
     * "natural")
     */
    var i = 0
    while (i < z.length) {
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      assert(math.abs(expected(i) - value) < 1e-5)
      i += 1
    }

  }

  test("NaturalAgainstScilab") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)
    val spline = CubicSpline.makeCubicSpline(x, y, new Natural, new Natural)
    val z = Array.ofDim[Double](21)
    val expected = Array[Double](0.0121951, 0.0056699, 0.0031232, 0.0085338, 0.0258802,
      0.0581919, 0.0914586, 0.1006235, 0.2793195, 0.7835630, 1.0,
      0.5539767, 0.2352723, 0.1163160, 0.0688776, 0.0487282,
      0.0346897, 0.0250792, 0.0188565, 0.0149455, 0.0121951)
    // obtained
    // from
    // Scilab.
    z(0) = x(0);
    var maxError = 0.0
    var i = 1
    while (i < 21) {
      z(i) = z(i - 1) + (18 / 20.0)
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      val newError = math.abs(value - expected(i))
      maxError = math.max(maxError, newError)
      assert(math.abs(expected(i) - value) < 1e-6)
      i += 1
    }
    println("maxError=" + maxError)
  }

  test("Hyman83AgainstR") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)
    val spline = CubicSpline.makeCubicSpline(x, y, new Natural, new Natural, MonotoneHyman83DerivativeFilter)
    val z = Array[Double](-9, -8.1, -4, -2.5, -1, 3, 5, 8, 10)
    val expected = Array[Double](0.0121951,
      0.012494290712500002,
      0.0741557,
      0.11250000000000002,
      0.7209582817230875,
      0.09483451181873953,
      0.04032051221768887,
      0.015301547829752635,
      0.00908865217024736)
    var i = 0
    while (i < z.length) {
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      assert(math.abs(expected(i) - value) < 1e-5)
      i += 1
    }
  }

  test("Hyman89AgainstR") {
    val x = Array[Double](-9, -5, -3, -2, 0, 1, 2, 4, 7, 9)
    val y = Array[Double](0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
      0.02, 0.0121951)
    val spline = CubicSpline.makeCubicSpline(x, y, new Natural, new Natural, MonotoneHyman89DerivativeFilter)
    val z = Array[Double](-9, -8.1, -4, -2.5, -1, 3, 5, 8, 10)
    val expected = Array[Double](0.0121951,
      0.012494290712500002,
      0.0741557,
      0.11250000000000002,
      0.7209582817230875,
      0.09483451181873953,
      0.04032051221768887,
      0.015301547829752635,
      0.00908865217024736)
    var i = 0
    while (i < z.length) {
      val value = spline.value(z(i))
      println(z(i) + " " + value + " " + expected(i))
      assert(math.abs(expected(i) - value) < 1e-5)
      i += 1
    }
  }

  test("HymanRPN15A") {
    // test against Fritsch-Carlson RPN 15A data
    val x = Array[Double](7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20)
    val y = Array[Double](0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740,
      0.998636, 0.999919, 0.999994)
    val spline = CubicSpline.makeCubicSpline(x, y, new NotAKnot, new NotAKnot, MonotoneHyman83DerivativeFilter)
    val z = Array.ofDim[Double](1000 + 1)
    var i = 0
    while (i <= 1000) {
      z(i) = (((20 - 7.99) * i) / 1000.0) + 7.99
      i += 1
    }
    println("x y")
    var previous = spline.value(z(0))
    i = 1
    while (i <= 1000) {
      val value = spline.value(z(i))
      println(z(i) + " " + value)
      assert(value >= previous)
      previous = value
      i += 1
    }
  }

}

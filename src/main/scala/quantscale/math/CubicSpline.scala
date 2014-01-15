package quantscale.math

import quantscale.fdm.Epsilon
import quantscale.fdm.TridiagonalMatrix
import java.util.Arrays

sealed trait SplineBoundary

case class NotAKnot extends SplineBoundary

case class FirstDerivative(val value: Double) extends SplineBoundary

case class SecondDerivative(val value: Double) extends SplineBoundary

class Natural extends SecondDerivative(0.0)

trait DerivativeFilter {
  def filter(x: Array[Double], y: Array[Double], yPrime: Array[Double])
}

object NoopDerivativeFilter extends DerivativeFilter {
  def filter(x: Array[Double], y: Array[Double], yPrime: Array[Double]) {}
}

/**
 * The monotonicity is based on James M. Hyman "Accurate Monotonicity Preserving
 * Cubic Interpolation" (1983).
 */
object MonotoneHyman83DerivativeFilter extends DerivativeFilter {
  def filter(x: Array[Double], y: Array[Double], yPrime: Array[Double]) = {
    val n = y.length - 1
    var i = 1
    val dx0 = x(1) - x(0)
    val S0 = (y(1) - y(0)) / dx0
    yPrime(0) = if (S0 >= 0) math.min(math.max(0, yPrime(0)), 3 * math.abs(S0))
    else math.max(math.min(0, yPrime(0)), -3 * math.abs(S0))

    while (i < n) {
      val dxi = x(i + 1) - x(i)
      val Si = (y(i + 1) - y(i)) / dxi
      val dxim = x(i) - x(i - 1)
      val Sim = (y(i) - y(i - 1)) / dxim
      val sig = if (Si * Sim > 0) Si else yPrime(i)
      yPrime(i) = if (sig >= 0) math.min(math.max(0, yPrime(i)), 3 * math.min(math.abs(Si), math.abs(Sim)))
      else math.max(math.min(0, yPrime(i)), -3 * math.min(math.abs(Si), math.abs(Sim)));
      i += 1
    }
    val dxnm1 = x(n) - x(n - 1)
    val Snm1 = (y(n) - y(n - 1)) / dxnm1
    yPrime(n) = if (Snm1 >= 0) math.min(math.max(0, yPrime(n)), 3 * math.abs(Snm1));
    else math.max(math.min(0, yPrime(n)), -3 * math.abs(Snm1));
  }
}

/**
 * monotonic local interpolation scheme based on piecewise cubic
 * polynomials. This follows closely M. Steffen
 * "A simple method for monotonic interpolation in one dimension" published in
 * Astronomy and Astrophysics vol 239, 443-450 (1990).
 */
object MonotoneSteffenDerivativeFilter extends DerivativeFilter {
  def filter(x: Array[Double], y: Array[Double], yPrime: Array[Double]) = {
    throw new UnsupportedOperationException()
  }
}

/**
 * The monotonicity is based on Dougherty, Edelman, Hyman "Nonnegativity,
 * Monotonicity, or Convexity Preserving Cubic and Quintic Hermite
 * Interpolation" (1989).
 */
object MonotoneHyman89DerivativeFilter extends DerivativeFilter {
  def filter(x: Array[Double], y: Array[Double], yPrime: Array[Double]) = {
    //FIXME test this, it seems like it breaks with x nearly uniform
    val n = y.length
    val tmp = yPrime

    var pm = 0.0
    var pu = 0.0
    var pd = 0.0
    var M = 0.0
    val dx0 = x(1) - x(0)
    val S0 = (y(1) - y(0)) / dx0
    var correction = if ((tmp(0) * S0) > 0) math.signum(tmp(0)) * math.min(Math.abs(tmp(0)), math.abs(3.0 * S0)) else 0.0
    if (correction != tmp(0)) {
      tmp(0) = correction
    }
    var i = 1
    while (i < (n - 1)) {
      val dxi = x(i + 1) - x(i)
      val Si = (y(i + 1) - y(i)) / dxi
      val dxim = x(i) - x(i - 1)
      val Sim = (y(i) - y(i - 1)) / dxim
      pm = (Sim * dxi + Si * dxim) / (dxim + dxi)
      M = 3.0 * math.min(math.min(math.abs(Si), math.abs(Sim)), math.abs(pm))
      if (i > 1) {
        val dxim2 = x(i - 1) - x(i - 2)
        val Sim2 = (y(i - 1) - y(i - 2)) / dxim2
        if (((Sim - Sim2) * (Si - Sim)) > 0.0) {
          pd = (Sim * (2.0 * dxim + dxim2) - Sim2 * dxim) / (dxim + dxim2)
          if (((pm * pd) > 0.0) && ((pm * (Sim - Sim2)) > 0.0)) {
            M = math.max(M, 1.5 * math.min(math.abs(pm), math.abs(pd)))
          }
        }
      }
      if (i < (n - 2)) {
        val dxip = x(i + 2) - x(i + 1)
        val Sip = (y(i + 2) - y(i + 1)) / dxip
        if (((Si - Sim) * (Sip - Si)) > 0.0) {
          pu = (Si * (2.0 * dxi + dxip) - Sip * dxi) / (dxi + dxip)
          if (((pm * pu) > 0.0) && ((-pm * (Si - Sim)) > 0.0)) {
            M = math.max(M, 1.5 * math.min(math.abs(pm), math.abs(pu)))
          }
        }
      }
      correction = if ((tmp(i) * pm) > 0.0) math.signum(tmp(i)) * math.min(math.abs(tmp(i)), M) else 0.0

      if (correction != tmp(i)) {
        tmp(i) = correction
      }
      i += 1
    }
    val Snm2 = (y(n - 1) - y(n - 2)) / (x(n - 1) - x(n - 2))
    if ((tmp(n - 1) * Snm2) > 0) {
      correction = math.signum(tmp(n - 1)) * math.min(math.abs(tmp(n - 1)), math.abs(3.0 * Snm2))
    } else {
      correction = 0.0
    }
    if (correction != tmp(n - 1)) {
      tmp(n - 1) = correction
    }

  }
}

object CubicSpline {
  def makeLinearPolynomial(x: Array[Double], y: Array[Double]): CubicPP = {
    val n = 2;
    val pp = new CubicPP(y, new Array[Double](2), new Array[Double](2), new Array[Double](2), x);

    val t = (y(1) - y(0));
    if (Math.abs(x(1) - x(0)) < Epsilon.MACHINE_EPSILON) {
      pp.b(0) = 0.0;
    } else {
      pp.b(0) = t / (x(1) - x(0));
    }
    pp.b(1) = pp.b(0);
    pp.c(0) = 0;
    pp.c(1) = 0;
    pp.d(0) = 0;
    pp.d(1) = 0;
    return pp;
  }

  def makeCubicSpline(x: Array[Double], y: Array[Double], leftBoundary: SplineBoundary = new Natural, rightBoundary: SplineBoundary = new Natural, filter: DerivativeFilter = NoopDerivativeFilter): CubicPP = {
    val yPrime = new Array[Double](y.length)
    computeC2FirstDerivative(x, y, yPrime, leftBoundary, rightBoundary)
    filter.filter(x, y, yPrime)
    return makeHermiteSpline(x, y, yPrime)
  }

  def computeC2FirstDerivative(x: Array[Double], y: Array[Double], yPrime: Array[Double], leftBoundary: SplineBoundary = new Natural, rightBoundary: SplineBoundary = new Natural) {
    val n = y.length;

    if (n == 2) {
      val slope = (y(1) - y(0)) / (x(1) - x(0))
      yPrime(0) = slope
      yPrime(1) = slope
    }
    val a = y;

    val dx = new Array[Double](n - 1);
    val S = new Array[Double](n - 1);

    var i = 0;
    while (i < n - 1) {
      dx(i) = x(i + 1) - x(i);
      if (Math.abs(dx(i)) < Epsilon.MACHINE_EPSILON) {
        throw new IllegalArgumentException("x(" + i + ")=x(" + (i + 1) + ")=" + x(i));
      }
      S(i) = (a(i + 1) - a(i)) / dx(i);
      i += 1;
    }

    val m = new TridiagonalMatrix(n);
    val alpha = new Array[Double](n);
    i = 1;
    while (i < n - 1) {
      m.lower(i) = dx(i);
      m.upper(i) = dx(i - 1);
      m.middle(i) = 2 * (dx(i) + dx(i - 1));
      alpha(i) = 3.0 * (dx(i) * S(i - 1) + dx(i - 1) * S(i));
      i += 1;
    }
    leftBoundary match {
      case NotAKnot() => {
        m.middle(0) = dx(1) * (dx(1) + dx(0));
        m.upper(0) = (dx(0) + dx(1)) * (dx(0) + dx(1));
        alpha(0) = S(0) * dx(1) * (2.0 * dx(1) + 3.0 * dx(0)) + S(1) * dx(0) * dx(0);
      }
      case FirstDerivative(leftValue) => {
        m.middle(0) = 1.0;
        m.upper(0) = 0;
        alpha(0) = leftValue;
      }
      case SecondDerivative(leftValue) => {
        m.middle(0) = 2.0;
        m.upper(0) = 1.0;
        alpha(0) = 3 * S(0) - leftValue * dx(0) / 2.0;

      }
      case _ =>
        throw new IllegalArgumentException("Boundary type " + leftBoundary
          + " is not supported");
    }
    rightBoundary match {
      case NotAKnot() => {
        m.lower(n - 1) = -(dx(n - 2) + dx(n - 3)) * (dx(n - 2) + dx(n - 3));
        m.middle(n - 1) = -dx(n - 3) * (dx(n - 3) + dx(n - 2));
        alpha(n - 1) = -S(n - 3) * dx(n - 2) * dx(n - 2) - S(n - 2) *
          dx(n - 3) * (3.0 * dx(n - 2) + 2.0 * dx(n - 3));
      }

      case FirstDerivative(rightValue) => {
        m.lower(n - 1) = 0;
        m.middle(n - 1) = 1.0;
        alpha(n - 1) = rightValue;
      }
      case SecondDerivative(rightValue) => {
        m.lower(n - 1) = 1.0;
        m.middle(n - 1) = 2.0;
        alpha(n - 1) = 3.0 * S(n - 2) + rightValue * dx(n - 2) / 2.0;
      }
      case _ =>
        throw new IllegalArgumentException("Boundary type " + rightBoundary
          + " is not supported");

    }
    solveTridiagonal(m.lower, m.middle, m.upper, alpha, yPrime, n);
  }

  def makeHermiteSpline(x: Array[Double], y: Array[Double], yPrime: Array[Double]): CubicPP = {
    val n = y.length;

    if (n == 2) {
      return makeLinearPolynomial(x, y)
    }
    val a = y;
    val b = yPrime;
    val c, d = new Array[Double](n)
    c(n - 1) = 0;
    d(n - 1) = 0;
    var i = 0;
    while (i < n - 1) {
      val dx = x(i + 1) - x(i)
      val S = (y(i + 1) - y(i)) / dx
      c(i) = (3 * S - b(i + 1) - 2.0 * b(i)) / dx;
      d(i) = (b(i + 1) + b(i) - 2.0 * S) / (dx * dx);
      i += 1
    }
    return new CubicPP(a, b, c, d, x);
  }

  def computeHarmonicFirstDerivativePCHIM(x: Array[Double], y: Array[Double], yPrime: Array[Double]) = {
    var h0, h1, del0, del1, dmax = 0.0
    val n = x.length;
    h0 = x(1) - x(0);
    del0 = (y(1) - y(0)) / h0;
    h1 = x(2) - x(1);
    del1 = (y(2) - y(1)) / h1;
    var hsum = h0 + h1;
    var w0 = (h0 + hsum) / hsum;
    var w1 = -h0 / hsum;
    val b = yPrime;
    b(0) = w0 * del0 + w1 * del1;

    // SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
    // SHAPE-PRESERVING.
    if (pchst(b(0), del0) <= 0) {
      b(0) = 0;
    } else if (pchst(del0, del1) < 0) {
      // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
      dmax = 3 * del0;
      if (Math.abs(b(0)) > Math.abs(dmax)) {
        b(0) = dmax;
      }
    }
    var i = 1
    while (i < n - 1) {
      if (i != 1) {
        h0 = h1;
        h1 = x(i + 1) - x(i);
        hsum = h0 + h1;
        del0 = del1;
        del1 = (y(i + 1) - y(i)) / h1;
      }
      // SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
      b(i) = 0;
      if (pchst(del0, del1) > 0) {
        // USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
        val hsumt3 = hsum + hsum + hsum;
        w0 = (hsum + h0) / hsumt3;
        w1 = (hsum + h1) / hsumt3;
        dmax = Math.max(Math.abs(del0), Math.abs(del1));
        val dmin = Math.min(Math.abs(del0), Math.abs(del1));
        val drat0 = del0 / dmax;
        val drat1 = del1 / dmax;
        b(i) = dmin / (w0 * drat0 + w1 * drat1);
      }
      i += 1
    }

    // SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
    // SHAPE-PRESERVING.
    w0 = -h1 / hsum;
    w1 = (h1 + hsum) / hsum;
    b(n - 1) = w0 * del0 + w1 * del1;
    if (pchst(b(n - 1), del1) <= 0) {
      b(n - 1) = 0;
    } else if (pchst(del0, del1) <= 0) {
      // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
      dmax = 3 * del1;
      if (Math.abs(b(n - 1)) > Math.abs(dmax)) {
        b(n - 1) = dmax;
      }
    }

  }

  private def pchst(a: Double, b: Double): Double = {
    Math.signum(a) * Math.signum(b)
  }

  //output in yPrime
  def computeParabolicFirstDerivative(x: Array[Double], y: Array[Double], yPrime: Array[Double]) = {
    val n = y.length;
    val dx = new Array[Double](n - 1);
    val S = new Array[Double](n - 1);

    var i = 0;
    while (i < n - 1) {
      dx(i) = x(i + 1) - x(i);
      if (Math.abs(dx(i)) < Epsilon.MACHINE_EPSILON) {
        throw new IllegalArgumentException("x(" + i + ")=x(" + (i + 1) + ")=" + x(i));
      }
      S(i) = (y(i + 1) - y(i)) / dx(i);
      i += 1;
    }

    i = 1
    while (i < n - 1) {
      yPrime(i) = (dx(i - 1) * S(i) + dx(i) * S(i - 1)) / (dx(i) + dx(i - 1));
      i += 1
    }
    // end points natural
    yPrime(0) = ((2.0 * dx(0) + dx(1)) * S(0) - dx(0) * S(1)) /
      (dx(0) + dx(1));
    yPrime(n - 1) = ((2.0 * dx(n - 2) + dx(n - 3)) * S(n - 2) - dx(n - 2) * S(n - 3)) /
      (dx(n - 2) + dx(n - 3));
  }

  def makeBesselSpline(
                        x: Array[Double], y: Array[Double]): CubicPP = {
    val b = new Array[Double](x.length)
    computeParabolicFirstDerivative(x, y, b)
    return makeHermiteSpline(x, y, b)
  }

  def makeHarmonicSpline(x: Array[Double], y: Array[Double]): CubicPP = {
    val b = new Array[Double](x.length)
    computeHarmonicFirstDerivativePCHIM(x, y, b)
    return makeHermiteSpline(x, y, b)
  }

  def makeHymanCubicSpline(x: Array[Double], y: Array[Double]): CubicPP = {
    val b = new Array[Double](x.length)
    computeC2FirstDerivative(x, y, b)
    MonotoneHyman83DerivativeFilter.filter(x, y, b)
    return makeHermiteSpline(x, y, b)
  }

  def makeHymanBesselSpline(x: Array[Double], y: Array[Double]): CubicPP = {
    val b = new Array[Double](x.length)
    computeParabolicFirstDerivative(x, y, b)
    MonotoneHyman83DerivativeFilter.filter(x, y, b)
    return makeHermiteSpline(x, y, b)
  }


  /**
   * Solve M . output = d where M is a tridiagonal matrix. Solving will modify
   * c and d and fill output.
   *
   * @param a
   * , lower diagonal defined from 1 to n-1, a[0] is ignored
   * @param b
   * , middle diagonal defined from 0 to n-1
   * @param c
   * , upper diagonal defined from 0 to n-2, c[n-2] is ignored
   * @param d
   * , right hand side
   * @param output
   * , solution of the tridiagonal system
   * @param n
   * , size on which to solve.
   */
  def solveTridiagonal(a: Array[Double], b: Array[Double], c: Array[Double], d: Array[Double],
                       output: Array[Double], n: Int) {
    // this code has been performance tested and is very likely to perform
    // better than any other alternative.
    c(0) /= b(0);
    d(0) /= b(0);
    var i = 1;
    while (i < n) {
      val l = 1.0 / (b(i) - a(i) * c(i - 1));
      c(i) *= l;
      d(i) = (d(i) - a(i) * d(i - 1)) * l;
      i += 1;
    }
    output(n - 1) = d(n - 1);
    i = n - 2;
    while (i >= 0) {
      output(i) = d(i) - c(i) * output(i + 1);
      i -= 1;
    }
  }

}

/**
 * Cubic Piecewise Polynomial
 */
class CubicPP(

               var a: Array[Double],

               /**
                * array of size n+1
                */
               var b: Array[Double],

               /**
                * array of size n+1
                */
               var c: Array[Double],

               /**
                * array of size n+1
                */
               var d: Array[Double],

               var x: Array[Double]) {

  def value(z: Double): Double = {
    if (z < x(0)) {
      return evalPolynomial(0, z);
    }
    if (z > x(x.length - 1)) {
      return evalPolynomial(x.length - 2, z);
    }
    // no extrapolation to do
    val i = locateLowerIndexBounded(z);
    return evalPolynomial(i, z);
  }

  protected def evalPolynomial(i: Int, z: Double): Double = {
    val h = z - x(i);
    return a(i) + h * (b(i) + h * (c(i) + h * d(i)));
  }

  private def locateLowerIndex(z: Double): Int = {
    if (z < x(0)) {
      return 0;
    }
    if (z > x(x.length - 1)) {
      return x.length - 2;
    }
    return locateLowerIndexBounded(z);
  }

  /**
   * locate lower index of interval containing z given that the precondition
   * is x[0] <= z <= x[n-1] is satisfied
   *
   * @param z
   * @return lower index of interval containing z
   */
  def locateLowerIndexBounded(z: Double): Int = {
    // WARNING the array must be sorted first, we assume that.
    val j = Arrays.binarySearch(x, z);
    return if (j < 0) -j - 2 else j;
  }

  /**
   * Compute the first derivative at a given point
   *
   * @param z
   * @return first derivative value at z
   */
  def derivative(z: Double): Double = {
    var j = locateLowerIndex(z);
    if (j > x.length - 2) {
      j = x.length - 2;
    }
    return evalDerivativePolynomial(j, z);
  }

  def evalDerivativePolynomial(j: Int, z: Double): Double = {
    val h = z - x(j);
    return b(j) + (2.0 * c(j) + 3.0 * d(j) * h) * h;
  }

  /**
   * Compute the second derivative at a given point
   *
   * @param z
   * @return second derivative value at z
   */
  def secondDerivative(z: Double): Double = {
    var j = locateLowerIndex(z);
    if (j > x.length - 2) {
      j = x.length - 2;
    }
    val h = z - x(j);
    return 2.0 * c(j) + 6.0 * h * d(j);
  }
}

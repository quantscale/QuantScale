package quantscale.fdm
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.payoff.AmericanFDPayoff
import java.util.Arrays

class ThomasTridiagonalSolver(var payoff: FDPayoff = null) extends TridiagonalSolver {

  private var c2: Array[Double] = null;
  private var d2: Array[Double] = null;
  var size = 0

  def copy(): TridiagonalSolver = {
    val c = new ThomasTridiagonalSolver(payoff)
    c.init(size)
    return c
  }
  /**
   * initialize internal arrays
   */
  def init(sizeV: Int) {
    size = sizeV
    c2 = new Array[Double](size);
    d2 = new Array[Double](size);
  }

  def lowerBound(): Array[Double] = {
    payoff match {
      case t: AmericanFDPayoff => return payoff.asInstanceOf[AmericanFDPayoff].lowerBound;
      case _                   => return null;
    }
  }

  def solve3(m: TridiagonalMatrix, d: Array[Double], x: Array[Double]): Unit = {
    val a = m.lower;
    val b = m.middle.clone();
    val c = m.upper;
    val v = d.clone();
    val n = c2.size;
    var i = 1;

    while (i < n) {
      var m = a(i) / b(i - 1);
      b(i) = b(i) - m * c(i - 1);
      v(i) = v(i) - m * v(i - 1);
      i += 1;
    }

    x(n - 1) = v(n - 1) / b(n - 1);
    i = n - 2;
    while (i >= 0) {
      x(i) = (v(i) - c(i) * x(i + 1)) / b(i);
      i -= 1;
    }
  }

  def solve2(m: TridiagonalMatrix, d: Array[Double], x: Array[Double]): Unit = {
    val a = m.upper;
    val b = m.middle;
    val c = m.lower.clone();
    val y = d.clone();
    c(0) /= b(0);
    y(0) /= b(0);
    var i = 1;
    val n = c2.size;
    while (i < n) {
      val l = 1 / (b(i) - a(i) * c(i - 1));
      c(i) *= l;
      y(i) = (y(i) - a(i) * y(i - 1)) * l;
      i += 1;
    }
    x(n - 1) = c(n - 1);
    i = n - 2;
    while (i >= 0) {
      x(i) = y(i) - c(i) * x(i + 1);
      i -= 1;
    }
  }
  /**
   * Solve M . x = d where M is a tridiagonal matrix.
   *
   * @param y
   *            , right hand side
   * @param x
   *            , solution of the tridiagonal system
   */

  def solve(m: TridiagonalMatrix, y: Array[Double], x: Array[Double]): Unit = {
    m.reduceBoundaries(y)
    val a = m.lower;
    val b = m.middle;
    val c = m.upper;

    //Solving will modify
    // c2 and d2 and fill output.

    c2(0) = c(0) / b(0);
    d2(0) = y(0) / b(0);
//    if (c2(0).isNaN || d2(0).isNaN) {
//      println(c2(0))
//    }
    var i = 1;
    val n = c2.size;
    while (i < n) {
      val l = 1 / (b(i) - a(i) * c2(i - 1));     
      c2(i) = c(i) * l;
      d2(i) = (y(i) - a(i) * d2(i - 1)) * l;
      i += 1;
    }
    x(n - 1) = d2(n - 1);
    i = n - 2;
    while (i >= 0) {
      x(i) = d2(i) - c2(i) * x(i + 1);
      i -= 1;
    }

    val lowerBoundV = lowerBound()
    if (lowerBoundV != null) {
      i = 0
      while (i < n) {
        x(i) = math.max(x(i), lowerBoundV(i))
        i += 1
      }
    }
  }

}
package quantscale.fdm
import quantscale.fdm.payoff.FDPayoff
import payoff.AmericanFDPayoff
import payoff.AmericanFDPayoff
import java.util.Arrays

class ElliotOckendonTridiagonalSolver(payoff: FDPayoff) extends TridiagonalSolver {

  var BRENNAN_SCHWARTZ: Boolean = true;
  var size = 0;
  private var minDoubleArray: Array[Double] = null

  def copy(): TridiagonalSolver = {
    val c = new ElliotOckendonTridiagonalSolver(payoff)
    c.init(size)
    return c
  }

  def init(sizeV: Int) {
    size = sizeV;
    minDoubleArray = new Array[Double](sizeV)
    Arrays.fill(minDoubleArray, 0.0)
  }

  def lowerBound(): Array[Double] = {
    payoff match {
      case t: AmericanFDPayoff => return payoff.asInstanceOf[AmericanFDPayoff].lowerBound;
      case _                   => return minDoubleArray
    }
  }

  def solve(m: TridiagonalMatrix, d: Array[Double], x: Array[Double]): Unit = {
    m.reduceBoundaries(d)
    val tmpLowerBound = lowerBound();
    val AtimesPayoff = new Array[Double](size);
    val a = m.lower;
    val b = m.middle;
    val c = m.upper;
    m.multiply(tmpLowerBound, AtimesPayoff);
    val v = d;
    val y = new Array[Double](size);
    var i = 0;
    while (i < size) {
      if (math.abs(v(i)) < Epsilon.MACHINE_EPSILON) {
        v(i) = 0d;
      }
      v(i) -= AtimesPayoff(i);
      i += 1
    }
    val isMonotonic = detectSignMonotonicity(v);
    if (!isMonotonic) {
      //            System.out.println("non monotonic input found");
      BRENNAN_SCHWARTZ = true;
    }
    val z = new Array[Double](size);
    // v sign should be monotonic.
    if (v(0) < 0) {
      //          we may simply reverse the order of the components of the vectors q and
      //          z, do likewise to the elements of the main-, sub-, and super-diagonal elements
      //          of M, interchange the sub- and super-diagonals, and then proceed

      val uc = new Array[Double](size);
      i = size - 1;
      var lb = b(i);
      y(i) = v(i) / b(i);
      while ((i > 0) && (BRENNAN_SCHWARTZ || i == 1 || y(i) >= (v(i - 1) / c(i - 1)))) {
        i -= 1;
        uc(i + 1) = a(i + 1) / lb;
        lb = b(i) - c(i) * uc(i + 1);
        y(i) = (v(i) - c(i) * y(i + 1)) / lb;
      }
      val k = i;
      i = 0;
      while (i < k) {
        z(i) = 0;
        i += 1;
      }
      z(k) = math.max(0, y(k));
      i = k + 1;
      while (i < size - 1) {
        z(i) = y(i) - uc(i) * z(i - 1);
        z(i) = math.max(0, z(i));
        i += 1;
      }
    } else {
      //typical for a call
      val uc = new Array[Double](size);
      i = 0;
      var lb = b(0);
      y(0) = v(0) / b(0);
      //the check at i == size-3 is useful for the boundary case
      //we process like Brennan Schwartz for the boundary case.
      //for other cases Elliot Ockendon cuts down the computations
      //compared to Brennan Schwartz
      while ((i < size - 1) && (BRENNAN_SCHWARTZ || i == size - 2 || y(i) >= (v(i + 1) / a(i + 1)))) {
        i += 1;
        uc(i - 1) = c(i - 1) / lb;
        lb = b(i) - a(i) * uc(i - 1);
        y(i) = (v(i) - a(i) * y(i - 1)) / lb;
      }
      val k = i;
      i = size - 1;
      while (i > k) {
        z(i) = 0;
        i -= 1;
      }
      z(k) = math.max(0, y(k));
      i = k - 1;
      while (i >= 0) {
        z(i) = y(i) - uc(i) * z(i + 1);
        z(i) = math.max(0, z(i));
        i -= 1;
      }
    }

    i = 0;
    while (i < size) {
      x(i) = z(i) + tmpLowerBound(i);
      i += 1;
    }
  }

  private def detectSignMonotonicity(v: Array[Double]): Boolean = {
    var signChange: Int = 0;
    var prev = (v(0) >= 0);
    var i = 1;
    while (signChange <= 1 && i < v.length - 1) {
      val isPositive = (v(i) >= 0);
      if (isPositive != prev) {
        signChange += 1;
        prev = isPositive;
      }
      i += 1;
    }
    return signChange <= 1;
  }

}
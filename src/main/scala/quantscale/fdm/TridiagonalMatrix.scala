package quantscale.fdm
import java.util.Arrays

class DifferentialCache(val x: Array[Double]) {
  var h: Array[Double] = new Array[Double](x.length)
  var hm: Array[Double] = new Array[Double](x.length)
  var hOverHm: Array[Double] = new Array[Double](x.length)
  var denom: Array[Double] = new Array[Double](x.length)
  for (i <- 1 until x.length - 1) {
    h(i) = x(i + 1) - x(i)
    hm(i) = (x(i) - x(i - 1))
    hOverHm(i) = h(i) / hm(i)
    denom(i) = 1.0 / (h(i) * hm(i) * (1 + hOverHm(i)))
  }
}

class OperatorLine(val iStart: Int, val iEnd: Int) {
  private val v = new Array[Double](iEnd - iStart)

  def fill(z: Double) {
    Arrays.fill(v, z)
  }

  def value(i: Int) = v(i - iStart)
  
  def update(i:Int, value: Double) = v(i-iStart)=value

  def plusD1Forward(i: Int, x: Array[Double], z: Double): OperatorLine = {
    val hi = x(i + 1) - x(i)
    v(i - iStart) += -z / hi
    v(i - iStart + 1) += z / hi
    return this
  }

  def plusD1Backward(i: Int, x: Array[Double], z: Double): OperatorLine = {
    val hi = x(i) - x(i - 1)
    v(i - iStart) += z / hi
    v(i - iStart - 1) += -z / hi
    return this
  }

  def plus(i: Int, z: Double): OperatorLine = {
    v(i - iStart) += z
    return this
  }

  def plusD2Backward(i: Int, x: Array[Double], z: Double): OperatorLine = {
    val hi = x(i) - x(i - 1)
    val him = x(i - 1) - x(i - 2)
    val hiOverHim = hi / him
    val denom = 1.0 / (hi * him * (1 + hiOverHim))
    v(i - iStart - 2) += 2 * hiOverHim * denom * z
    v(i - iStart - 1) += -2 * (1 + hiOverHim) * denom * z
    v(i - iStart) += 2 * denom * z
    return this
  }

  //TODO test this
  def plusD2Forward(i: Int, x: Array[Double], z: Double): OperatorLine = {
    val hi = x(i + 2) - x(i + 1)
    val him = x(i + 1) - x(i)
    val hiOverHim = hi / him
    val denom = 1.0 / (hi * him * (1 + hiOverHim))
    v(i - iStart + 2) += 2 * hiOverHim * denom * z
    v(i - iStart + 1) += -2 * (1 + hiOverHim) * denom * z
    v(i - iStart) += 2 * denom * z
    return this
  }

  //TODO test this
  def plusD1Forward2(i: Int, x: Array[Double], z: Double): OperatorLine = {
    val hi = x(i + 1) - x(i)
    val him = x(i + 2) - x(i + 1)
    val hiOverHim = hi / him
    //    val denom = 1.0 / (hi * him * (1 + hiOverHim))
    val alpha3 = -hiOverHim / (him * (1 + hiOverHim))
    val alpha2 = +(1 / him + 1 / hi)
    val alpha1 = -alpha2 - alpha3
    v(i - iStart + 2) += alpha3 * z
    v(i - iStart + 1) += alpha2 * z
    v(i - iStart) += alpha1 * z
    return this
  }

  def plusD1Backward2(i: Int, x: Array[Double], z: Double): OperatorLine = {
    val hi = x(i) - x(i - 1)
    val him = x(i - 1) - x(i - 2)
    val hiOverHim = hi / him
    //    val denom = 1.0 / (hi * him * (1 + hiOverHim))
    val alpha3 = hiOverHim / (him * (1 + hiOverHim))
    val alpha2 = -(1 / him + 1 / hi)
    val alpha1 = -alpha2 - alpha3
    v(i - iStart - 2) += alpha3 * z
    v(i - iStart - 1) += alpha2 * z
    v(i - iStart) += alpha1 * z
    return this
  }

}

object TridiagonalMatrix {
  def identity(size: Int): TridiagonalMatrix = {
    val m = new TridiagonalMatrix(size)
    Arrays.fill(m.middle, 1.0)
    return m
  }

}

class TridiagonalMatrix(
  val size: Int,
  val lower: Array[Double], //public getter, no setter
  val middle: Array[Double],
  val upper: Array[Double],
  var firstLine: OperatorLine = null,
  var lastLine: OperatorLine = null) {

  def this(sizeV: Int) = this(sizeV, new Array[Double](sizeV), new Array[Double](sizeV), new Array[Double](sizeV))

  def copy(): TridiagonalMatrix = {
    return new TridiagonalMatrix(size, lower.clone, middle.clone, upper.clone, firstLine, lastLine)
  }

  /**
   * Modifies this matrix to add z*M
   * firstLine and lastLine of M are taken into account
   */
  def plus(M: TridiagonalMatrix, z: Double): TridiagonalMatrix = {
    var i = size - 2
    while (i >= 1) {
      lower(i) += z * M.lower(i)
      middle(i) += z * M.middle(i)
      upper(i) += z * M.upper(i)
      i -= 1
    }
    if (M.firstLine == null) {
      lower(0) += z * M.lower(0)
      middle(0) += z * M.middle(0)
      upper(0) += z * M.upper(0)
    } else {
      firstLine = new OperatorLine(M.firstLine.iStart, M.firstLine.iEnd)
      i = M.firstLine.iStart
      firstLine(i) = middle(i) + z * M.firstLine.value(i)
      i += 1
      firstLine(i) = upper(i-1) + z * M.firstLine.value(i)
      i += 1
      while (i < firstLine.iEnd) {
        firstLine(i) = z * M.firstLine.value(i)
        i += 1
      }
    }
    if (M.lastLine == null) {
      lower(size - 1) += z * M.lower(size - 1)
      middle(size - 1) += z * M.middle(size - 1)
      upper(size - 1) += z * M.upper(size - 1)
    } else {
      lastLine = new OperatorLine(M.lastLine.iStart, M.lastLine.iEnd)
      i = M.lastLine.iEnd-1
      lastLine(i) = middle(i) + z * M.lastLine.value(i)
      i -= 1
      lastLine(i) = lower(i+1) + z * M.lastLine.value(i)
      i -= 1
      while (i >= lastLine.iStart) {
        lastLine(i) = z * M.lastLine.value(i)
        i -= 1
      }
    }
    return this
  }

  def multiply(x: Array[Double], y: Array[Double]) {
    var i = 1;
    while (i < size - 1) {
      y(i) = lower(i) * x(i - 1) + middle(i) * x(i) + upper(i) * x(i + 1);
      i += 1;
    }
    if (firstLine == null) {
      y(0) = upper(0) * x(1) + middle(0) * x(0);
    } else {
      i = firstLine.iStart
      var sum = 0.0
      while (i < firstLine.iEnd) {
        sum += firstLine.value(i) * x(i)
        i += 1
      }
      y(0) = sum
    }
    if (lastLine == null) {
      y(size - 1) = lower(size - 1) * x(size - 2) + middle(size - 1) * x(size - 1);
    } else {
      i = lastLine.iStart
      var sum = 0.0
      while (i < lastLine.iEnd) {
        sum += lastLine.value(i) * x(i)
        i += 1
      }
      y(size - 1) = sum
    }
  }

  override def toString(): String = {
    return "lower=" + Arrays.toString(lower) + "\n" + "middle=" + Arrays.toString(middle) + "\n" + "upper=" + Arrays.toString(upper)
  }

  def fill(z: Double) {
    Arrays.fill(lower, z)
    Arrays.fill(middle, z)
    Arrays.fill(upper, z)
  }

  def parabolicOperator(iStart: Int, iEnd: Int, x: DifferentialCache, m2: Array[Double], d2: Double, m1: Array[Double], d1: Double, d0: Double) {
    for (i <- iStart until iEnd) {
      val him = x.hm(i)
      val hiOverHim = x.hOverHm(i)
      val denom = x.denom(i)
      lower(i) = 2 * hiOverHim * denom * m2(i) * d2 -
        hiOverHim * hiOverHim * denom * him * m1(i) * d1
      middle(i) = -2 * (1 + hiOverHim) * denom * m2(i) * d2 -
        (1 - hiOverHim * hiOverHim) * denom * him * m1(i) * d1 +
        d0
      upper(i) = 2 * denom * m2(i) * d2 +
        him * m1(i) * denom * d1
    }
  }
  /**
   * second derivative operator
   */
  def plusD2(iStart: Int, iEnd: Int, x: DifferentialCache, m: Array[Double], d: Double): TridiagonalMatrix = {
    for (i <- iStart until iEnd) {
      val hiOverHim = x.hOverHm(i)
      val denom = x.denom(i)
      lower(i) += 2 * hiOverHim * denom * m(i) * d
      middle(i) += -2 * (1 + hiOverHim) * denom * m(i) * d
      upper(i) += 2 * denom * m(i) * d
    }
    return this
  }

  /**
   * first derivative operator
   */
  def plusD1Central(iStart: Int, iEnd: Int, x: DifferentialCache, m: Array[Double], d: Double): TridiagonalMatrix = {
    for (i <- iStart until iEnd) {
      val hi = 1.0 / (x.x(i + 1) - x.x(i - 1))
      lower(i) += -hi * m(i) * d
      upper(i) += hi * m(i) * d
    }
    return this
  }
  /**
   * first derivative operator
   */
  def plusD1(iStart: Int, iEnd: Int, x: DifferentialCache, m: Array[Double], d: Double): TridiagonalMatrix = {
    for (i <- iStart until iEnd) {
      val hiOverHimSq = x.hOverHm(i) * x.hOverHm(i)
      val denom = x.denom(i) * x.hm(i)
      lower(i) += -hiOverHimSq * denom * m(i) * d
      middle(i) += -(1 - hiOverHimSq) * denom * m(i) * d
      upper(i) += m(i) * denom * d
    }
    return this
  }

  def plusD0(iStart: Int, iEnd: Int, m: Array[Double], d: Double): TridiagonalMatrix = {
    for (i <- iStart until iEnd) {
      middle(i) += d * m(i)
    }
    return this
  }

  def plusD0(iStart: Int, iEnd: Int, m: Array[Double], d: Double, z: Double): TridiagonalMatrix = {
    for (i <- iStart until iEnd) {
      middle(i) += m(i) * d + z
    }
    return this
  }

  def setBoundaries(firstLine: OperatorLine, lastLine: OperatorLine) {
    this.firstLine = firstLine
    this.lastLine = lastLine
  }

  def reduceBoundaries(y: Array[Double]) {
    val a = lower;
    val b = middle;
    val c = upper;
    if (firstLine != null) {
      //      b(0) = 1
      if (firstLine.iEnd == 3 && math.abs(firstLine.value(2)) > Epsilon.MACHINE_EPSILON_SQRT) {
        val factor = c(1) / firstLine.value(2)
        b(0) = a(1) - factor * firstLine.value(0)
        c(0) = b(1) - factor * firstLine.value(1)
        y(0) = y(1) - factor * y(0)
      } else {
        b(0) = firstLine.value(0)
        c(0) = firstLine.value(1)
      }
    }
    if (lastLine != null) {
      val n = b.length
      //      b(n-1) = 1
      if (lastLine.iStart == n - 3 && math.abs(lastLine.value(n - 3)) > Epsilon.MACHINE_EPSILON_SQRT) {
        val factor = a(n - 2) / lastLine.value(n - 3)
        a(n - 1) = b(n - 2) - factor * lastLine.value(n - 2)
        b(n - 1) = c(n - 2) - factor * lastLine.value(n - 1)
        y(n - 1) = (y(n - 2) - factor * y(n - 1))
      } else {
        a(n - 1) = lastLine.value(n - 2)
        b(n - 1) = lastLine.value(n - 1)
      }
    }
  }
}
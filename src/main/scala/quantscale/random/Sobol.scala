package quantscale.random

import java.util.Scanner
import java.io.FileOutputStream
import java.io.File
import java.io.ObjectOutputStream
import java.io.ObjectInputStream

object SobolFactory extends RandomSequenceFactory {
  def makeRandomSequenceGenerator(dimension: Int): RandomSequenceGenerator = {
    return new Sobol(dimension, 0, TextFileSobolDirectionNumbers)
  }
}

trait SobolInitializationData {
  def directionNumbers(dimension: Int): Array[Int]
}

object TextFileSobolDirectionNumbers extends SobolInitializationData {
  val defaultFileName = "new-joe-kuo-6.21201"
  val defaultBinFileName = "new-joe-kuo-6.21201.bin"
  private var allDirectionNumbers = new Array[Array[Int]](21202)

  readAll(defaultFileName)
  //  readAllBinary(defaultBinFileName)

  def readAllBinary(fileName: String) {
    val r = this.getClass().getResourceAsStream(fileName)
    val in = new ObjectInputStream(r)
    allDirectionNumbers = in.readObject().asInstanceOf[Array[Array[Int]]]
    in.close()
  }

  def readAll(fileName: String) {
    val url = this.getClass().getResourceAsStream(fileName)
    // Turn the resource into a File object
    //    val f = new File(url.toURI())
    var scan = new Scanner(url)
    var line = scan.nextLine()
    while (scan.hasNextLine()) {
      line = scan.nextLine()
      val (d, aAndm) = lineToDirectionNumbers(line)
      allDirectionNumbers(d) = aAndm
    }
  }

  def writeBinaryFile(fileName: String) {
    val f = new File(fileName)
    val out = new ObjectOutputStream(new FileOutputStream(f))
    out.writeObject(allDirectionNumbers)
    out.close()
  }

  private def lineToDirectionNumbers(line: String): (Int, Array[Int]) = {
    //    var scanner = new Scanner(line);
    //    var d = scanner.nextInt()
    //    var s = scanner.nextInt()
    //    var a = scanner.nextInt()
    //    var m = Array.ofDim[Int](s)
    //    var i = 0
    //    while (i < s) {
    //      m(i) = scanner.nextInt()
    //      i += 1
    //    }
    var strs = line.trim().split("\\s+")
    var d = strs(0).toInt
    var s = strs(1).toInt
    var a = strs(2).toInt
    var m = Array.ofDim[Int](s + 1)
    var i = 0
    m(i) = a
    while (i < s) {
      m(i + 1) = strs(3 + i).toInt
      i += 1
    }
    return (d, m)
  }

  def directionNumbers(dimension: Int): Array[Int] = {
    allDirectionNumbers(dimension + 1)
  }

  def main(args: Array[String]) {
    writeBinaryFile(defaultFileName + ".bin")
  }
}

//class CachedSobolDirectionNumbers(delegate: SobolDirectionNumbers) extends SobolDirectionNumbers {
//  private var _cache = Map[Int, Array[Int]]()
//
//  def numbersByDimension(dimension: Int): Array[Int] = {
//    var o = _cache.get(dimension)
//    if (o.isEmpty) {
//      var m = delegate.numbersByDimension(dimension)
//      _cache += (dimension -> m)
//      return m
//    }
//    return o.get
//  }
//}

/**
 * Sobol quasi random generator that can work with Joe-Kuo direction numbers for up to 21000 dimensions.
 *
 */
class Sobol(D: Int, N: Int, directions: SobolInitializationData) extends RandomSequenceGenerator {
  private val NORM = 2.32830643653869628906e-10;
  private val MASK_UNSIGNED_INT = 0xffffffffL;

  private var L, C, counter = 0
  private var X: Array[Int] = null
  private var V: Array[Array[Int]] = null
  private var skipLen = 0
  private val LOG2 = Math.log(2.0);
  private val WORDSIZE = 32;

  init()

  private def init() {
    if (N == 0) {
      this.L = WORDSIZE;
    } else {
      this.L = math.ceil(Math.log(N) / LOG2).toInt;
    }
    this.counter = 0;
    this.C = 1;
    this.X = Array.ofDim[Int](D)
    var i = 0
    while (i < X.length) {
      X(i) = 0;
      i += 1
    }
    V = Array.ofDim[Int](L + 1, D);
    computeV()
  }

  private def computeV() {
    var i = 1
    while (i <= L) {
      V(i)(0) = 1 << (WORDSIZE - i);
      i += 1
    }
    var j = 1
    while (j <= D - 1) {

      val dn = directions.directionNumbers(j)
      var a = dn(0)
      var s = dn.length - 1

      if (L <= s) {
        i = 1
        while (i <= L) {
          V(i)(j) = dn(i) << (WORDSIZE - i);
          i += 1
        }

      } else {
        i = 1
        while (i <= s) {
          V(i)(j) = dn(i) << (WORDSIZE - i);
          i += 1
        }
        i = s + 1
        while (i <= L) {
          V(i)(j) = V(i - s)(j) ^ (V(i - s)(j) >>> s);
          var k = 1
          while (k <= s - 1) {
            V(i)(j) ^= (((a >>> (s - 1 - k)) & 1) * V(i - k)(j));
            k += 1
          }
          i += 1
        }
      }
      j += 1
    }
  }

  /**
   * skip n-th sample in the low-discrepancy sequence e.g. skipTo(10)
   * is equivalent to calling 10 times next().
   *
   * This supports only skip of length < 2^31
   */
  def skip(skip: Int) {

    var N = (skip + counter) << 1;

    // Convert to Gray code
    var G = N ^ (N >>> 1);
    var j = 0
    while (j < D) {
      X(j) = 0;
      var index = 1
      while (index < 32) { //from 1 to L inclusive
        if (((G >>> index) & 1) != 0) {
          X(j) ^= V(index)(j);
        }
        index += 1
      }
      j += 1
    }
    counter = skip;
    C = ffzJava(counter);
  }

  def nextIntSequence(): Array[Int] = {
    counter += 1;

    val V_C = V(C);
    var j = 0
    while (j < D) {
      X(j) ^= V_C(j);
      j += 1
    }

    var value = counter;
    C = ffzJava(value);
    return X;
  }

  def nextSequence(POINTS: Array[Double]) = {
    counter += 1;
    var V_C = V(C);
    var j = 0
    while (j < D) {
      X(j) ^= V_C(j)
      POINTS(j) = (X(j) & MASK_UNSIGNED_INT) * NORM
      j += 1
    }
    var value = counter;
    C = ffzJava(value);
  }

  def ffzSlow(v: Int): Int = {
    var i = 1;
    var value = v
    while ((value & 1) != 0) {
      value >>>= 1;
      i += 1;
    }
    return i;
  }

  def ffzJava(value: Int): Int = {
    return Integer.numberOfTrailingZeros(Integer.lowestOneBit(~value)) + 1;
  }

  def ffzFast(value: Int): Int = {
    // as sometimes done at the OS level
    var x = ~value;
    var r = 1;

    if (x == 0)
      return 0;
    if ((x & 0xffff) == 0) {
      x >>= 16;
      r += 16;
    }
    if ((x & 0xff) == 0) {
      x >>= 8;
      r += 8;
    }
    if ((x & 0xf) == 0) {
      x >>= 4;
      r += 4;
    }
    if ((x & 3) == 0) {
      x >>= 2;
      r += 2;
    }
    if ((x & 1) == 0) {
      x >>= 1;
      r += 1;
    }
    return r;

  }

  def dimension(): Int = {
    return this.D;
  }
}
package quantscale.random

object Scrambling extends Enumeration {
  type Scrambling = Value
  val NONE, OWEN, FAURE_TEZUKA, OWEN_FAURE_TEZUKA = Value
}
import Scrambling._


object ScrambledSobolFactory extends RandomSequenceFactory {
  def makeRandomSequenceGenerator(dimension: Int): RandomSequenceGenerator = {
    return new ScrambledSobol(dimension, 0, TextFileSobolDirectionNumbers, OWEN_FAURE_TEZUKA)
  }
}
/**
 * A Sobol generator with scrambling.
 *
 * OWEN, FAURE_TEZUKA, OWEN_FAURE_TEZUKA scrambling from Algorithm 823 COLLECTED ALGORITHMS FROM ACM. TRANSACTIONS ON
 * MATHEMATICAL SOFTWARE, VOL. 29, NO. 2, June, 2003, P. 95--109.
 *
 * NONE scrambling will produce the same numbers as Sobol
 *
 */
class ScrambledSobol(D: Int, N: Int, directions: SobolInitializationData, scrambling: Scrambling) extends RandomSequenceGenerator {
  private val NORM = 2.32830643653869628906e-10;
  private val MASK_UNSIGNED_INT = 0xffffffffL;

  private var L, C, counter = 0
  private var X: Array[Int] = null
  private var V: Array[Array[Int]] = null
  private var skipLen = 0
  private val LOG2 = Math.log(2.0)
  private val WORDSIZE = 32
  private var shift: Array[Int] = null
  private var ssobol_seedi = 0
  private var ssobol_seedj = 0
  private var ssobol_seedcarry = 0.0
  private var ssobol_seedseeds = Array.ofDim[Double](24)
  private var maxd = WORDSIZE - 2
  private var norm: Double = NORM;

  init()

  private def init() {
    if (N == 0) {
      this.L = WORDSIZE - 2;
    } else {
      this.L = math.ceil(Math.log(N) / LOG2).toInt;
    }
    this.counter = 0;
    this.C = 1;

    V = Array.ofDim[Int](L + 1, D)

    computeV();
    scramble();
    this.X = Array.ofDim[Int](D)
    var i = 0
    while (i < X.length) {
      X(i) = 0;
      i += 1
    }
  }

  private def scramble() {
    shift = Array.ofDim(D);
    if (scrambling == Scrambling.NONE) {
      norm = math.pow(2, -L);
    } else {
      val scrambledV = Array.ofDim[Int](L + 1, D);
      maxd = 30;
      if (scrambling == Scrambling.OWEN || scrambling == Scrambling.OWEN_FAURE_TEZUKA) {
        norm = Math.pow(2, -maxd);

        var lsm = genscrml(shift);
        var i = 1
        while (i <= D) {
          var j = 1
          while (j <= L) {
            var l = 1;
            var temp2 = 0;
            var p = maxd;
            while (p >= 1) { // or L?
              var temp1 = 0;
              var k = 1
              while (k <= L) {
                temp1 += lbitbits(lsm(p)(i - 1), k - 1, 1) * lbitbits(V(j)(i - 1), k - 1, 1);
                k -= 1
              }
              temp1 %= 2;
              temp2 += temp1 * l;
              l <<= 1;
              p -= 1
            }
            scrambledV(j)(i - 1) = temp2;
            j += 1
          }
          i += 1
        }
      }
      if (scrambling == Scrambling.FAURE_TEZUKA || scrambling == Scrambling.OWEN_FAURE_TEZUKA) {
        val ushift = Array.ofDim[Int](L)
        val usm = genscrmu(ushift);
        val tv = Array.ofDim[Int](D, L, L)
        var maxx = if (scrambling == Scrambling.FAURE_TEZUKA) L else maxd
        var i = 1
        while (i <= D) {
          var j = 1
          while (j <= L) {
            var p = maxx;
            var k = 1
            while (k <= maxx) {
              if (scrambling == Scrambling.FAURE_TEZUKA) {
                tv(i - 1)(p - 1)(j - 1) = lbitbits(V(j)(i - 1), k - 1, 1);
              } else {
                tv(i - 1)(p - 1)(j - 1) = lbitbits(scrambledV(j)(i - 1), k - 1, 1);
              }
              p -= 1;
              k += 1
            }
            j += 1
          }
          var pp = 1
          while (pp <= L) {
            var temp2 = 0;
            var temp4 = 0;
            var l = 1;
            var j = maxx
            while (j >= 1) {
              var temp1 = 0;
              var temp3 = 0;
              var p = 1
              while (p <= L) {
                temp1 += tv(i - 1)(j - 1)(p - 1) * usm(p - 1)(pp - 1);
                if (pp == 1) {
                  temp3 += tv(i - 1)(j - 1)(p - 1) * ushift(p - 1);
                }
                p += 1
              }
              temp1 %= 2;
              temp2 += temp1 * l;
              if (pp == 1) {
                temp3 %= 2;
                temp4 += temp3 * l;
              }
              l <<= 1;
              j -= 1
            }
            scrambledV(pp)(i - 1) = temp2;
            if (pp == 1) {
              if (scrambling == Scrambling.OWEN_FAURE_TEZUKA) {
                shift(i - 1) = temp4 ^ shift(i - 1); // exor
              } else {
                shift(i - 1) = temp4;
              }
            }
            pp += 1
          }
          i += 1
        }
        norm = math.pow(2.0, -maxx)

      }
      V = scrambledV;
    }
  }

  private def genscrml(shift: Array[Int]): Array[Array[Int]] = {
    var lsm = Array.ofDim[Int](maxd + 1, D);
    var p = 1
    while (p <= D) {
      shift(p - 1) = 0;
      var l = 1;
      var i = maxd
      while (i >= 1) {
        lsm(i)(p - 1) = 0;
        var stemp = (unirnd() * 1e3f).toInt % 2;
        shift(p - 1) += stemp * l;
        l <<= 1;
        var ll = 1;
        var j = L;
        while (j >= 1) {
          var temp = 0
          if (j == i) {
            temp = 1;
          } else if (j < i) {
            temp = (unirnd() * 1e3f).toInt % 2;
          } else {
            temp = 0;
          }
          lsm(i)(p - 1) += temp * ll;
          ll <<= 1;
          j -= 1
        }
        i -= 1
      }
      p += 1
    }
    return lsm;
  }

  private def genscrmu(ushift: Array[Int]): Array[Array[Int]] = {
    var usm = Array.ofDim[Int](L, L)
    var i = 1
    var temp = 0
    while (i <= L) {
      var stemp = (unirnd() * 1e3f).toInt % 2;
      ushift(i - 1) = stemp;
      var j = 1
      while (j <= L) {
        if (j == i) {
          temp = 1;
        } else if (j > i) {
          temp = (unirnd() * 1e3f).toInt % 2;
        } else {
          temp = 0;
        }
        usm(i - 1)(j - 1) = temp;
        j += 1
      }
      i += 1
    }
    return usm;
  }

  def resetSeed() {
    val seeds = Array(.8804418, .2694365, .0367681, .4068699, .4554052, .2880635, .1463408, .2390333, .6407298,
      .1755283, .713294, .4913043, .2979918, .1396858, .3589528, .5254809, .9857749, .4612127, .2196441,
      .7848351, .40961, .9807353, .2689915, .5140357)
    setSeed(seeds);
  }

  def setSeed(seeds: Array[Double]) {
    var i = 0;
    ssobol_seedi = 24;
    ssobol_seedj = 10;
    ssobol_seedcarry = 0.;
    while (i < 24) {
      ssobol_seedseeds(i) = seeds(i)
      i += 1
    }
  }

  def unirnd(): Double = {
    var ret_val = ssobol_seedseeds(ssobol_seedi - 1) - ssobol_seedseeds(ssobol_seedj - 1) - ssobol_seedcarry;
    if (ret_val < 0.) {
      ret_val += 1;
      ssobol_seedcarry = 5.9604644775390625e-8;
    } else {
      ssobol_seedcarry = 0.;
    }
    ssobol_seedseeds(ssobol_seedi - 1) = ret_val;
    ssobol_seedi = 24 - (25 - ssobol_seedi) % 24;
    ssobol_seedj = 24 - (25 - ssobol_seedj) % 24;
    return ret_val;
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

  def skip(skip: Int) {

    var N = (skip + counter) << 1;

    // Convert to Gray code
    var G = N ^ (N >>> 1);
    var j = 0
    while (j < D) {
      X(j) = 0;
      var index = 1
      while (index < 32) { //could restrict from 1 to L only
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

  /**
   * generate quasi random numbers in (0,1), 0 and 1 excluded.
   */
  def nextSequence(POINTS: Array[Double]) = {
    counter += 1;
    var V_C = V(C);
    var j = 0
    while (j < D) {
      X(j) ^= V_C(j)
      val X_j = X(j)
      POINTS(j) = if (X_j == 0) norm / 2 else (X_j & MASK_UNSIGNED_INT) * norm;
      j += 1
    }
    var value = counter;
    C = ffzJava(value);
  }

  def lbitbits(a: Int, b: Int, len: Int): Int = {
    /* Assume 2's complement arithmetic */

    var x: Long = a;
    var y: Long = -1L;
    x >>>= b;
    y <<= len;
    return (x & ~y).toInt;
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
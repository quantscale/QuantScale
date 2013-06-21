package quantscale.random

class Well1024a extends RandomGeneratorDouble {
  //2^-32
  private val NORM = 2.32830643653869628906e-10;
  private val MASK_UNSIGNED_INT = 0xffffffffL;
  private val MASK: Int = 0x1F;
  private val W: Int = 32;
  private val R = 32;
  private val P = 0;
  private val M1 = 3;
  private val M2 = 24;
  private val M3 = 10;

  private val A1 = 0xDB79FB31;
  private val B1 = 0x0FDAA2C1;
  private val B2 = 0x27C5A58D;
  private val C1 = 0x71E993FE;

  private var state_i = 0
  private var STATE = Array.ofDim[Int](R)

  init()

  private def init() {
    state_i = 0
    STATE(0) = 19650218
    var j = 1
    while (j < R) {
      STATE(j) = (1812433253 * (STATE(j - 1) ^ (STATE(j - 1) >>> 30)) + j)
      j += 1
    }
  }

  final def nextInt32(): Int = {
    // this code is an adaptation of the one in the paper "Improved
    // Long-Period
    // Generators Based on Linear Recurrences Modulo 2"
    val z0 = STATE((state_i + 31) & MASK);
    val z1 = STATE(state_i) ^ MAT3POS(8, STATE((state_i + M1) & MASK));
    val z2 = MAT3NEG(-19, STATE((state_i + M2) & MASK)) ^ MAT3NEG(-14, STATE((state_i + M3) & MASK));
    STATE(state_i) = z1 ^ z2;
    STATE((state_i + 31) & MASK) = MAT3NEG(-11, z0) ^ MAT3NEG(-7, z1) ^ MAT3NEG(-13, z2);
    state_i = (state_i + 31) & MASK;
    return STATE(state_i);
  }

  def nextDouble(): Double = {
    return (nextInt32() & MASK_UNSIGNED_INT) * NORM
  }
  def nextDoubleOpen(): Double = {
    ((nextInt32() & MASK_UNSIGNED_INT) + 0.5) * NORM
  }

  private def MAT3POS(t: Int, v: Int): Int = {
    return (v ^ (v >>> t));
  }

  private def MAT3NEG(t: Int, v: Int): Int = {
    return (v ^ (v << (-t)));
  }

}
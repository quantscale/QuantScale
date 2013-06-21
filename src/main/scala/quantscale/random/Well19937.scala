package quantscale.random

class Well19937(tempering: Boolean = true) extends RandomGeneratorDouble {
  private val NORM = 2.32830643653869628906e-10;
  private val MASK_UNSIGNED_INT = 0xffffffffL;

  private val CP_DEGREE = 19937;

  // k = 1024, w = 32, r = 32, p = 0
  //
  // M1 M3(8) M3(-19) M3(-14) 0
  // 3 24 10 M3(-11) M3(-7) M3(-13) M0 407

  private val W = 32;

  private val R = 624;
  private val P = 31;
  private val MASKU = (0xffffffffL >>> (W - P)).toInt;
  private val MASKL = (~(0xffffffffL >>> (W - P))).toInt;

  private val M1 = 70;
  private val M2 = 179;
  private val M3 = 449;

  private val TEMPERB = 0xe46e1700;
  private val TEMPERC = 0x9b868000;

  //    // temporary cache for sliding window jumps
  //   private val                   QQ        = 6;
  //    // =2^QQ
  //   private val                   LL        = 1 << QQ;
  //
  //    private var                            vec_h     = Array.ofDim[Int](LL,R + 1)
  //
  //    // gray codes table for the cache
  //    private var                              h         = Array.ofDim[Int](LL)

  private def MAT0POS(t: Int, v: Int): Int = {
    return (v ^ (v >>> t));
  }

  private def MAT0NEG(t: Int, v: Int): Int = {
    return (v ^ (v << (-t)));
  }

  private def MAT3POS(t: Int, v: Int): Int = {
    return (v >>> t);
  }

  private var z0: Int = 0
  private var z1: Int = 0
  private var z2: Int = 0
  private var state_i: Int = 0
  private var STATE = Array.ofDim[Int](R)
  private var caseNumber = 1;

  init()

  /**
   *
   * @param tempering
   *            , if true, this is equivalent to Well19937c and is better
   *            equidistributed
   */
  def init() {
    state_i = 0;

    STATE(0) = 19650218;
    var j = 1
    while (j < R) {
      STATE(j) = (1812433253 * (STATE(j - 1) ^ (STATE(j - 1) >>> 30)) + j);
      j += 1
    }
  }

  /**
   *
   * @param init
   *            array of size 32 of used as initial seed.
   */
  def this(init: Array[Int], tempering: Boolean) {
    this(tempering);
    var j = 0;
    state_i = 0;

    while (j < R) {
      STATE(j) = init(j);
      j += 1
    }
  }

  private def case1(): Int = {
    // case1
    val VRm1Under = STATE(state_i + R - 1);
    val VRm2Under = STATE(state_i + R - 2);
    val V0 = STATE(state_i);
    val VM1 = STATE(state_i + M1);
    val VM2 = STATE(state_i + M2);
    val VM3 = STATE(state_i + M3);

    // state_i == 0
    z0 = (VRm1Under & MASKL) | (VRm2Under & MASKU);
    z1 = MAT0NEG(-25, V0) ^ MAT0POS(27, VM1);
    z2 = MAT3POS(9, VM2) ^ MAT0POS(1, VM3);
    STATE(state_i) = z1 ^ z2;
    STATE(state_i - 1 + R) = (z0) ^ MAT0NEG(-9, z1) ^ MAT0NEG(-21, z2) ^ MAT0POS(21, STATE(state_i));
    state_i = R - 1;
    caseNumber = 3;
    if (tempering) {
      var y = STATE(state_i) ^ ((STATE(state_i) << 7) & TEMPERB);
      y = y ^ ((y << 15) & TEMPERC);
      return (y);
    } else {
      return STATE(state_i);
    }
  }

  private def case2(): Int = {
    // case1
    val VRm1 = STATE(state_i - 1);
    val VRm2Under = STATE(state_i + R - 2);
    val V0 = STATE(state_i);
    val VM1 = STATE(state_i + M1);
    val VM2 = STATE(state_i + M2);
    val VM3 = STATE(state_i + M3);

    // state_i == 0
    z0 = (VRm1 & MASKL) | (VRm2Under & MASKU);
    z1 = MAT0NEG(-25, V0) ^ MAT0POS(27, VM1);
    z2 = MAT3POS(9, VM2) ^ MAT0POS(1, VM3);
    STATE(state_i) = z1 ^ z2;
    STATE(state_i - 1) = (z0) ^ MAT0NEG(-9, z1) ^ MAT0NEG(-21, z2) ^ MAT0POS(21, STATE(state_i));
    state_i = 0;
    caseNumber = 1;
    if (tempering) {
      var y = STATE(state_i) ^ ((STATE(state_i) << 7) & TEMPERB);
      y = y ^ ((y << 15) & TEMPERC);
      return (y);
    } else {
      return STATE(state_i);
    }

  }

  private def case3(): Int = {
    // case1
    val VRm1 = STATE(state_i - 1);
    val VRm2 = STATE(state_i - 2);
    val V0 = STATE(state_i);
    val VM1Over = STATE(state_i + M1 - R);
    val VM2Over = STATE(state_i + M2 - R);
    val VM3Over = STATE(state_i + M3 - R);

    z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
    z1 = MAT0NEG(-25, V0) ^ MAT0POS(27, VM1Over);
    z2 = MAT3POS(9, VM2Over) ^ MAT0POS(1, VM3Over);
    STATE(state_i) = z1 ^ z2;
    STATE(state_i - 1) = (z0) ^ MAT0NEG(-9, z1) ^ MAT0NEG(-21, z2) ^ MAT0POS(21, STATE(state_i));
    state_i -= 1;
    if (state_i + M1 < R) {
      caseNumber = 5;
    }
    if (tempering) {
      var y = STATE(state_i) ^ ((STATE(state_i) << 7) & TEMPERB);
      y = y ^ ((y << 15) & TEMPERC);
      return (y);
    } else {
      return STATE(state_i);
    }
  }

  private def case4(): Int = {
    // case1
    val VRm1 = STATE(state_i - 1);
    val VRm2 = STATE(state_i - 2);
    val V0 = STATE(state_i);
    val VM1 = STATE(state_i + M1);
    val VM2 = STATE(state_i + M2);
    val VM3Over = STATE(state_i + M3 - R);

    // state_i == 0
    z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
    z1 = MAT0NEG(-25, V0) ^ MAT0POS(27, VM1);
    z2 = MAT3POS(9, VM2) ^ MAT0POS(1, VM3Over);
    STATE(state_i) = z1 ^ z2;
    STATE(state_i - 1) = (z0) ^ MAT0NEG(-9, z1) ^ MAT0NEG(-21, z2) ^ MAT0POS(21, STATE(state_i));
    state_i -= 1;
    if (state_i + M3 < R) {
      caseNumber = 6;
    }
    if (tempering) {
      var y = STATE(state_i) ^ ((STATE(state_i) << 7) & TEMPERB);
      y = y ^ ((y << 15) & TEMPERC);
      return (y);
    } else {
      return STATE(state_i);
    }
  }

  private def case5(): Int = {
    // case1
    val VRm1 = STATE(state_i - 1);
    val VRm2 = STATE(state_i - 2);
    val V0 = STATE(state_i);
    val VM1 = STATE(state_i + M1);
    val VM2Over = STATE(state_i + M2 - R);
    val VM3Over = STATE(state_i + M3 - R);

    // state_i == 0
    z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
    z1 = MAT0NEG(-25, V0) ^ MAT0POS(27, VM1);
    z2 = MAT3POS(9, VM2Over) ^ MAT0POS(1, VM3Over);
    STATE(state_i) = z1 ^ z2;
    STATE(state_i - 1) = (z0) ^ MAT0NEG(-9, z1) ^ MAT0NEG(-21, z2) ^ MAT0POS(21, STATE(state_i));
    state_i -= 1;
    if (state_i + M2 < R) {
      caseNumber = 4;
    }
    if (tempering) {
      var y = STATE(state_i) ^ ((STATE(state_i) << 7) & TEMPERB);
      y = y ^ ((y << 15) & TEMPERC);
      return (y);
    } else {
      return STATE(state_i);
    }
  }

  private def case6(): Int = {
    // case1
    val VRm1 = STATE(state_i - 1);
    val VRm2 = STATE(state_i - 2);
    val V0 = STATE(state_i);
    val VM1 = STATE(state_i + M1);
    val VM2 = STATE(state_i + M2);
    val VM3 = STATE(state_i + M3);

    // state_i == 0
    z0 = (VRm1 & MASKL) | (VRm2 & MASKU);
    z1 = MAT0NEG(-25, V0) ^ MAT0POS(27, VM1);
    z2 = MAT3POS(9, VM2) ^ MAT0POS(1, VM3);
    STATE(state_i) = z1 ^ z2;
    STATE(state_i - 1) = (z0) ^ MAT0NEG(-9, z1) ^ MAT0NEG(-21, z2) ^ MAT0POS(21, STATE(state_i));
    state_i -= 1;
    if (state_i == 1) {
      caseNumber = 2;
    }
    if (tempering) {
      var y = STATE(state_i) ^ ((STATE(state_i) << 7) & TEMPERB);
      y = y ^ ((y << 15) & TEMPERC);
      return (y);
    } else {
      return STATE(state_i);
    }
  }

  def nextInt32(): Int = {
    caseNumber match {
      case 1 => return case1();
      case 2 =>
        return case2();
      case 3 =>
        return case3();
      case 4 =>
        return case4();
      case 5 =>
        return case5();
      case _ =>
        return case6();
    }

  }

  def nextDouble(): Double = {
    return (nextInt32() & MASK_UNSIGNED_INT) * NORM
  }
  
  def nextDoubleOpen(): Double = {
    ((nextInt32() & MASK_UNSIGNED_INT) + 0.5) * NORM
  }
}
package quantscale.random

/**
 * Mersenne Twister on 64 bits. 
 * This allows fast generation of random numbers with machine epsilon precision.
 * 
 * The code is a direct port from Matsumoto C code.
 */
class MersenneTwister19937x64 extends RandomGeneratorDouble {
  private val NN = 312
  private val MM = 156
  private val MATRIX_A = 0xB5026F5AA96619E9L
  private val UM = 0xFFFFFFFF80000000L /* Most significant 33 bits */
  private val LM = 0x7FFFFFFFL /* Least significant 31 bits */

  /* The array for the state vector */
  private var mt = Array.ofDim[Long](NN)

  private val mag01 = Array(0L, MATRIX_A)

  /* mti==NN+1 means mt[NN] is not initialized */
  private var mti = NN + 1

  /* initializes mt[NN] with a seed */
  def init(seed: Long) {
    mt(0) = seed
    mti = 1
    while (mti < NN) {
      mt(mti) = (6364136223846793005L * (mt(mti - 1) ^ (mt(mti - 1) >>> 62)) + mti)
      mti += 1
    }
  }

  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  def initByArray64(init_key: Array[Long],
                    key_length: Long) {
    init(19650218L);
    var i = 1; var j = 0;
    var k = if (NN > key_length) NN else key_length;
    while (k > 0) {
      mt(i) = (mt(i) ^ ((mt(i - 1) ^ (mt(i - 1) >>> 62)) * 3935559000370003845L)) + init_key(j) + j /* non linear */
      i += 1
      j += 1
      if (i >= NN) { mt(0) = mt(NN - 1); i = 1; }
      if (j >= key_length) j = 0;
      k -= 1
    }
    k = NN - 1
    while (k > 0) {
      mt(i) = (mt(i) ^ ((mt(i - 1) ^ (mt(i - 1) >>> 62)) * 2862933555777941757L)) - i /* non linear */
      i += 1
      if (i >= NN) { mt(0) = mt(NN - 1); i = 1; }
      k -= 1
    }

    mt(0) = 1L << 63; /* MSB is 1; assuring non-zero initial array */
  }

  /* generates a random number on [0, 2^64-1]-interval */
  def nextLong(): Long =
    {
      var i = 0
      var x: Long = 0

      if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN + 1)
          init(5489L);
        i = 0
        while (i < NN - MM) {
          x = (mt(i) & UM) | (mt(i + 1) & LM);
          mt(i) = mt(i + MM) ^ (x >>> 1) ^ mag01((x & 1L).toInt);
          i += 1
        }
        while (i < NN - 1) {
          x = (mt(i) & UM) | (mt(i + 1) & LM);
          mt(i) = mt(i + (MM - NN)) ^ (x >>> 1) ^ mag01((x & 1L).toInt);
          i += 1
        }
        x = (mt(NN - 1) & UM) | (mt(0) & LM);
        mt(NN - 1) = mt(MM - 1) ^ (x >>> 1) ^ mag01((x & 1L).toInt);

        mti = 0;
      }

      x = mt(mti);
      mti += 1

      x ^= (x >>> 29) & 0x5555555555555555L;
      x ^= (x << 17) & 0x71D67FFFEDA60000L;
      x ^= (x << 37) & 0xFFF7EEE000000000L;
      x ^= (x >>> 43);

      return x;
    }

  /* generates a random number on [0, 2^63-1]-interval */
  def nextLong63: Long =
    {
      return (nextLong() >>> 1);
    }

  /* generates a random number on [0,1]-real-interval */
  def nextDoubleClosed(): Double =
    {
      return (nextLong() >>> 11) * (1.0 / 9007199254740991.0);
    }

  /* generates a random number on [0,1)-real-interval */
  def nextDouble(): Double =
    {
      return (nextLong() >>> 11) * (1.0 / 9007199254740992.0);
    }

  /* generates a random number on (0,1)-real-interval */
  def nextDoubleOpen(): Double =
    {
      return ((nextLong() >>> 12) + 0.5) * (1.0 / 4503599627370496.0);
    }

}
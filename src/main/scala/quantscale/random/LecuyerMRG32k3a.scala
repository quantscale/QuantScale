package quantscale.random

/**
 * MRG32k3a random number generator. A uniform random number generator that can
 * skip ahead.
 *
 * See P. l'Ecuyer. "Good parameter sets for combined multiple recursive random
 * number generators." Operations Research, 47(1):159-164, 1999.
 *
 * The skip ahead logic follows from P. l'Ecuyer "AN OBJECT-ORIENTED
 * RANDOM-NUMBER PACKAGE WITH MANY LONG STREAMS AND SUBSTREAMS" Operations
 * Research, 2001.
 *
 */
class LecuyerMRG32k3a extends RandomGeneratorDouble {

  private val m1: Long = 4294967087l
  private val m2: Long = 4294944443l
  private val a12: Long = 1403580l
  private val a13n: Long = 810728l
  private val a21: Long = 527612l
  private val a23n: Long = 1370589l
  private val norm = 2.328306549295727688e-10

  /**
   * The same seed is used for every instance. Use skipTo if you want an
   * independent instance.
   */
  var s0: Long = 12345
  var s1: Long = 12345
  var s2: Long = 12345
  var s3: Long = 12345
  var s4: Long = 12345
  var s5: Long = 12345

  def nextLong32Open(): Long = {
    // this code is directly from P. l'Ecuyer. Good parameter sets for
    // combined multiple recursive random
    // number generators. (PS file) Operations Research, 47(1):159-164, 1999

    // replaced double with long as long is faster on 64-bit machines.
    /* Component 1 */
    var p1: Long = a12 * s1 - a13n * s0;
    var k: Int = (p1 / m1).toInt;
    p1 -= k * m1;
    if (p1 < 0.0)
      p1 += m1;
    s0 = s1;
    s1 = s2;
    s2 = p1;
    /* Component 2 */
    var p2: Long = a21 * s5 - a23n * s3;
    k = (p2 / m2).toInt;
    p2 -= k * m2;
    if (p2 < 0.0)
      p2 += m2;
    s3 = s4;
    s4 = s5;
    s5 = p2;
    /* Combination */
    return if (p1 > p2) (p1 - p2) else (p1 - p2 + m1);
  }

  def nextLong32(): Long = {
    // this code is directly from P. l'Ecuyer. Good parameter sets for
    // combined multiple recursive random
    // number generators. (PS file) Operations Research, 47(1):159-164, 1999

    // replaced double with long as long is faster on 64-bit machines.
    /* Component 1 */
    var p1: Long = a12 * s1 - a13n * s0;
    var k: Int = (p1 / m1).toInt;
    p1 -= k * m1;
    if (p1 < 0.0)
      p1 += m1;
    s0 = s1;
    s1 = s2;
    s2 = p1;
    /* Component 2 */
    var p2: Long = a21 * s5 - a23n * s3;
    k = (p2 / m2).toInt;
    p2 -= k * m2;
    if (p2 < 0.0)
      p2 += m2;
    s3 = s4;
    s4 = s5;
    s5 = p2;
    /* Combination */
    //Quantscale modification : allow 0 - we replace p1 > p2, with p1 >= p2
    return if (p1 >= p2) (p1 - p2) else (p1 - p2 + m1);

  }

  def nextInt32(): Int = {
    return nextLong32().toInt;
  }

  /**
   * This is not thread-safe: use skipTo() and multiple instances to
   * distribute accross threads.
   *
   * @return a uniformly distributed random number in [0,1)
   */
  def nextDouble(): Double = {
    return nextLong32() * norm;
  }

  /**
   * This is not thread-safe: use skipTo() and multiple instances to
   * distribute accross threads.
   *
   * @return a uniformly distributed random number in (0,1)
   */
  def nextDoubleOpen(): Double = {
    return nextLong32Open() * norm;
  }
}
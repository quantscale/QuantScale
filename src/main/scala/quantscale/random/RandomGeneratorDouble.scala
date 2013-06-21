package quantscale.random

trait RandomGeneratorDouble {
  /**
   * return a double in 0.0 (inclusive) to 1.0 (exclusive)
   */
  def nextDouble(): Double

  /**
   * return a double in 0.0 (exclusive) to 1.0 (exclusive)
   */
  def nextDoubleOpen(): Double
}
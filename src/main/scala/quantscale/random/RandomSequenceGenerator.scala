package quantscale.random

trait RandomSequenceGenerator {
  /**
   * return a list of uniform random numbers in (0,1) (0 and 1 excluded)
   */
  def nextSequence(sequence: Array[Double])

  def dimension(): Int
}
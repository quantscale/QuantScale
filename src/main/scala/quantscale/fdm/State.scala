package quantscale.fdm

import java.util.Arrays

class State(val stateDimensions: Int, val size: Int, var values: Array[Array[Double]]) {
  def price: Array[Double] = { return values(0) }
  def price_=(v: Array[Double]) { values(0) = v }

  def this(stateDimensions: Int, size: Int) = this(stateDimensions, size, Array.ofDim[Double](stateDimensions, size))

  override def toString(): String = {
    var s = ""
    for (i <- 0 until stateDimensions) {
      s += "state " + i + Arrays.toString(values(i))
      s += "\n"
    }
    return s
  }

  def copy(): State = {
    val s = new State(stateDimensions, size)
    var i = 0
    while (i < stateDimensions) {
      System.arraycopy(values(i), 0, s.values(i), 0, size)
      i+=1
    }
    return s
  }
}
package quantscale.mc

class MCPath(val time: Array[Double], val dimension: Int) {
  private val _values = Array.ofDim[Double](time.length, dimension)

  def values = _values

}

class AdjointMCPath(override val time: Array[Double], override val dimension: Int) extends MCPath(time, dimension) {
  private val _values_b = Array.ofDim[Double](time.length, dimension)

  def adjointValues = _values_b
  private val _delta = Array.ofDim[Double](time.length, dimension)
  private val _vega = Array.ofDim[Double](time.length, dimension)
  def delta = _delta
  def vega = _vega

}
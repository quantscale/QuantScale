package quantscale.mc

class MCPath(val time: Array[Double], val dimension: Int) {
  private var _values = Array.ofDim[Double](time.length, dimension)

  def values = _values
  
}
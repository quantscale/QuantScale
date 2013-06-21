package quantscale.mc

trait BrownianTransformFactory {
    def makeBrownianTransform(t : Array[Double]) : BrownianTransform
}

trait BrownianTransform {
  def transform(input: Array[Double], dimension: Int,
                output: Array[Array[Double]])
}
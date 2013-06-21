package quantscale.mc

trait PathGenerator {
  def init(evaluationTimes: Array[Double])
  def nextPath(): MCPath
}

trait PathPricer {
  def eval(path: MCPath): Double
  def evaluationTimes: Array[Double]
}
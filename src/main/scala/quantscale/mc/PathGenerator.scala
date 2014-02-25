package quantscale.mc

trait PathGenerator {
  def init(evaluationTimes: Array[Double])
  def nextPath(): MCPath
  def nextAdjointPath() : AdjointMCPath
  def assetDimension() : Int
}

trait PathPricer {
  def eval(path: MCPath): Double
  def evalAdjoint(path: AdjointMCPath) : Double
  def evaluationTimes: Array[Double]
}
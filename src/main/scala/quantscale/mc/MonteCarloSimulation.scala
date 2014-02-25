package quantscale.mc

import quantscale.random.RandomSequenceGenerator
import quantscale.mc.payoff.MCPayoff

class MonteCarloSimulation(pathGenerator: PathGenerator, pathPricer: PathPricer) {
  def evaluateForward(numSimulations: Int): Double = {
    var i = 0
    var sum = 0.0
    var path: MCPath = null
    pathGenerator.init(pathPricer.evaluationTimes)
    while (i < numSimulations) {
      path = pathGenerator.nextPath()
      sum += pathPricer.eval(path)
      i += 1
    }
    return sum / numSimulations
  }

  def evaluateAdjointForward(numSimulations: Int): (Double, Array[Double]) = {
    var i = 0
    var path: AdjointMCPath = null
    pathGenerator.init(pathPricer.evaluationTimes)
    var sum = 0.0
    val deltaSum = Array.ofDim[Double](pathGenerator.assetDimension())
    var iAsset = 0
    while (i < numSimulations) {
      path = pathGenerator.nextAdjointPath()
      sum += pathPricer.evalAdjoint(path)
      val path_b = path.adjointValues
      var tIndex = path_b.length-1
      val pathDelta = path.delta
      while (tIndex >= 0) {
        iAsset = deltaSum.length-1
        while (iAsset >= 0 ) {
          deltaSum(iAsset) += path_b(tIndex)(iAsset)*pathDelta(tIndex)(iAsset)
          iAsset -= 1
        }
        tIndex -= 1
      }
      i += 1
    }
    iAsset = 0
    while (iAsset < deltaSum.length) {
      deltaSum(iAsset) /= numSimulations
      iAsset +=1
    }
    return (sum / numSimulations,deltaSum)
  }

}
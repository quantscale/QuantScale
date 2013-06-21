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

}
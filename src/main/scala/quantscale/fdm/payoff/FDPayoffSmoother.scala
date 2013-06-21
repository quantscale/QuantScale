package quantscale.fdm.payoff

abstract class FDPayoffSmoother {
  def makeSmooth(v : FDPayoff);

}
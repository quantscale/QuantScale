package quantscale.fdm.payoff

abstract class FDPayoffSmoother {
  def makeSmooth(v: FDPayoff);

}

class NoFDPayoffSmoother extends FDPayoffSmoother {
  def makeSmooth(v: FDPayoff) {}
}
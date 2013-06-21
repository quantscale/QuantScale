package quantscale.fdm.payoff

abstract class AmericanFDPayoff() extends FDPayoff() {
    def lowerBound : Array[Double];
}
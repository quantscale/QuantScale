package quantscale.fdm.sabr

trait SABRDensitySolver {
  def solve()

  def price(isCall: Boolean, strike: Double): Double
}

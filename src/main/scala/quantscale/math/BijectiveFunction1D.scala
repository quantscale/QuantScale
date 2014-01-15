package quantscale.math

trait BijectiveFunction1D extends Function1D {
  def inverseValue(y: Double): Double
}

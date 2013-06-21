package quantscale.fdm {

  object Epsilon {
    val MACHINE_EPSILON = 2.2204460492503131e-016
    val MACHINE_EPSILON_SQRT = Math.sqrt(MACHINE_EPSILON)
    val MACHINE_EPSILON_FOURTH_ROOT = Math.sqrt(MACHINE_EPSILON_SQRT)
  }
}
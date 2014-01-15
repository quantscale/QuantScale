package quantscale.fdm.method

import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.OperatorLine
import quantscale.fdm.DifferentialCache
import quantscale.fdm.State

/**
 * Solves the log transformed problem: the space contains log(assetPrice)
 */
class TRBDF2CentralParabolic1DMethod(payoff: FDPayoff = null) extends TRBDF2Parabolic1DMethod(payoff) {

  override def initLeftHandSide(t: Double, dt: Double) {
    computeDiscount(ex, t, dt, discountVector)
    computeDrift(ex, t, dt, driftVector)
    computeVariance(ex, t, dt, driftVector, varianceVector)
    val multiplier = -backcoeff * dt
    val A = tridiagonal
    A.fill(0.0)
    A.plusD2(1, ex.length - 1, diffCache, varianceVector, 0.5 * multiplier)
      .plusD1Central(1, ex.length - 1, diffCache, driftVector, multiplier)
      .plusD0(1, ex.length - 1, discountVector, -multiplier, 1.0)
    //    A.parabolicOperator(1, x.length - 1, diffCache, varianceVector, 0.5 * multiplier, driftVector, multiplier, 1.0 - multiplier * r)
  }

}
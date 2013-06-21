package quantscale.fdm.method
import quantscale.fdm.TridiagonalSolver
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State

/**
 * PDE solver for dX = mu dt + sigma dW ( use mu*X  sigma*X for for lognormal)
 */
abstract class Parabolic1DMethod {
  
  def spec : Parabolic1DFDSpec 
  
  def initSystem(specV: Parabolic1DFDSpec)

  
  def solve(currentTime: Double, dt: Double, f: State)

  var solver : TridiagonalSolver = null
  
  var smearingReducer: VarianceFilter = new ScharfetterGummelVarianceFilter

  var lowerBoundary : Parabolic1DBoundaryFactory = ForwardLinearOrder1Parabolic1DBoundaryFactory 
  var upperBoundary : Parabolic1DBoundaryFactory = BackwardLinearOrder1Parabolic1DBoundaryFactory 

  def copy() : Parabolic1DMethod
  
  def shutdown() {}
  
  def computeDrift(x: Array[Double], t: Double, dt: Double, driftVector: Array[Double]) {
    if (!spec.bIsStateDependent) {
      val mu0 = spec.b(t, dt, 0)
      for (i <- 0 until x.length) {
        driftVector(i) = mu0
      }
    } else {
      for (i <- 0 until x.length) {
        driftVector(i) = spec.b(t, dt, x(i))
      }
    }
  }
  
  def computeDiscount(x: Array[Double], t: Double, dt: Double, discountVector: Array[Double]) {
    if (!spec.cIsStateDependent) {
      val mu0 = -spec.c(t, dt, 0)
      for (i <- 0 until x.length) {
        discountVector(i) = mu0
      }
    } else {
      for (i <- 0 until x.length) {
        discountVector(i) = -spec.c(t, dt, x(i))
      }
    }
  }

  def computeVariance(x: Array[Double], t: Double, dt: Double, driftVector: Array[Double], varianceVector: Array[Double]) {
    if (smearingReducer != null) {
      for (i <- 1 until x.length - 1) {
        var sigmasq = spec.a(t, dt, x(i))
        sigmasq = smearingReducer.filter(2 * sigmasq, driftVector(i), x(i) - x(i - 1), x(i + 1) - x(i))
        varianceVector(i) = sigmasq
      }
      var sig0 = spec.a(t, dt, x(0))
      sig0 = smearingReducer.filter(2 * sig0, driftVector(0), x(1) - x(0), x(1) - x(0))
      varianceVector(0) = sig0
      sig0 = spec.a(t, dt, x(x.length - 1))
      sig0 = smearingReducer.filter(2 * sig0, driftVector(x.length - 1), x(x.length - 1) - x(x.length - 2), x(x.length - 1) - x(x.length - 2))
      varianceVector(x.length - 1) = sig0
    } else {
      for (i <- 0 until x.length) {
        var sigma = spec.a(t, dt, x(i))
        varianceVector(i) = 2 * sigma
      }
    }
  }

}


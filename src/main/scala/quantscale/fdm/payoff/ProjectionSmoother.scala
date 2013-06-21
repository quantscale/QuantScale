package quantscale.fdm.payoff
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.TridiagonalSolver
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.Epsilon
import java.util.Arrays

//supposes that discontinuities are aligned with the space points
class ProjectionSmoother(discontinuities: Array[Double]) extends FDPayoffSmoother {

  def makeSmooth(payoff: FDPayoff) {
    val space = payoff.space
    val size = payoff.space.length

    val originalState = payoff.state.price;
    val originalSpace = payoff.space;
    val ooSpace = payoff.originalSpace;

    var M = new TridiagonalMatrix(size)
    val F = new Array[Double](size)
    val isDiscontinuity = new Array[Boolean](size)
    var i = 1
    val factor = 1.0 / 6.0
    M.middle(0) = factor * 2 * (space(1) - space(0))
    M.upper(0) = factor * (space(1) - space(0))
    F(0) = integrateLower(payoff, space(0), space(1))
    var k = 0
    while (i < size - 1) {
      M.lower(i) = factor * (space(i) - space(i - 1))
      M.middle(i) = factor * 2 * (space(i + 1) - space(i - 1))
      M.upper(i) = factor * (space(i + 1) - space(i))

      while (k < discontinuities.length && space(i) > (discontinuities(k) + Epsilon.MACHINE_EPSILON_SQRT)) {
        k += 1
      }
      isDiscontinuity(i) = k < discontinuities.length && Math.abs(discontinuities(k) - space(i)) < Epsilon.MACHINE_EPSILON_SQRT
      F(i) = integrate(payoff, space(i - 1), space(i), space(i + 1), isDiscontinuity(i))
      i += 1
    }
    F(size - 1) = integrateUpper(payoff, space(size - 2), space(size - 1))
    M.middle(size - 1) = factor * 2 * (space(size - 1) - space(size - 2))
    M.lower(size - 1) = factor * (space(size - 1) - space(size - 2))
    val solver = new ThomasTridiagonalSolver()
    solver.init(size)
    //    val x : Array[Double] = new Array[Double]
//    println(M)
//    println("F="+Arrays.toString(F))
    solver.solve(M, F, originalState)

    payoff.space = originalSpace;
    payoff.originalSpace = ooSpace;
    payoff.state.price = originalState;
//    println("state="+Arrays.toString(payoff.state))
  }

  def integrateLower(payoff: FDPayoff, lower: Double, upper: Double): Double = {
    //TODO evaluate all payoff at Si, Si+0.5, Si+1 at the same time: Array size = 2*space.length-1
    val tmpSpace = Array(lower, 0.5 * (lower + upper))
    payoff.initState(tmpSpace)
    payoff.eval()
    val fmid = payoff.state.price(1)
    val flower = payoff.state.price(0)
    return (upper - lower) * 0.5 * (4.0 / 6.0 * fmid + 1.0 / 3.0 * flower)
  }

  def integrateUpper(payoff: FDPayoff, lower: Double, upper: Double): Double = {
    //TODO evaluate all payoff at Si, Si+0.5, Si+1 at the same time: Array size = 2*space.length-1
    val tmpSpace = Array(0.5 * (lower + upper), upper)
    payoff.initState(tmpSpace)
    payoff.eval()
    val fmid = payoff.state.price(0)
    val fupper = payoff.state.price(1)
    return (upper - lower) * 0.5 * (4.0 / 6.0 * fmid + 1.0 / 3.0 * fupper)
  }

  def integrate(payoff: FDPayoff, lower: Double, middle: Double, upper: Double, isDiscontinuity: Boolean): Double = {
    var integral = 0.0

    if (isDiscontinuity) {
      val eps = Epsilon.MACHINE_EPSILON_SQRT
      val tmpSpace: Array[Double] = Array((lower + middle) * 0.5, middle - eps, middle + eps, 0.5 * (middle + upper))
       payoff.initState(tmpSpace)
      payoff.eval()
      val flower = payoff.state.price(0)
      val fmidlow = payoff.state.price(1)
      val fmidup = payoff.state.price(2)
      val fupper = payoff.state.price(3)
//           println("tmpSpace="+Arrays.toString(tmpSpace)+" f="+flower+","+fmidlow+","+fupper)

      integral = (middle - lower) * 0.5 * (4.0 / 6.0 * flower + 1.0 / 3.0 * fmidlow) + (upper - middle) * 0.5 * (4.0 / 6.0 * fupper + 1.0 / 3.0 * fmidup)
    } else {
      val tmpSpace = Array((lower + middle) * 0.5, middle, 0.5 * (middle + upper))
      payoff.initState(tmpSpace)
      payoff.eval()
      val flower = payoff.state.price(0)
      val fmiddle = payoff.state.price(1)
      val fupper = payoff.state.price(2)
//      println("tmpSpace="+Arrays.toString(tmpSpace)+" f="+flower+","+fmiddle+","+fupper)
      integral = (middle - lower) * 0.5 * (4.0 / 6.0 * flower + 1.0 / 3.0 * fmiddle) + (upper - middle) * 0.5 * (4.0 / 6.0 * fupper + 1.0 / 3.0 * fmiddle)
    }

    return integral
  }
}
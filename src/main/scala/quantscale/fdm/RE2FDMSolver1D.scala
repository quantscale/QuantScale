package quantscale.fdm

import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.mesh.DividedMesh2D
import java.util.ArrayList
import scala.concurrent.forkjoin.ForkJoinPool
import java.util.concurrent.Callable
import java.util.concurrent.Executors

/**
 * Richardson Extrapolation in time of order 2 based on order 1 scheme
 */
class RE2FDMSolver1D(val originalSolver: FDMSolver1D) extends FDMSolver1D(originalSolver.spec, originalSolver.method, originalSolver.solver) {

  override def solve(payoff: FDPayoff) {
    originalSolver.solve(payoff)
    price = payoff.state.price.clone

    // half step
    var gridHalf = new DividedMesh2D(originalSolver.grid, 2)
    originalSolver.grid = gridHalf
    originalSolver.solve(payoff)
    var priceHalf = originalSolver.price.clone()
    var i = 0
    while (i < price.length) {
      price(i) = 2 * priceHalf(i) - price(i)
      i += 1
    }
  }
  
  def solveParallel(payoff: FDPayoff) {
    var solvers = new Array[FDMSolver1D](2)

    var fjPool = new ForkJoinPool()
    solvers(0) = originalSolver.copy()
    solvers(1) = originalSolver

    // half step
    var gridHalf = new DividedMesh2D(originalSolver.grid, 2)
    originalSolver.grid = gridHalf

    val list = new ArrayList[Callable[Array[Double]]]()
    for (d <- 0 until solvers.length) {
      val c = new Callable[Array[Double]]() {
        def call(): Array[Double] = {
          solvers(d).solve(payoff.copy())
          return solvers(d).price
        }
      }
      list.add(c)
    }
    var listResult = fjPool.invokeAll(list)

    var priceFull = listResult.get(0).get()
    var priceHalf = listResult.get(1).get()
    var i = 0
    price = priceFull
    while (i < price.length) {
      price(i) = 2 * priceHalf(i) - priceFull(i)
      i += 1
    }
  }

}
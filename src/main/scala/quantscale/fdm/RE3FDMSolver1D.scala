package quantscale.fdm

import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.mesh.DividedMesh2D

/**
 * Richardson Extrapolation in time of order 3 based on order 1 scheme
 */
class RE3FDMSolver1D(val originalSolver : FDMSolver1D) extends FDMSolver1D(originalSolver.spec, originalSolver.method, originalSolver.solver) {

  override def solve(payoff : FDPayoff) {
      originalSolver.solve(payoff)
      price = payoff.state.price.clone
      
      // half step
      var gridHalf = new DividedMesh2D(originalSolver.grid, 2)
      var gridThird = new DividedMesh2D(originalSolver.grid, 3)
      originalSolver.grid = gridHalf
      originalSolver.solve(payoff)
      var priceHalf = originalSolver.price.clone()
      
       // third step
      originalSolver.grid = gridThird
      originalSolver.solve(payoff)
      var priceThird = originalSolver.price
      
      var i = 0
      while (i< price.length) {
        price(i) = 4.5*priceThird(i)-4*priceHalf(i)+0.5*price(i)
        i+=1
      }
  }
  
}
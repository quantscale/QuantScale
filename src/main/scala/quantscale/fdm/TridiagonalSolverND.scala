package quantscale.fdm

import java.util.ArrayList
import scala.concurrent.forkjoin.ForkJoinPool
import scala.collection.mutable.ArrayBuffer
import java.util.concurrent.Callable

class TridiagonalSolverND(private val solver1d: TridiagonalSolver, val n: Int) {

  private var solvers = new Array[TridiagonalSolver](n)
  private var fjPool = new ForkJoinPool()
  for (d <- 0 until n) {
    solvers(d) = solver1d.copy()
  }

  def solve(m: TridiagonalMatrix, d: Array[Array[Double]], x: Array[Array[Double]]) = solveSerial(m, d, x)
  
  def solveSerial(m: TridiagonalMatrix, d: Array[Array[Double]], x: Array[Array[Double]]) {
    //      val executor = Executors.newFixedThreadPool(8)
    //      
    //    for (i <- 0 until n) {
    //      executor.submit(new Runnable() { def run() {solvers(i).solve(m, d(i), x(i))}})
    //    }
    //      executor.shutdown()
    //      executor.awaitTermination(100, TimeUnit.SECONDS)
    for (i <- 0 until n) {
        solvers(i).solve(m, d(i), x(i))
    }
  }
  
  def solveParallel(m: TridiagonalMatrix, d: Array[Array[Double]], x: Array[Array[Double]]) {
    //      val executor = Executors.newFixedThreadPool(8)
    //      
    //    for (i <- 0 until n) {
    //      executor.submit(new Runnable() { def run() {solvers(i).solve(m, d(i), x(i))}})
    //    }
    //      executor.shutdown()
    //      executor.awaitTermination(100, TimeUnit.SECONDS)
    val list = new ArrayList[Callable[Object]]()
    for (i <- 0 until n) {
      list.add(new Callable[Object]() {
        def call(): Object = {
          solvers(i).solve(m, d(i), x(i))
          return null
        }
      })
    }
    fjPool.invokeAll(list)
  }
}
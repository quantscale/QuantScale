package quantscale.fdm.method;

import org.slf4j.LoggerFactory
import quantscale.fdm.TridiagonalMatrix
import java.util.Arrays
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.Parabolic1DFDSpec
import quantscale.fdm.State
import quantscale.fdm.OperatorLine
import quantscale.fdm.DifferentialCache
import quantscale.fdm.TridiagonalSolverND
import scala.concurrent.forkjoin.ForkJoinPool
import java.util.ArrayList
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.ExecutorService

/**
 * a must be between 1/2 and 2-sqrt(2) or over 2+sqrt(2)
 */
class TGAParabolic1DMethod(payoff: FDPayoff, private var _a: Double = 0.54) extends Parabolic1DMethod {
  final val logger = LoggerFactory.getLogger(getClass());

  private var fjPool = new ForkJoinPool(4)
  private var tridiagonal1, tridiagonal2: TridiagonalMatrix = null
  private var rhs1, rhs2: State = null

  private var specialIndex: Int = 0;
  private var r1, r2, s1, s2: Double = 0.0;
  private var _spec: Parabolic1DFDSpec = null;
  private var x: Array[Double] = null;
  var ex: Array[Double] = null;
  var driftVector: Array[Double] = null
  var varianceVector: Array[Double] = null
  var discountVector: Array[Double] = null
  var diffCache: DifferentialCache = null
  private var firstLine1, firstLine2: OperatorLine = null
  private var lastLine1, lastLine2: OperatorLine = null
  private var solverND1, solverND2: TridiagonalSolverND = null

  private var pool : ExecutorService = Executors.newFixedThreadPool(2)
  private var solverThread1, solverThread2 : SolverThread= null
  
  def copy(): Parabolic1DMethod = {
    return new TGAParabolic1DMethod(payoff, _a)
  }

  override def spec = _spec

  override def initSystem(specV: Parabolic1DFDSpec) {
    _spec = specV
    tridiagonal1 = new TridiagonalMatrix(spec.grid.spaceVector.size)
    tridiagonal2 = new TridiagonalMatrix(spec.grid.spaceVector.size)

    x = spec.grid.spaceVector;
    driftVector = new Array[Double](x.length)
    discountVector = new Array[Double](x.length)
    varianceVector = new Array[Double](x.length)
    ex = spec.grid.spaceTransform.transform(x);
    diffCache = new DifferentialCache(x)
    firstLine1 = new OperatorLine(0, 3)
    lastLine1 = new OperatorLine(x.length - 3, x.length)
    firstLine2 = new OperatorLine(0, 3)
    lastLine2 = new OperatorLine(x.length - 3, x.length)
    initFactors()
  }

  private def initFactors() {
    val sqrt = math.sqrt(_a * _a - 4 * _a + 2)
    r1 = (2 * _a - 1) / (_a + sqrt)
    r2 = (2 * _a - 1) / (_a - sqrt)
    s1 = (1 - _a + r1) / (r1 - r2)
    s2 = (1 - _a + r2) / (r2 - r1)
    //    println(_a+" "+sqrt+" "+r1+" "+r2+" "+s1+" "+s2)
  }

  def a_=(Value: Double): Unit = {
    _a = Value
    initFactors()
  }

  def initLeftHandSide(t: Double, dt: Double) {
    computeDiscount(ex, t, dt, discountVector)
    computeDrift(ex, t, dt, driftVector)
    computeVariance(ex, t, dt, driftVector, varianceVector)
    var A = tridiagonal1
    var m = -r1 * dt
    A.fill(0.0)
    A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5 * m)
      .plusD1(1, x.length - 1, diffCache, driftVector, 1.0 * m)
      .plusD0(1, x.length - 1, discountVector, -1.0 * m, 1.0)

    A = tridiagonal2
    m = -r2 * dt
    A.fill(0.0)
    A.plusD2(1, x.length - 1, diffCache, varianceVector, 0.5 * m)
      .plusD1(1, x.length - 1, diffCache, driftVector, 1.0 * m)
      .plusD0(1, x.length - 1, discountVector, -1.0 * m, 1.0)

    //    var i = 0
    //    val ratio = r2 / r1
    //    while (i < A.size) {
    //      tridiagonal2.lower(i) = tridiagonal1.lower(i) * ratio
    //      tridiagonal2.upper(i) = tridiagonal1.upper(i) * ratio
    //      tridiagonal2.middle(i) = 1 + (tridiagonal1.middle(i) - 1) * ratio
    //      i += 1
    //    }

  }

  def initBoundaries(t: Double, dt: Double, f: State) {
    val multiplier1 = -r1 * dt
    //other ode factory without +1 and mult
    lowerBoundary.makeLine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), multiplier1, firstLine1)
    val m = x.length - 1
    upperBoundary.makeLine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), multiplier1, lastLine1)

    tridiagonal1.setBoundaries(firstLine1, lastLine1)

    val multiplier2 = -r2 * dt
    //other ode factory without +1 and mult
    lowerBoundary.makeLine(0, x, varianceVector(0) * 0.5, driftVector(0), -discountVector(0), multiplier2, firstLine2)
    upperBoundary.makeLine(m, x, varianceVector(m) * 0.5, driftVector(m), -discountVector(m), multiplier2, lastLine2)

    tridiagonal2.setBoundaries(firstLine2, lastLine2)

  }

  def initRightHandSide(f: State) {
    if (rhs1 == null) {

      rhs1 = new State(f.stateDimensions, f.size)
      rhs2 = new State(f.stateDimensions, f.size)
    }
    for (d <- 0 until rhs1.stateDimensions) {
      val rhs1Values = rhs1.values(d)
      val rhs2Values = rhs2.values(d)
      val fValues = f.values(d)
      System.arraycopy(fValues, 0, rhs1Values, 0, x.length)
      System.arraycopy(fValues, 0, rhs2Values, 0, x.length)

      //      var j: Int = x.length - 1;
      //      while (j >= 0) {
      //        rhs1Values(j) = fValues(j)
      //        rhs2Values(j) = fValues(j)
      //        j -= 1
      //      }
    }
  }

  override def solve(currentTime: Double, dt: Double, f: State) {

    initLeftHandSide(currentTime, dt)
    initRightHandSide(f)
    initBoundaries(currentTime, dt, f)
    if (solverND1 == null) {
      solverND1 = new TridiagonalSolverND(solver, f.stateDimensions)
      solverND2 = new TridiagonalSolverND(solver, f.stateDimensions)
      solverThread1 = new SolverThread(solverND1, tridiagonal1, rhs1, rhs1)
      solverThread2 = new SolverThread(solverND2, tridiagonal2, rhs2, rhs2)
    }
//              solverND1.solve(tridiagonal1, rhs1.values, rhs1.values)
//              solverND2.solve(tridiagonal2, rhs2.values, rhs2.values)

//    val list = new ArrayList[Callable[Object]]()
//          list.add(new Callable[Object]() {
//            def call(): Object = {
//              solverND1.solve(tridiagonal1, rhs1.values, rhs1.values)
//              return null
//            }
//          })
//            list.add(new Callable[Object]() {
//            def call(): Object = {
//              solverND2.solve(tridiagonal2, rhs2.values, rhs2.values)
//              return null
//            }
//          })
//        val fList = fjPool.invokeAll(list)
//
//        var i = fList.size()-1
//        
//        while (i >= 0) {
//          fList.get(i).get()
//          i-=1
//        }
        
//    solverThread1 = new SolverThread(solverND1, tridiagonal1, rhs1, rhs1)
//      solverThread2 = new SolverThread(solverND2, tridiagonal2, rhs2, rhs2)
//         solverThread1.start()
//         solverThread2.start()
//         solverThread1.join()
//         solverThread2.join()
    val f1 = fjPool.submit(solverThread1)
    val f2 = fjPool.submit(solverThread2)
    f1.get()
    f2.get()
//    solverThread1.run()
//    solverThread2.run()
    
    for (d <- 0 until f.stateDimensions) {
      val fValues = f.values(d)
      val rhs1Values = rhs1.values(d)
      val rhs2Values = rhs2.values(d)
      var i = 0
      while (i < fValues.length) {
        fValues(i) = s1 * rhs1Values(i) + s2 * rhs2Values(i)
        i += 1
      }
    }
  }

  override def shutdown() {
    fjPool.shutdown()
    pool.shutdown
//    println("shutdown")
  }

}

class SolverThread(solverND: TridiagonalSolverND, var tridiagonal: TridiagonalMatrix, var rhs: State, var result: State) extends Thread {

  override def run() {
      solverND.solve(tridiagonal, rhs.values, result.values)
//      println(result.values(0)(100))
  }
}
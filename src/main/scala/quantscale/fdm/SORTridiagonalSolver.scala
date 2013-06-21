package quantscale.fdm
import payoff.AmericanFDPayoff
import payoff.FDPayoff
import java.util.Arrays

class SORTridiagonalSolver(payoff: FDPayoff, val absoluteTolerance : Double = 1e-8) extends TridiagonalSolver {

  
  val relativeTolerance = 0
  val maxIterations: Int = 1000
  val omega: Double = 1.0 //1.0 = GaussSiedel, > 1.0 overrelaxation. Converges for 0<omega<2
  var size = 0

  private var minDoubleArray: Array[Double] = null

  var thomasSolver = new ThomasTridiagonalSolver()

  def copy() : TridiagonalSolver = {
    val t = new ThomasTridiagonalSolver()
    val c = new SORTridiagonalSolver(payoff)
    c.thomasSolver = t
    c.init(size)
    return c
  }
  
  def init(sizeV: Int) {
    size = sizeV
    thomasSolver.init(sizeV)
    minDoubleArray = new Array[Double](size)
    Arrays.fill(minDoubleArray, Double.MinValue)
  }

  def lowerBound(): Array[Double] = {
    payoff match {
      case t: AmericanFDPayoff => return payoff.asInstanceOf[AmericanFDPayoff].lowerBound;
      case _                   => return minDoubleArray;
    }
  }

  def solve(m: TridiagonalMatrix, d: Array[Double], x: Array[Double]): Unit = {
    var loopCount = 0

    val lowerBoundV = lowerBound()

    thomasSolver.solve(m, d, x)
    var i = 0
    while (i < size) {
      x(i) = Math.max(x(i), lowerBoundV(i))
      i += 1
    }

    val a = m.lower
    val b = m.middle
    val c = m.upper
    val toleranceL2 = absoluteTolerance * absoluteTolerance * size
    var errorDenominator = 0.0
    var error = 0.0
    //TODO transform to while loop (esp for europeans) 
    do {
      loopCount += 1
      error = 0.0
      i = 0
      while (i < size) {
        var y = 0.0
        if (i == 0) {
          if (m.firstLine != null) {
            var back = 0.0
            for (j <- 1 until m.firstLine.iEnd) {
              back += m.firstLine.value(j) * x(j)
            }
            y = 1.0 / m.firstLine.value(i) * (d(i) - back)
          } else {
            val back = c(i) * x(i + 1)
            y = 1.0 / b(i) * (d(i) - back)
          }
        } else if (i == size - 1) {
          if (m.lastLine != null) {
            var back = 0.0
            for (j <- m.lastLine.iStart until i) {
              back += m.lastLine.value(j) * x(j)
            }
            y = 1.0 / m.lastLine.value(i) * (d(i) - back)
          } else {
            val back = a(i) * x(i - 1)
            y = 1.0 / b(i) * (d(i) - back)
          }
        } else {
          val back = c(i) * x(i + 1) + a(i) * x(i - 1)
          y = 1.0 / b(i) * (d(i) - back)
        }
        y = x(i) + omega * (y - x(i))
        y = math.max(y, lowerBoundV(i))
        val errorOne = y - x(i)
        error += errorOne * errorOne
        if (relativeTolerance >= 0) {
          errorDenominator += y * y;
        }
        x(i) = y
        i += 1
      }

    } while (error > toleranceL2 &&
      (error > errorDenominator * relativeTolerance) &&
      loopCount < maxIterations);

    if (loopCount == maxIterations) {
      throw new RuntimeException("Solver did not converge enough, error="
        + Math.sqrt(error / size));
    }
  }
}
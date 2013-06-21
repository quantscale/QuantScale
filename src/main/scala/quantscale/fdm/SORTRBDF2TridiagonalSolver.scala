package quantscale.fdm
import payoff.FDPayoff
import payoff.AmericanFDPayoff
import quantscale.fdm.method.TRBDF2SingleCentralBSM1FMethod

class SORTRBDF2TridiagonalSolver(payoff: FDPayoff) extends TridiagonalSolver {

  val tolerance = 1e-6
  val maxIterations: Int = 1000
  val omega: Double = 1.0 //1.0 = GaussSiedel, > 1.0 overrelaxation. Converges for 0<omega<2
  var size = 0
  var thomasSolver = new ThomasTridiagonalSolver()
  private var fTr : Array[Double] = null
  private var f,d2 : Array[Double] = null
  var method : TRBDF2SingleCentralBSM1FMethod = null
  
  def copy() : TridiagonalSolver = {
    val c = new SORTridiagonalSolver(payoff)
    c.init(size)
    return c
  }
  
  def init(sizeV: Int) {
    size = sizeV
    thomasSolver.init(sizeV)
    fTr = new Array[Double](sizeV)
    d2 = new Array[Double](sizeV)
    f = new Array[Double](sizeV)
  }

  def lowerBound(): Array[Double] = {
    payoff match {
      case t: AmericanFDPayoff => return payoff.asInstanceOf[AmericanFDPayoff].lowerBound;
      case _                   => return new Array[Double](size);
    }
  }

  def solve(m: TridiagonalMatrix, d: Array[Double], x: Array[Double]): Unit = {
    var loopCount = 0

    val lowerBoundV = lowerBound()

    thomasSolver.solve(m, d, fTr)
    var i = 0
    while (i < size) {  
      fTr(i) = Math.max(fTr(i), lowerBoundV(i))
      f(i) = x(i)
      i += 1
    }

    val a = m.lower
    val b = m.middle
    val c = m.upper
    val toleranceL2 = tolerance * tolerance * size

    var error = 0.0
    //TODO transform to while loop (esp for europeans) 
    do {
      loopCount += 1
      error = 0.0
      i = 0
      while (i < size) {
        var back = 0.0
        if (i == 0) {
          back = c(i) * fTr(i + 1)
        } else if (i == size - 1) {
          back = a(i) * fTr(i - 1)
        } else {
          back = c(i) * fTr(i + 1) + a(i) * fTr(i - 1)
        }
        var y = 1.0 / b(i) * (d(i) - back)
        y = fTr(i) + omega * (y - fTr(i))
        y = Math.max(y, lowerBoundV(i))
        fTr(i) = y
        i += 1
      }
      method.computeBDF2Rhs(fTr, f, d2);
      i = 0
      while (i < size) {
        var back = 0.0
        if (i == 0) {
          back = c(i) * x(i + 1)
        } else if (i == size - 1) {
          back = a(i) * x(i - 1)
        } else {
          back = c(i) * x(i + 1) + a(i) * x(i - 1)
        }
        var y = 1.0 / b(i) * (d2(i) - back)
        y = x(i) + omega * (y - x(i))
        y = Math.max(y, lowerBoundV(i))
        val errorOne = y - x(i)
        error += errorOne * errorOne
        x(i) = y
        i += 1
      }

    } while (error > toleranceL2 && loopCount < maxIterations);

    if (loopCount == maxIterations) {
      throw new RuntimeException("Solver did not converge enough, error="
        + Math.sqrt(error / size));
    }
  }
}
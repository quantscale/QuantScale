package quantscale.fdm

abstract class TridiagonalSolver {

  def copy() : TridiagonalSolver
  
  def init(size: Int)
  /**
   * 
   * Solve m.x = d
   * The unknown is x, and will be populated.
   */
  def solve(m: TridiagonalMatrix, d: Array[Double], x:Array[Double])
}
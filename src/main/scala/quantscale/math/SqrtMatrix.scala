package quantscale.math

import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.SingularValueDecomposition
import org.apache.commons.math3.linear.CholeskyDecomposition

trait SqrtMatrix {
def sqrt(m: RealMatrix): RealMatrix
}

class SVDSqrt {
  def sqrt(m: RealMatrix): RealMatrix = {
    val svd = new SingularValueDecomposition(m)
    val u = svd.getU()
    val s = svd.getSingularValues()
    var j = 0
    while (j < u.getRowDimension()) {
      var k = 0
      while (k < u.getColumnDimension()) {
        val sqrt = math.sqrt(math.max(0, s(k)))
        u.multiplyEntry(j, k, sqrt)
        k += 1
      }
      j += 1
    }
    return u
  }
}

class CholeskySqrt {
  def sqrt(m: RealMatrix): RealMatrix = {
    val decomp = new CholeskyDecomposition(m)
    return decomp.getL()
  }
}
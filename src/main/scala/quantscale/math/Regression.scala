package quantscale.math

import org.apache.commons.math3.linear.QRDecomposition
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector

trait Regression {
  def regress(x: Array[Array[Double]], y: Array[Double])
  def value(x: Array[Double]) : Double
}

class QRRegression extends Regression {
  var beta: RealVector = null

  def regress(x: Array[Array[Double]], y: Array[Double]) {
    var decomp = new QRDecomposition(new Array2DRowRealMatrix(x));
    beta = decomp.getSolver().solve(new ArrayRealVector(y));
  }

  def value(x: Array[Double]): Double = {
    var sum = beta.dotProduct(new ArrayRealVector(x))
    return sum;
  }
}
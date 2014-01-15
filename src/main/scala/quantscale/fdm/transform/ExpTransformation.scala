package quantscale.fdm.transform

import quantscale.math.BijectiveFunction1D

class ExpTransformation(specialPoint: Double) extends CoordinateTransformation with BijectiveFunction1D {

  def value(x: Double): Double = {
    return specialPoint * Math.exp(x)
  }

  def inverseValue(x: Double): Double = {
    return Math.log(x / specialPoint)
  }

  def transform(x: Array[Double]): Array[Double] = {
    var y = new Array[Double](x.length);
    var i = 0;
    while (i < x.length) {
      y(i) = specialPoint * Math.exp(x(i));
      i += 1;
    }
    return y;
  }

  def inverseTransform(x: Array[Double]): Array[Double] = {
    var y = new Array[Double](x.length);
    var i = 0;
    while (i < x.length) {
      y(i) = Math.log(x(i) / specialPoint);
      i += 1;
    }
    return y;
  }
}
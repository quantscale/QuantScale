package quantscale.fdm.transform

abstract class CoordinateTransformation {
  def transform(x: Array[Double]): Array[Double];
}

object IdentityTransformation extends CoordinateTransformation {
  def transform(x: Array[Double]): Array[Double] = {
    return x;
  }
}
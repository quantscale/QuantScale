package quantscale.fdm.transform

class ExpTransformation(specialPoint : Double) extends CoordinateTransformation {

  def transform(x: Array[Double]): Array[Double] = {
    var y = new Array[Double](x.length);
    var i = 0;
    while (i < x.length) {
      y(i) = specialPoint*Math.exp(x(i));
      i += 1;
    }
    return y;
  }
  
  def inverseTransform(x : Array[Double]): Array[Double] =  {
    var y = new Array[Double](x.length);
    var i = 0;
    while (i < x.length) {
      y(i) = Math.log(x(i)/specialPoint);
      i += 1;
    }
    return y;    
  }
}
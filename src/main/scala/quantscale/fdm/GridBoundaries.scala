package quantscale.fdm

object GridBoundaries {
  def makeBoundariesWithStdDev(tte: Double, strike: Double, vol: Double, driftDf: Double = 1.0, numberOfStdDev: Int = 5): GridBoundaries = {
    val e = Math.exp(numberOfStdDev * vol *Math.sqrt( tte));
    return new GridBoundaries(0, tte,
      strike / driftDf / e,
      strike / driftDf * e);
  }
}

class GridBoundaries(tMin: Double, tMax: Double, xMin: Double, xMax: Double) {

  var firstTime = tMin;
  var lastTime = tMax;
  var bottomSpace = xMin;
  var topSpace = xMax;

}
package quantscale.fdm.mesh

object MeshBoundaries {
 
  def makeBoundariesWithStdDev(tte: Double, strike: Double, vol: Double, driftDf: Double = 1.0, numberOfStdDev: Int = 5): MeshBoundaries = {
    val e = math.exp(numberOfStdDev * vol *Math.sqrt( tte));
    return new MeshBoundaries(0, tte,
      if (driftDf > 1) strike /driftDf / e else strike / e,
      if (driftDf > 1) strike * e else strike/driftDf*e);
  }
}

class MeshBoundaries(val firstTime: Double, val lastTime: Double, val bottomSpace: Double, val topSpace: Double) {
  
  val spaceBoundaries : Mesh1DBoundaries = new Mesh1DBoundaries(bottomSpace, topSpace)
  val timeBoundaries : Mesh1DBoundaries = new Mesh1DBoundaries(firstTime, lastTime)
}

class Mesh1DBoundaries(val min: Double, val max : Double) {
}
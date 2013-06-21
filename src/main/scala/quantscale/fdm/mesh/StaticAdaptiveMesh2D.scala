package quantscale.fdm.mesh

class StaticAdaptiveMesh2D(spaceTransform: Mesh1D, timeTransform: Mesh1D) extends Mesh2D {
  private var _lastTime = 0.0
  
  private val t = timeTransform.x
  _lastTime = t(t.length-1)
  
  override def copy() : Mesh2D = {
    return new StaticAdaptiveMesh2D(spaceTransform, timeTransform)
  }
  override def spaceSize = spaceTransform.size;
  override def spaceVector = spaceTransform.x;
  override def lastTime = _lastTime;
  override def timeIterator = t.reverse.iterator;
}
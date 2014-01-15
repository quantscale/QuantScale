package quantscale.fdm.mesh

class StaticAdaptiveMesh2D(spaceMesh: Mesh1D, timeTransform: Mesh1D) extends Mesh2D {
  private var _lastTime = 0.0

  private val t = timeTransform.x
  _lastTime = t(t.length - 1)

  override def copy(): Mesh2D = {
    return new StaticAdaptiveMesh2D(spaceMesh, timeTransform)
  }

  override def spaceSize = spaceMesh.size;

  override def spaceVector = spaceMesh.x;

  override def lastTime = _lastTime;

  override def timeIterator = t.reverse.iterator;
}
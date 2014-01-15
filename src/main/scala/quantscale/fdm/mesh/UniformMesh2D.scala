package quantscale.fdm.mesh {

import scala.collection.mutable.ArrayBuffer
import java.util.Arrays
import quantscale.fdm.Epsilon
import scala.util.Sorting

class UniformMesh2D(spaceSizeV: Int, timeSize: Int, bounds: MeshBoundaries, specialPoint: Double, isMiddle: Boolean = false) extends Mesh2D {

  private var _lastTime = 0.0

  private val timeMesh = new UniformMesh1D(timeSize, new Mesh1DBoundaries(bounds.firstTime, bounds.lastTime))
  private val t = timeMesh.x
  _lastTime = t(t.length - 1)

  private val spaceMesh = new UniformMesh1D(spaceSizeV, new Mesh1DBoundaries(bounds.bottomSpace, bounds.topSpace), specialPoint, isMiddle)

  override def spaceSize = spaceMesh.size;

  override def spaceVector = spaceMesh.x;

  override def lastTime = _lastTime;

  override def timeIterator = t.reverse.iterator;


  override def copy(): Mesh2D = {
    return new UniformMesh2D(spaceSizeV, timeSize, bounds, specialPoint, isMiddle)
  }

}

}
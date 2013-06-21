package quantscale.fdm.listener
import quantscale.fdm.Epsilon
import scala.collection.mutable.ArrayBuffer
import java.io.PrintWriter
import quantscale.fdm.State

class Price2DListener() extends FDMListener1D {
  var price: ArrayBuffer[Array[Double]] = new ArrayBuffer[Array[Double]](50)
  var x: Array[Double] = null
  var allT = new ArrayBuffer[Double](50)

  override def update(t: Double, f: State, x: Array[Double]): Unit = {
    if (this.x == null) this.x = x
    allT += t
    price += f.price.clone
  }

  def print(pw: PrintWriter) {
    var buf = new StringBuilder()
    for (j <- 0 until x.length) {
      buf.append(x(j))
      if (j != x.length - 1) {
        buf.append(" ")
      }
    }
    pw.println(buf.toString())
    for (i <- 0 until allT.length) {
      buf = new StringBuilder()
      for (j <- 0 until x.length) {
        buf.append(price(i)(j))
        if (j != x.length - 1) {
          buf.append(" ")
        }
      }
      pw.println(buf)
    }
  }
}
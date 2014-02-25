package quantscale.math

/**
 * solve y'(x)=f(x,y(x))
 */

class RungeKuttaODESolver {

  //y(n+1) = y(n)+1.0/6.0*(k1+2*k2+2*k3+k4)
  //k1 = h*f(x(n),y(n))
  //k2 = h*f(x(n)+h/2.0,y(n)+k1/2.0)
  //k3 = h*f(x(n)+h/2.0,y(n)+k2/2.0)
  //k4 = h*f(x(n)+h,y(n)+k3)

  //Boundary contains y0=y(x0)
  //x is x range, for use [0,1] discretized
  def solve(f: Function2D, x: Array[Double], y0: Double, y: Array[Double]) {
    var n = 0

    y(0) = y0
    val nMax = x.length - 1
    while (n < nMax) {
      val h = x(n + 1) - x(n)
      val k1 = h * f.value(x(n), y(n))
      val k2 = h * f.value(x(n) + h / 2.0, y(n) + k1 / 2.0)
      val k3 = h * f.value(x(n) + h / 2.0, y(n) + k2 / 2.0)
      val k4 = h * f.value(x(n) + h, y(n) + k3)
      y(n + 1) = y(n) + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
      n += 1
    }
  }
}
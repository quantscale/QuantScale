package quantscale.vol
import quantscale.random.LecuyerMRG32k3a
import quantscale.math._

object ButterflyArbitrageFinder extends App {

  class QuadraticPoint(var w: Double, var dwdy: Double, var d2wdy: Double) {
    def this() = this(0, 0, 0)

    override def toString(): String = {
      return w + " " + dwdy + " " + d2wdy
    }
  }

  def computeButterfly(y: Double, w: Double, dwdy: Double, d2wdy: Double): Double = {
    val winv = 1 / w
    return 1 - y * winv * (dwdy) + 0.25 * (-0.25 - winv + y * y * winv * winv) * dwdy * dwdy + 0.5 * d2wdy
  }

  def computeButterfly(y: Double, p: QuadraticPoint): Double = {
    return computeButterfly(y, p.w, p.dwdy, p.d2wdy)
  }

  def isIncreasing(w0: Double, w1: Double): Boolean = {
    return w1 >= w0
  }

  def hasButterflyArbitrage(y: Double, p: QuadraticPoint): Boolean = {
    return computeButterfly(y, p) <= 0
  }

  def interpolateLinearly(t: Double, w0: Double, w1: Double): Double = {
    return (1 - t) * w0 + t * w1
  }

  def interpolatePointLinearly(t: Double, p0: QuadraticPoint, p1: QuadraticPoint, pOut: QuadraticPoint) = {
    pOut.w = interpolateLinearly(t, p0.w, p1.w)
    pOut.dwdy = interpolateLinearly(t, p0.dwdy, p1.dwdy)
    pOut.d2wdy = interpolateLinearly(t, p0.d2wdy, p1.d2wdy)
  }

  def makeRandomPoint(rng: LecuyerMRG32k3a, pMin: QuadraticPoint, pMax: QuadraticPoint, pOut: QuadraticPoint) = {
    pOut.w = (pMax.w - pMin.w) * rng.nextDoubleOpen() + pMin.w
    pOut.dwdy = (pMax.dwdy - pMin.dwdy) * rng.nextDoubleOpen() + pMin.dwdy
    pOut.d2wdy = (pMax.d2wdy - pMin.d2wdy) * rng.nextDoubleOpen() + pMin.d2wdy
  }

  def findIt(spaceSize: Integer, timeSize: Integer, y: Double, pMin: QuadraticPoint, pMax: QuadraticPoint) = {

    val rng = new LecuyerMRG32k3a()
    var hasArbitrage = false
    var countSpace: Integer = 0
    var countSpaceOk: Integer = 0
    var t = 0.0
    val p0 = new QuadraticPoint()
    val p1 = new QuadraticPoint()
    val pTemp = new QuadraticPoint()
    val p = new QuadraticPoint()
    pTemp.w = pMin.w
    pTemp.d2wdy = pMin.d2wdy
    pTemp.dwdy = pMin.dwdy
    while (!hasArbitrage && countSpace < spaceSize) {
      makeRandomPoint(rng, pMin, pMax, p0)
      pTemp.w = p0.w //increasing variance
      makeRandomPoint(rng, pTemp, pMax, p1)
        p1.dwdy = 0
      p1.d2wdy = 0
      if (!hasButterflyArbitrage(y, p0) && !hasButterflyArbitrage(y, p1)) {
        var timeIndex = 1
        while (timeIndex < timeSize) {
          t = timeIndex.toDouble / timeSize
          interpolatePointLinearly(t, p0, p1, p)
          val winv = 1 / p.w
          //          val aaa = computeButterfly(y,p)-t*computeButterfly(y,p1)-(1-t)*computeButterfly(y,p0)
          //            // (1-t)*t*(p0.dwdy*p0.dwdy*p1.w/p0.w+p1.dwdy*p1.dwdy*p0.w/p1.w-2*p0.dwdy*p1.dwdy)
          //            hasArbitrage =
          //            aaa < 0 //linear convexity approx
          //          if (aaa < 0) Console.println(aaa);
          val aaa = hasButterflyArbitrage(y, p)
          if (aaa) {
            Console.println(t + " " + computeButterfly(y, p));
            Console.println(p);
          }
          hasArbitrage = hasArbitrage || aaa
          //second derivative
          //          val mid = computeButterfly(y,p)
          //          val midg = mid*p.w
          //          val eps = 1e-3
          //          interpolatePointLinearly(t+eps, p0, p1, p)
          //          val up = computeButterfly(y,p)
          //          val upg = up*p.w
          //          interpolatePointLinearly(t-eps, p0, p1, p)
          //          val down = computeButterfly(y,p)
          //          val downg = down*p.w
          //          val second = (up-2*mid+down)/eps/eps
          //          if (second > 0) {
          //              Console.println("found convexity at "+p+ " countok="+countSpaceOk+ " second="+second)
          //              Console.println(p0 +" "+ hasButterflyArbitrage(y,p0));
          //              Console.println(p1+" "+ hasButterflyArbitrage(y,p1));
          //              Console.println("t=" + t);
          //          }
          //          val secondg= (upg-2*midg+downg)/eps/eps
          //          if (secondg > 0) {
          //              Console.println("found g convexity  at "+p+ " countok="+countSpaceOk+ " second="+secondg)
          //              Console.println(p0 +" "+ hasButterflyArbitrage(y,p0));
          //              Console.println(p1+" "+ hasButterflyArbitrage(y,p1));
          //              Console.println("t=" + t);
          //          }
          timeIndex += 1
        }
        countSpaceOk += 1
        //        Console.println(countSpaceOk);
        //        
        //        Console.println(p0 +" "+ hasButterflyArbitrage(y,p0));
        //      Console.println(p1+" "+ hasButterflyArbitrage(y,p1));
      }
      countSpace += 1
    }
    Console.println("hasArbitrage=" + hasArbitrage + " countSpace=" + countSpace + " countSpaceOk=" + countSpaceOk)
    if (hasArbitrage) {
      Console.println(p0 + " " + hasButterflyArbitrage(y, p0));
      Console.println(p1 + " " + hasButterflyArbitrage(y, p1));
      Console.println("t=" + t);

    }
  }

  def findArbitrageInSpline(maxIterations: Integer, spaceSize: Integer, timeSize: Integer, yMin: Double, yMax: Double, wMin: Double, wMax: Double, numberOfNodes: Integer) {
    val rng = new LecuyerMRG32k3a()
    var hasArbitrage = false
    var countSpace: Integer = 0
    var countSpaceOk: Integer = 0

    var y = new Array[Double](numberOfNodes)
    for (i <- 0 until y.length) {
      y(i) = (yMax - yMin) * i / (y.length - 1) + yMin
    }

    while(countSpace < maxIterations) {
    var w0 = new Array[Double](numberOfNodes)
    var spline0 : CubicPP = null
    hasArbitrage = true
    while (hasArbitrage && countSpace < maxIterations) {
      for (i <- 0 until numberOfNodes) {
        w0(i) = rng.nextDoubleOpen() * (wMax - wMin) + wMin
      }
      spline0 = CubicSpline.makeCubicSpline(y, w0, new Natural(), new Natural())
      //verify that spline0 has no arbitrage between yMin and yMax
      var i: Integer = 0
      while (i < spaceSize && !hasArbitrage) {
        val currentY = (yMax - yMin) * i / (spaceSize - 1) + yMin
        val arbValue = computeButterfly(currentY, spline0.value(currentY), spline0.derivative(currentY), spline0.secondDerivative(currentY))
        hasArbitrage = hasArbitrage || arbValue < 0
        i += 1
      }
      countSpace += 1
    }
    if (!hasArbitrage) Console.println("spline ok "+w0)
    var w1 = new Array[Double](numberOfNodes)
    hasArbitrage = true
    var spline1: CubicPP = null
    while (hasArbitrage && countSpace < maxIterations) {
      for (i <- 0 until numberOfNodes) {
        w1(i) = rng.nextDoubleOpen() * (wMax - w0(i)) + w0(i)
      }
      spline1 = CubicSpline.makeCubicSpline(y, w1, new Natural(), new Natural())
      //verify that spline0 has no arbitrage between yMin and yMax
      var i: Integer = 0
      while (i < spaceSize && !hasArbitrage) {
        val currentY = (yMax - yMin) * i / (spaceSize - 1) + yMin
        val arbValue = computeButterfly(currentY, spline1.value(currentY), spline1.derivative(currentY), spline1.secondDerivative(currentY))
        hasArbitrage = hasArbitrage || arbValue < 0
        i += 1
      }
      countSpace += 1
    }
    if (!hasArbitrage) {
        countSpaceOk += 1
    
        var t = 1
        while (t < timeSize) {
            val currentTime = t/(timeSize.doubleValue())
            var i: Integer = 0
            while (i < spaceSize && !hasArbitrage) {
                val currentY = (yMax - yMin) * i / (spaceSize - 1) + yMin
                val wInterp = spline0.value(currentY)*(1-currentTime)+currentTime*spline1.value(currentY)
                val dwInterp = spline0.derivative(currentY)*(1-currentTime)+currentTime*spline1.derivative(currentY)
                val d2wInterp = spline0.secondDerivative(currentY)*(1-currentTime)+currentTime*spline1.secondDerivative(currentY)
                val arbValue = computeButterfly(currentY, wInterp, dwInterp, d2wInterp)
                hasArbitrage = hasArbitrage || arbValue < 0
                if (hasArbitrage) {
                  Console.println("countSpace="+countSpace+" countOk="+countSpaceOk)
                  Console.println(arbValue+" Arbitrage at t="+currentTime+" y="+currentY)
                  Console.println(w0)
                  Console.println(w1)
                  System.exit(0)
                }
                i += 1
            }
            t+=1
        }
    }
    }
                 Console.println("countSpace="+countSpace+" countOk="+countSpaceOk)

  }

  //w = w(yMin)+dw/dy(yMin)*(y-yMin)= a + b y => a and b are random
  def findArbitrageInLine(spaceSize: Integer, timeSize: Integer, yMin: Double, yMax: Double, pMin: QuadraticPoint, pMax: QuadraticPoint) = {

    val rng = new LecuyerMRG32k3a()
    var hasArbitrage = false
    var countSpace: Integer = 0
    var countSpaceOk: Integer = 0
    var t = 0.0
    val p0 = new QuadraticPoint()
    val p1 = new QuadraticPoint()
    val p = new QuadraticPoint()
    var y = 0.0
    while (!hasArbitrage && countSpace < spaceSize) {
      var w0 = rng.nextDoubleOpen() * (pMax.w - pMin.w) + pMin.w

//      var dw0=0.0
      var dw0 = rng.nextDoubleOpen() * (pMax.dwdy - pMin.dwdy) + pMin.dwdy
      if (dw0 < 0) {
        w0 -= dw0 * (yMax - yMin) //keep w0 > wMin over yMin, yMax
      }

      var w1 = rng.nextDoubleOpen() * (pMax.w - w0) + w0

      var dw1 = rng.nextDoubleOpen() * (pMax.dwdy - pMin.dwdy) + pMin.dwdy
      if (dw1 < 0) {
        w1 -= dw1 * (yMax - yMin) //keep w0 > wMin over yMin, yMax
      }
      var hasBoundaryArb = false
      var yIndex = 0
      while (yIndex <= 20 && !hasBoundaryArb) {
        y = yMin + (yMax - yMin) * yIndex / 20.0
        p0.w = w0 + (y - yMin) * dw0
        p0.dwdy = dw0
        p1.w = w1 + (y - yMin) * dw1
        p1.dwdy = dw1
        hasBoundaryArb = hasButterflyArbitrage(y, p0) || hasButterflyArbitrage(y, p1)
        yIndex += 1
      }
      if (!hasBoundaryArb) {
        var timeIndex = 1
        while (timeIndex < timeSize && !hasArbitrage) {
          t = timeIndex.toDouble / timeSize
          yIndex = 0
          while (yIndex <= 20 && !hasArbitrage) {
            y = yMin + (yMax - yMin) * yIndex / 20.0
            p0.w = w0 + (y - yMin) * dw0
            p0.dwdy = dw0
            p1.w = w1 + (y - yMin) * dw1
            p1.dwdy = dw1
            interpolatePointLinearly(t, p0, p1, p)
            hasArbitrage = hasButterflyArbitrage(y, p)
            if (hasArbitrage) {
              Console.println("found " + p);
              Console.println("w0=" + w0 + " dw0=" + dw0);
              Console.println("w1=" + w1 + " dw1=" + dw1);
              Console.println("t =" + t + " y=" + y);
            }
            yIndex += 1
          }
          timeIndex += 1
        }

        countSpaceOk += 1
      }
      countSpace += 1

    }
    Console.println("hasArbitrage=" + hasArbitrage + " countSpace=" + countSpace + " countSpaceOk=" + countSpaceOk)
    if (hasArbitrage) {
      Console.println(p0 + " " + hasButterflyArbitrage(y, p0));
      Console.println(p1 + " " + hasButterflyArbitrage(y, p1));
      Console.println("t =" + t + " y=" + y);

    }
  }

  val t0 = System.nanoTime()
    findIt(1000*1000*1, 10, 0.5, new QuadraticPoint(0.1,-10,0), new QuadraticPoint(10, 10, 0))
//  findArbitrageInLine(1000 * 1000 * 10, 10, 0.6, 0.8, new QuadraticPoint(0.1, -10, 0), new QuadraticPoint(10, 10, 0))
//  findArbitrageInSpline(1000*1000*10, 10, 5, 0.6, 0.8, 0.1, 1.0, 3)
  val t1 = System.nanoTime()
  println("Elapsed time: " + (t1 - t0) * 1e-9 + "s")
}
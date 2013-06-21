package test.quantscale.fdm
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import quantscale.fdm.listener.GammaListener
import quantscale.fdm.method.TRBDF2UncertainBSM1FMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.ThetaUncertainBSM1FMethod
import quantscale.fdm.payoff.ButterflyFDPayoff
import quantscale.fdm.payoff.DigitalFDPayoff
import quantscale.fdm.payoff.SimpsonIntegralSmoother
import quantscale.fdm.Epsilon
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.UncertainBSM1FFDSpec
import quantscale.math.JacobianConcentration
import quantscale.math.RungeKuttaODESolver
import org.scalatest.junit.JUnitRunner
import quantscale.fdm.mesh.UniformMesh2D
import quantscale.fdm.mesh.MeshBoundaries

@RunWith(classOf[JUnitRunner])
class UncertainVolatilityEulerSuite extends FunSuite {

  test("Butterfly Spread Gamma First Step 1") {
    var spaceSize = 400
    var timeSize = 20
    var spot = 100.0

    val tte = 0.25
    val r = 0.1
    val mu = 0.1
    val K1 = 90
    val K2 = 110
    val volMin = 0.15
    val volMax = 0.25

    val xMin = 0
    val xMax = 200

    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      new MeshBoundaries(0, tte, xMin, xMax),
      K1)

    val spec = new UncertainBSM1FFDSpec(
      grid,
      volMin,
      volMax,
      mu,
      r,
      false)

    val payoff = new ButterflyFDPayoff(true, K1, K2, tte);
    val smoother = null //new SimpsonIntegralSmoother(Array(K1,(K1+K2)/2,K2));
    val solver = new ThomasTridiagonalSolver(payoff);
    val method = new ThetaUncertainBSM1FMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
    val pricer = new FDMSolver1D(spec, method, solver);
    pricer.smoother = smoother;
    val timeIterator = grid.timeIterator
    timeIterator.next();
    val firstTime = timeIterator.next();
    assert(Math.abs(firstTime - (tte - tte / timeSize)) < Epsilon.MACHINE_EPSILON_SQRT, "first time step is not correct")
    val gammaListener = new GammaListener(Array(firstTime))
    pricer.listener = gammaListener;
    pricer.solve(payoff);
    val price = pricer.price(spot);
    val refPrice = 2.29769 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"
    println(" priceRef=" + refPrice + " price=" + price);

    gammaListener.print()
    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));

  }
  test("Butterfly Spread Price 1") {
    var spaceSize = 3600
    var timeSize = 20
    var spot = 100.0

    val tte = 0.25
    val r = 0.1
    val mu = 0.1
    val K1 = 90
    val K2 = 110
    val volMin = 0.15
    val volMax = 0.25

    val xMin = 0
    val xMax = 200

    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      new MeshBoundaries(0, tte, xMin, xMax),
      K1)

    val spec = new UncertainBSM1FFDSpec(
      grid,
      volMin,
      volMax,
      mu,
      r,
      false)

    val payoff = new ButterflyFDPayoff(true, K1, K2, tte);
    val smoother = null //new SimpsonIntegralSmoother(Array(K1,(K1+K2)/2,K2));
    val solver = new ThomasTridiagonalSolver(payoff);
    val method = new ThetaUncertainBSM1FMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
    val pricer = new FDMSolver1D(spec, method, solver);
    pricer.smoother = smoother;

    pricer.solve(payoff);
    val price = pricer.price(spot);
    val refPrice = 2.29769 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"
    println(" priceRef=" + refPrice + " price=" + price);
    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));

  }

  test("Butterfly Spread 1 Convergence") {
    var spaceSize = Array(60, 120, 240, 480, 960, 960 * 2, 960 * 4)
    var timeSize =
      //     Array(10, 20, 40, 80, 160, 320, 640)
      Array(25, 50, 100, 200, 400, 800, 1600) //, 1600*2)
    var spot = 100.0

    val tte = 0.25
    val r = 0.1
    val mu = 0.1
    val K1 = 90
    val K2 = 110
    val volMin = 0.15
    val volMax = 0.25

    val xMin = 0
    val xMax = 600
    var k = 0
    var previousPrice = Double.NaN
    var previousChange = Double.NaN
    val refPrice = 2.29769 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"

    println(" priceRef=" + refPrice);
    while (k < spaceSize.length) {
      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(k),
        timeSize(k),
        new MeshBoundaries(0, tte, xMin, xMax),
        K1)

      val spec = new UncertainBSM1FFDSpec(
        grid,
        volMin,
        volMax,
        mu,
        r,
        false)

      val payoff = new ButterflyFDPayoff(true, K1, K2, tte);
      val smoother = null //new SimpsonIntegralSmoother(Array(K1,(K1+K2)/2,K2));
      val solver = new ThomasTridiagonalSolver(payoff);
      val method = new ThetaUncertainBSM1FMethod(
        ThetaParabolic1DMethod.THETA_IMPLICIT)
      //        ThetaParabolic1DMethod.THETA_CRANK_NICOLSON);
      val pricer = new FDMSolver1D(spec, method, solver);
      pricer.smoother = smoother;

      pricer.solve(payoff);
      val price = pricer.price(spot);
      val change = price - previousPrice
      val ratio = previousChange / change
      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(k), timeSize(k), price, change, ratio, 0.0);
      println(lineStr);
      //      assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));
      previousPrice = price
      previousChange = change
      k += 1
    }

  }

  def createUniformArray(n: Int, xMin: Double, xMax: Double): Array[Double] = {
    var x = new Array[Double](n)
    var i = 0
    val dx = (xMax - xMin) / (n - 1)
    x(0) = xMin
    while (i < n - 1) {
      x(i + 1) = x(i) + dx
      i += 1
    }
    return x
  }

  def createConcentratedArray(n: Int, Smin: Double, Smax: Double, Sstar: Array[Double]): Array[Double] = {
    val alpha = 0.01 * (Smax - Smin)

    val x = createUniformArray(n, 0.0, 1.0)
    var S = new Array[Double](n)
    var k = 0
    var yMax = 0.0
    var A0 = 1.0
    var f = new JacobianConcentration(alpha, Sstar, A0)
    val solver = new RungeKuttaODESolver()

    solver.solve(f, x, Smin, S)
    var yMax0 = S(x.length - 1) - Smax
    var A1 = A0
    var yMax1 = yMax0
    //bracket: find yMax0 & yMax1, A0 & A1
    if (yMax0 > Smax) {
      do {
        A0 = A0 * 0.5
        var f = new JacobianConcentration(alpha, Sstar, A0)
        solver.solve(f, x, Smin, S)
        yMax0 = S(x.length - 1) - Smax

      } while (yMax0 > 0)
    } else {
      do {
        A1 = A1 * 2.0
        var f = new JacobianConcentration(alpha, Sstar, A1)
        solver.solve(f, x, Smin, S)
        yMax1 = S(x.length - 1) - Smax

      } while (yMax1 < 0)
    }
    //bisection
    do {
      var m = 0.5 * (A0 + A1) //bisection
      //      var m = (yMax1*A0-yMax0*A1)/(yMax1-yMax0) //regula falsi= worse because derivative is not significant
      f = new JacobianConcentration(alpha, Sstar, m)
      solver.solve(f, x, Smin, S)
      yMax = S(x.length - 1) - Smax
      if (yMax0 * yMax <= 0) {
        A1 = m
        yMax1 = yMax
      } else {
        A0 = m
        yMax0 = yMax
      }
      k += 1
    } while (k < 50 && Math.abs(yMax) > 1e-5)
    k = 1
    var i = 0
    while (false && k < (S.length - 1) && (i < Sstar.length)) {
      if (i < Sstar.length && S(k) > Sstar(i)) {
        val prevDiff = Math.abs(S(k - 1) - Sstar(i))
        val diff = Math.abs(S(k) - Sstar(i))
        if (prevDiff < diff) {
          S(k - 1) = Sstar(i)
          //          val dS = Math.min(S(k-1)-S(k-4),S(k+2)-S(k-1))
          //          S(k-3) = S(k-1)-dS*2.0/3.0
          //          S(k-4) = S(k-1)-dS
          //          S(k+1) = S(k-1)+dS*2.0/3.0
          //          S(k+2) = S(k-1)+dS
        } else {
          S(k) = Sstar(i)
          //         val dS = Math.min(S(k)-S(k-3),S(k+3)-S(k))
          //          S(k-2) = S(k)-dS*2.0/3.0
          //          S(k-3) = S(k)-dS
          //          S(k+2) = S(k)+dS*2.0/3.0
          //          S(k+3) = S(k)+dS
        }
        //        println("S changed, new S="+S(k-1)+" "+S(k))
        i += 1

      }
      k += 1
    }
    //    println("k=" + k + "S=" + Arrays.toString(S))
    return S
  }

  test("Butterfly Spread 1 Convergence Non Uniform") {
    var spaceSize = Array(60, 120, 240, 480, 960, 960 * 2, 960 * 4)
    var timeSize =
      //     Array(10, 20, 40, 80, 160, 320, 640)
      Array(25, 50, 100, 200, 400, 800, 1600) //, 1600*2)
    var spot = 100.0

    val tte = 0.25
    val r = 0.1
    val mu = 0.1
    val K1 = 90
    val K2 = 110
    val volMin = 0.15
    val volMax = 0.25

    val xMin = 0
    val xMax = 600
    var k = 0
    var previousPrice = Double.NaN
    var previousChange = Double.NaN
    val refPrice = 2.29769 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"

    println(" priceRef=" + refPrice);
    while (k < spaceSize.length) {
      val S = createConcentratedArray(spaceSize(k), xMin, xMax, Array(K1, (K1 + K2) * 0.5, K2))
      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(k),
        timeSize(k),
        new MeshBoundaries(0, tte, xMin, xMax),
        K1) {
        override def spaceVector = S;
      }

      val spec = new UncertainBSM1FFDSpec(
        grid,
        volMin,
        volMax,
        mu,
        r,
        false)

      val payoff = new ButterflyFDPayoff(true, K1, K2, tte);
      val smoother = new SimpsonIntegralSmoother(Array(K1,(K1+K2)/2,K2));
      val solver = new ThomasTridiagonalSolver(payoff);
      val method =
//        new TRBDF2UncertainBSM1FMethod(payoff)
              new ThetaUncertainBSM1FMethod(
              ThetaParabolic1DMethod.THETA_IMPLICIT)
      //              ThetaParabolic1DMethod.THETA_CRANK_NICOLSON);
      val pricer = new FDMSolver1D(spec, method, solver);
      pricer.smoother = smoother;

      pricer.solve(payoff);
      val price = pricer.price(spot);
      val change = price - previousPrice
      val ratio = previousChange / change
      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(k), timeSize(k), price, change, ratio, 0.0);
      println(lineStr);
      //      assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));
      previousPrice = price
      previousChange = change
      k += 1
    }      
  }

  test("Digital Convergence Non Uniform", SpecificTest) {
    var spaceSize = Array(60, 120, 240, 480, 960, 960 * 2, 960 * 4)
    var timeSize =
      //     Array(10, 20, 40, 80, 160, 320, 640)
      Array(25, 50, 100, 200, 400, 800, 1600) //, 1600*2)
    var spot = 100.0

    val tte = 0.25
    val r = 0.1
    val mu = 0.1
    val K = 100
    val volMin = 0.15
    val volMax = 0.25

    val xMin = 0
    val xMax = 300
    var k = 0
    var previousPrice = Double.NaN
    var previousChange = Double.NaN
    val refPrice = 0.4419641 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"

    println(" priceRef=" + refPrice);
    while (k < spaceSize.length) {
      val S = createConcentratedArray(spaceSize(k), xMin, xMax, Array(K))
      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(k),
        timeSize(k),
        new MeshBoundaries(0, tte, xMin, xMax),
        K) {
        override def spaceVector = S;
      }
//      grid.insertTime(Array(0.25*(1.0-2.0/10000),0.25*(1.0-1.0/10000)))
      var timeIte = grid.timeIterator
      val spec = new UncertainBSM1FFDSpec(
        grid,
        volMin,
        volMax,
        mu,
        r,
        false)

      val payoff = new DigitalFDPayoff(true, K, 1.0, tte);
      val smoother =
        //        new SplineSmoother(Array(K))
        new SimpsonIntegralSmoother(Array(K));
      val solver = new ThomasTridiagonalSolver(payoff);
      val method =
        new TRBDF2UncertainBSM1FMethod(payoff)
//              new ThetaUncertainBSM1FMethod(
//              ThetaParabolic1DMethod.THETA_IMPLICIT)
      //              ThetaParabolic1DMethod.THETA_CRANK_NICOLSON);
      val pricer = new FDMSolver1D(spec, method, solver);
      pricer.smoother = smoother;
      val timeIt = grid.timeIterator
      timeIt.next()
      timeIt.next()
      val gammaListener = new GammaListener(Array(timeIt.next()))
      pricer.listener = gammaListener;
      pricer.solve(payoff);
      val price = pricer.price(spot);
      val change = price - previousPrice
      val ratio = previousChange / change
      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(k), timeSize(k), price, change, ratio, 0.0);
      println(lineStr);
//      gammaListener.print()
      //      assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));
      previousPrice = price
      previousChange = change
      k += 1
    }

  }

  test("Butterfly Spread Price 2") {
    var spaceSize = 300 * 2
    var timeSize = 1000
    var spot = 100.0

    val tte = 0.5
    val r = 0.04
    val mu = 0.04
    val K1 = 95
    val K2 = 105
    val volMin = 0.30
    val volMax = 0.45

    val xMin = 0
    val xMax = 300

    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      new MeshBoundaries(0, tte, xMin, xMax),
      K1)

    val spec = new UncertainBSM1FFDSpec(
      grid,
      volMin,
      volMax,
      mu,
      r,
      false)

    val payoff = new ButterflyFDPayoff(true, K1, K2, tte);
    val smoother = null //new SimpsonIntegralSmoother(Array(K1,(K1+K2)/2,K2));
    val solver = new ThomasTridiagonalSolver(payoff);
    val method = new ThetaUncertainBSM1FMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);

    val pricer = new FDMSolver1D(spec, method, solver);
    pricer.smoother = smoother;

    pricer.solve(payoff);
    val price = pricer.price(spot);
    val refPrice = 0.125954 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"
    println(" priceRef=" + refPrice + " price=" + price);
    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));

  }

  test("Butterfly Spread 2 Convergence") {
    var spaceSize = Array(100, 200, 400, 800, 1600, 3200) //, 960*4)
    var timeSize =
      Array(20, 40, 80, 160, 320, 640)
    //      Array(100, 200, 400, 800, 1600, 3200) //, 1600*2)
    var spot = 100.0

    val tte = 0.5
    val r = 0.04
    val mu = 0.04
    val K1 = 95
    val K2 = 105
    val volMin = 0.30
    val volMax = 0.45

    val xMin = 0
    val xMax = 500
    var k = 0
    var previousPrice = Double.NaN
    var previousChange = Double.NaN
    val refPrice = 0.125954 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"

    println(" priceRef=" + refPrice);
    while (k < spaceSize.length) {
      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(k),
        timeSize(k),
        new MeshBoundaries(0, tte, xMin, xMax),
        K1)

      val spec = new UncertainBSM1FFDSpec(
        grid,
        volMin,
        volMax,
        mu,
        r,
        false)

      val payoff = new ButterflyFDPayoff(true, K1, K2, tte);
      val smoother = null //new SimpsonIntegralSmoother(Array(K1, (K1 + K2) / 2, K2));
      val solver = new ThomasTridiagonalSolver(payoff);
      val method = new ThetaUncertainBSM1FMethod(
        ThetaParabolic1DMethod.THETA_CRANK_NICOLSON //          ThetaParabolic1DMethod.THETA_IMPLICIT
        );

      val pricer = new FDMSolver1D(spec, method, solver);
      pricer.smoother = smoother;

      pricer.solve(payoff);
      val price = pricer.price(spot);
      val change = price - previousPrice
      val ratio = previousChange / change
      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(k), timeSize(k), price, change, ratio, 0.0);
      println(lineStr);
      //      assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));
      previousPrice = price
      previousChange = change
      k += 1
    }

  }
}
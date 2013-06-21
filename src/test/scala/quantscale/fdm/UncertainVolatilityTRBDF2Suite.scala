package test.quantscale.fdm
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.payoff.SimpsonIntegralSmoother
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.payoff.VanillaBermudanFDPayoff
import quantscale.fdm.UncertainBSM1FFDSpec
import quantscale.fdm.method.TRBDF2UncertainBSM1FMethod
import quantscale.fdm.payoff.ButterflyFDPayoff
import quantscale.fdm.listener.GammaListener
import quantscale.fdm._
import quantscale.fdm.mesh.UniformMesh2D
import quantscale.fdm.mesh.MeshBoundaries

@RunWith(classOf[JUnitRunner])
class UncertainVolatilityTRBDF2Suite extends FunSuite {
  test("Butterfly Spread Gamma 1") {
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
    val method = new TRBDF2UncertainBSM1FMethod(payoff);
    val pricer = new FDMSolver1D(spec, method, solver);
    pricer.smoother = smoother;
    val timeIterator = grid.timeIterator
    val gammaListener = new GammaListener(Array(0.0))
    pricer.listener = gammaListener;
    pricer.solve(payoff);
    val price = pricer.price(spot);
    val refPrice = 2.29769 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"
    println(" priceRef=" + refPrice + " price=" + price);

    gammaListener.print()
    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));

  }
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
    val method = new TRBDF2UncertainBSM1FMethod(payoff);
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
    val method = new TRBDF2UncertainBSM1FMethod(payoff);
    val pricer = new FDMSolver1D(spec, method, solver);
    pricer.smoother = smoother;

    pricer.solve(payoff);
    val price = pricer.price(spot);
    val refPrice = 2.29769 //from "Numerical Convergence Properties of Option Pricing PDEs with Uncertain Volatility"
    println(" priceRef=" + refPrice + " price=" + price);
    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - refPrice));

  }

  test("Butterfly Spread 1 Convergence") {
    var spaceSize = Array(60, 120, 240, 480, 960, 960 * 2, 960*4)
    var timeSize =
      //Array(10, 20, 40, 80, 160, 320, 640)
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
    val xMax = 300
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
      val smoother = null// new SimpsonIntegralSmoother(Array(K1,(K1+K2)/2,K2));
      val solver = new ThomasTridiagonalSolver(payoff);
      val method = new TRBDF2UncertainBSM1FMethod(payoff);
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
    val method = new TRBDF2UncertainBSM1FMethod(payoff);
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
      val smoother = null//new SimpsonIntegralSmoother(Array(K1, (K1 + K2) / 2, K2));
      val solver = new ThomasTridiagonalSolver(payoff);
      val method = new TRBDF2UncertainBSM1FMethod(payoff);
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
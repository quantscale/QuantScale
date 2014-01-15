package test.quantscale.fdm

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.Tag
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.LRE3Parabolic1DMethod
import quantscale.fdm.method.RannacherCentralBSM1FMethod
import quantscale.fdm.method.SmoothCrankNicolsonCentralBSM1FMethod
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.method.TRBDF2SingleCentralBSM1FMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.payoff._
import quantscale.fdm.transform.ExpTransformation
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.ConstantLogBSM1FFDSpec
import quantscale.fdm.ElliotOckendonTridiagonalSolver
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.SORTridiagonalSolver
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.TridiagonalSolver
import quantscale.math.CubicPP
import quantscale.math.CubicSpline
import org.scalatest.junit.JUnitRunner
import quantscale.fdm.mesh._
import quantscale.fdm.listener.Price2DListener
import quantscale.fdm.listener.Price2DListener
import java.io.PrintWriter
import java.io.File
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.TGAParabolic1DMethod
import quantscale.fdm.method.LMG3Parabolic1DMethod
import quantscale.fdm.method.TGASerialParabolic1DMethod

object SpecificTest extends Tag("SpecificTest") {

}

@RunWith(classOf[JUnitRunner])
class ConvergenceSuite extends FunSuite {

  class PriceLine(m: String, s: String, xLen: Int, tLen: Int, p: Double, e: Double, t: Double) {
    var method: String = m;
    var solver: String = s;
    var price: Double = p;
    var error: Double = e;
    var time: Double = t;
    var spaceSteps: Int = xLen;
    var timeSteps: Int = tLen;

    override def toString(): String = {
      var name = method;
      if (solver.length() > 0) {
        name = solver + "_" + name;
      }
      return "%d %s %.7f %.2e %.5f".format(timeSteps, name, price, error, time);
    }
  }

  def priceOSullivan(line: PriceLine, boundaries: MeshBoundaries, strike: Double, payoff: FDPayoff, spot: Double, vol: Double, mu: Double, r: Double, priceRef: Double) {
    val grid: UniformMesh2D = new UniformMesh2D(
      line.spaceSteps,
      line.timeSteps,
      boundaries,
      strike);
    val spec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r)
    val smoother = null;

    var solver: TridiagonalSolver = null;

    if (line.solver.equals("EO") || line.solver.equals("BS")) {
      solver = new ElliotOckendonTridiagonalSolver(payoff);
    } else if (line.solver.equals("SOR")) {
      solver = new SORTridiagonalSolver(payoff);
    } else {
      solver = new ThomasTridiagonalSolver(payoff);
    }
    var method: Parabolic1DMethod = null;
    if (line.method.equals("RAN")) {
      method = new RannacherCentralBSM1FMethod(Array(boundaries.lastTime));
    } else if (line.method.equals("C_RAN")) {
      val ranTimes = new ArrayBuffer[Double]();
      val it = grid.timeIterator;
      var i = 0;
      while (it.hasNext) {
        ranTimes += it.next();
        if (it.hasNext) it.next();
      }
      method = new RannacherCentralBSM1FMethod(ranTimes.reverse.toArray[Double]);
    } else if (line.method.equals("CN")) {
      method = new ThetaParabolic1DMethod();
    } else if (line.method.equals("TRBDF2")) {
      method = new TRBDF2Parabolic1DMethod(payoff);
    } else if (line.method.equals("TGA")) {
      method = new TGAParabolic1DMethod(payoff, 1.9999 - math.sqrt(2)) //1.25-0.5*math.sqrt(2));
    } else if (line.method.equals("TGAS")) {
      method = new TGASerialParabolic1DMethod(payoff, 2 - math.sqrt(2)) //1.25-0.5*math.sqrt(2));
    } else if (line.method.equals("EULER")) {
      method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
    } else if (line.method.equals("LMG2")) {
      method = new LMG2Parabolic1DMethod(payoff);
    } else if (line.method.equals("LMG3")) {
      method = new LMG3Parabolic1DMethod(payoff);
    } else if (line.method.equals("LRE3")) {
      method = new LRE3Parabolic1DMethod(payoff);
    } else if (line.method.equals("SCN")) {
      method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(boundaries.lastTime), payoff);
    }
    method.smearingReducer = null
    var pricer = new FDMSolver1D(spec, method, solver);

    pricer.smoother = smoother;
    var startTime = System.nanoTime();
    pricer.solve(payoff);

    if (line.solver.equals("RE")) {
      val reLine = new PriceLine(line.method, "", line.spaceSteps, line.timeSteps / 2, Double.NaN, Double.NaN, Double.NaN);
      priceOSullivan(reLine, boundaries, strike, payoff, spot, vol, mu, r, priceRef);
      line.time = (System.nanoTime() - startTime) * 1e-9;
      line.price = 2 * pricer.price(spot) - reLine.price;
      line.error = Math.abs(line.price - priceRef);
    } else {
      line.time = (System.nanoTime() - startTime) * 1e-9;
      line.price = pricer.price(spot);
      line.error = Math.abs(line.price - priceRef);
    }
    //
    //      var changeCN = Math.abs(priceCN - previousPriceCN);
    //      var ratioCN = previousChangeCN / changeCN;
    //      previousPriceCN = priceCN;
    //      previousChangeCN = changeCN;
    //      System.out.println(spaceSize(i) +
    //        " " + timeSize(i) +
    //        " " + priceCN +
    //        " " + changeCN +
    //        " " + ratioCN +
    //        " " + timeCN);

  }

  def priceOSullivanLog(line: PriceLine, boundaries: MeshBoundaries, strike: Double, payoff: FDPayoff, spot: Double, vol: Double, mu: Double, r: Double, priceRef: Double) {
    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime, Math.log(boundaries.bottomSpace / strike), Math.log(boundaries.topSpace / strike))
    val grid: UniformMesh2D = new UniformMesh2D(
      line.spaceSteps,
      line.timeSteps,
      logBounds,
      0);
    val spec = new ConstantLogBSM1FFDSpec(
      grid,
      vol,
      mu,
      r)
    grid.spaceTransform = new ExpTransformation(strike);
    val smoother = null; //new SimpsonIntegralSmoother(Array(strike));

    var solver: TridiagonalSolver = null;

    if (line.solver.equals("EO")) {
      solver = new ElliotOckendonTridiagonalSolver(payoff);
    } else if (line.solver.equals("SOR")) {
      solver = new SORTridiagonalSolver(payoff);
    } else {
      solver = new ThomasTridiagonalSolver(payoff);
    }
    var method: Parabolic1DMethod = null;
    if (line.method.equals("RAN")) {
      method = new RannacherCentralBSM1FMethod(Array(boundaries.lastTime));
    } else if (line.method.equals("CN")) {
      method = new ThetaParabolic1DMethod();
    } else if (line.method.equals("TRBDF2")) {
      method = new TRBDF2Parabolic1DMethod();
    } else if (line.method.equals("EULER")) {
      method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
    } else if (line.method.equals("LMG2")) {
      method = new LMG2Parabolic1DMethod(payoff);
    } else if (line.method.equals("LMG3")) {
      method = new LRE3Parabolic1DMethod(payoff);
    } else if (line.method.equals("SCN")) {
      method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(boundaries.lastTime), payoff);
    }
    method.smearingReducer = null
    var pricer = new FDMSolver1D(spec, method, solver);

    pricer.smoother = smoother;
    var startTime = System.nanoTime();
    pricer.solve(payoff);
    line.time = (System.nanoTime() - startTime) * 1e-9;
    line.price = pricer.price(0.0);
    line.error = Math.abs(line.price - priceRef);
    //
    //      var changeCN = Math.abs(priceCN - previousPriceCN);
    //      var ratioCN = previousChangeCN / changeCN;
    //      previousPriceCN = priceCN;
    //      previousChangeCN = changeCN;
    //      System.out.println(spaceSize(i) +
    //        " " + timeSize(i) +
    //        " " + priceCN +
    //        " " + changeCN +
    //        " " + ratioCN +
    //        " " + timeCN);

  }

  test("osullivan-american") {
    //don't smooth because this is not a real american price but the price given a specific discretization
    val reference = 6.0874933186;
    val tte = 1.0;
    val vol = 0.2;
    val r = 0.05;
    val mu = 0.05;
    val strike = 100;
    val spot = 100;
    val spaceSteps = 501;
    val boundaries = new MeshBoundaries(0, tte, 0, spaceSteps - 1);
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);
    val lines = new HashMap[String, ArrayBuffer[PriceLine]];
    for (loop <- 0 until 1) {
      var timeSteps = 20;

      System.out.println("loop " + loop);
      while (timeSteps <= 20060) {
        var line = new PriceLine("CN", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("RAN", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("SCN", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TRBDF2", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TGA", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TGAS", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LMG2", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LMG3", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LRE3", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("CN", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("RAN", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("SCN", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);

        line = new PriceLine("C_RAN", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TRBDF2", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TGA", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TGAS", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TRBDF2", "SOR", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LMG2", "SOR", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LMG3", "SOR", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);

        line = new PriceLine("LMG2", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LMG3", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("LRE3", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("TRBDF2", "RE", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        line = new PriceLine("EULER", "RE", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(line);
        timeSteps *= 2;
      }
    }
  }

  test("forsyth-log-american") {
    //smoothing impact convergence, grid size as well. see with 4 vs 5 which one is best.
    //without smoothing => better convergence ratio
    //with smoothing => lower ratio, but accuracy actually higher.
    var tte = 0.25;
    var r = 0.1;
    var mu = 0.1;
    var strike = 100.0;
    var spot = 100.0;
    val vol = 0.2;
    var priceRef = 14.67882;
    val timeSize = Array(25, 50, 100, 200, 400, 800);
    var spaceSize = Array(68, 135, 269, 537, 1073, 1073 * 2);
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);

    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new ThetaParabolic1DMethod(), "Thomas");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new ThetaParabolic1DMethod(), "EO");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new TRBDF2Parabolic1DMethod(), "Thomas");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new TRBDF2Parabolic1DMethod(), "EO");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new TRBDF2Parabolic1DMethod(), "SOR");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new RannacherCentralBSM1FMethod(Array(tte)), "Thomas");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new RannacherCentralBSM1FMethod(Array(tte)), "EO");
  }

  test("forsyth-log-european") {
    //smoothing impact convergence, grid size as well. see with 4 vs 5 which one is best.
    //without smoothing => better convergence ratio
    //with smoothing => lower ratio, but accuracy actually higher.
    var tte = 0.25;
    var r = 0.1;
    var mu = 0.1;
    var strike = 100.0;
    var spot = 100.0;
    val vol = 0.8;
    var priceRef = 14.45191;
    val timeSize = Array(25, 50, 100, 200, 400, 800);
    var spaceSize = Array(68, 135, 269, 537, 1073, 1073 * 2);
    val payoff = new VanillaFDPayoff(false, strike, tte);

    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new ThetaParabolic1DMethod(), "Thomas");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new TRBDF2Parabolic1DMethod(), "Thomas");
    testForsythLog(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, new RannacherCentralBSM1FMethod(Array(tte)), "Thomas");

  }

  def testForsythLog(tte: Double, r: Double, mu: Double, spot: Double, vol: Double, payoff: FDPayoff, strike: Double, priceRef: Double, timeSize: Array[Int], spaceSize: Array[Int], method: Parabolic1DMethod, solverType: String) {
    var previousPriceCN = Double.NaN;
    var previousChangeCN = Double.NaN;
    val zMin = mu * tte - 4 * vol * Math.sqrt(tte);
    val zMax = mu * tte + 4 * vol * Math.sqrt(tte);
    for (i <- 0 until timeSize.length) {

      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(i) * 2,
        timeSize(i),
        new MeshBoundaries(0, tte, zMin, zMax),
        0.0); //less will give higher european error (maybe because of order1 boundaries)
      grid.spaceTransform = new ExpTransformation(spot);
      val spec = new ConstantLogBSM1FFDSpec(
        grid,
        vol,
        mu,
        r)

      val smoother = null; //new SimpsonIntegralSmoother(Array(0.0));

      var solver: TridiagonalSolver = null;
      if (solverType.equals("EO")) {
        solver = new ElliotOckendonTridiagonalSolver(payoff);
      } else if (solverType.equals("SOR")) {
        solver = new SORTridiagonalSolver(payoff);
      } else {
        solver = new ThomasTridiagonalSolver();
      }
      var pricer = new FDMSolver1D(spec, method, solver);

      pricer.smoother = smoother;
      var startTime = System.nanoTime();
      pricer.solve(payoff);
      val timeCN = (System.nanoTime() - startTime) * 1e-9;
      val priceCN = pricer.price(0.0);
      val errorCN = Math.abs(priceCN - priceRef);

      var changeCN = Math.abs(priceCN - previousPriceCN);
      var ratioCN = previousChangeCN / changeCN;
      previousPriceCN = priceCN;
      previousChangeCN = changeCN;
      System.out.println(spaceSize(i) +
        " " + timeSize(i) +
        " " + priceCN +
        " " + changeCN +
        " " + ratioCN +
        " " + timeCN);

    }
  }

  //  test("forsyth-american") {
  test("forsyth-american", SpecificTest) {
    //smoothing impact convergence, grid size as well. see with 4 vs 5 which one is best.
    //without smoothing => better convergence ratio
    //with smoothing => lower ratio, but accuracy actually higher.
    var tte = 0.25;
    var r = 0.01;
    var mu = -0.1;
    var strike = 100.0;
    var spot = 100.0;
    val vol = 0.8;
    var priceRef = 14.67882;
    val timeSize = Array(25, 50, 100, 200, 400, 800);
    var spaceSize = Array(68, 135, 269, 537, 1073, 1073 * 2 + 1);
    val payoff = new VanillaAmericanFDPayoff(true, strike, 0, tte);

    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULER", "BS");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULERG", "BS");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2", "BS");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2G", "BS");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULER", "");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULERG", "");

    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "CN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "CN", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2single", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "SOR");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "RAN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "RAN", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "C_RAN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "C_RAN", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "SCN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "SCN", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "C_SCN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "C_SCN", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG3", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG3", "BS");

  }

  test("forsyth-european", SpecificTest) {
    //smoothing impact convergence, grid size as well. see with 4 vs 5 which one is best.
    var tte = 0.25;
    var r = 0.1;
    var mu = 0.1;
    var strike = 100.0;
    var spot = 100.0;
    val vol = 0.8;
    var priceRef = BlackScholesVanillaEuropean.priceEuropeanVanilla(false, strike, spot,
      BlackScholesVanillaEuropean.variance(vol, tte),
      BlackScholesVanillaEuropean.df(mu, tte),
      BlackScholesVanillaEuropean.df(r, tte));
    System.out.println("reference=" + priceRef);
    val timeSize = Array(25, 50, 100, 200, 400, 800);
    var spaceSize = Array(68, 135, 269, 539, 1073, 1072 * 2 + 1);
    val payoff = new VanillaFDPayoff(false, strike, tte);

    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULER", "");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULERG", "");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULER", "");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "EULERG", "");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2", "");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2G", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "CN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "BS");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "SOR");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "RAN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "SCN", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG2", "");
    //    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "LMG3", "");
  }

  test("forsythEuropeanCallBad") {
    var tte = 0.25;
    var r = 0.1;
    var mu = -0.1;
    var strike = 100.0;
    var spot = 100.0;
    val vol = 0.8;
    var priceRef = BlackScholesVanillaEuropean.priceEuropeanVanilla(true, strike, spot,
      BlackScholesVanillaEuropean.variance(vol, tte),
      BlackScholesVanillaEuropean.df(mu, tte),
      BlackScholesVanillaEuropean.df(r, tte));
    System.out.println("reference=" + priceRef);
    val timeSize = Array(25, 50, 100, 200, 400, 800);
    var spaceSize = Array(68, 135, 269, 539, 1073, 1072 * 2 + 1);
    val payoff = new VanillaFDPayoff(true, strike, tte);

    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "RAN", "");
  }

  def testForsyth(tte: Double, r: Double, mu: Double, spot: Double, vol: Double, payoff: FDPayoff, strike: Double, priceRef: Double, timeSize: Array[Int], spaceSize: Array[Int], methodType: String, solverType: String) {
    var previousPriceCN = Double.NaN;
    var previousChangeCN = Double.NaN;
    var scheme = methodType;
    if (solverType.length > 0) {
      scheme = solverType + "_" + scheme;
    }
    System.out.println(scheme);
    for (i <- 0 until timeSize.length) {
      val boundaries = new MeshBoundaries(0, tte, 0, 200); //350 orig

      //      val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3);

      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(i),
        timeSize(i),
        boundaries,
        strike); //less will give higher european error (maybe because of order1 boundaries)

      //      val x = grid.spaceVector;
      //      var j = 0;
      //      val index200 = ((spaceSize(i)-1)/1.240741).toInt;
      //      while (j<= index200) {
      //        x(j) = 200.0*j/index200;
      //        j+=1;
      //      }
      //      while (j< spaceSize(i)) {
      //        x(j) = 1000.0*j/(spaceSize(i)-1-index200)
      //        j+=1;
      //      }
      //      val dsSquare = (grid.spaceVector(grid.spaceVector.size-1)-grid.spaceVector(grid.spaceVector.size-2))*(grid.spaceVector(grid.spaceVector.size-1)-grid.spaceVector(grid.spaceVector.length-2))/(grid.spaceVector(grid.spaceVector.length-1)*grid.spaceVector(grid.spaceVector.length-1))
      //      println(dsSquare+" "+ tte/timeSize(i)+" " +(dsSquare-tte/timeSize(i)))
      val spec = new ConstantBSM1FFDSpec(
        grid,
        vol,
        mu,
        r)

      val smoother = null // new SimpsonIntegralSmoother(Array(strike));

      var solver: TridiagonalSolver = null;
      if (solverType.equals("EO") || solverType.equals("BS")) {
        solver = new ElliotOckendonTridiagonalSolver(payoff);
      } else if (solverType.equals("SOR")) {
        solver = new SORTridiagonalSolver(payoff);
      } else {
        solver = new ThomasTridiagonalSolver();
      }
      var method: Parabolic1DMethod = null;
      if (methodType.equals("RAN")) {
        method = new RannacherCentralBSM1FMethod(Array(boundaries.lastTime), null);
      } else if (methodType.equals("C_RAN")) {
        val ranTimes = new ArrayBuffer[Double]();
        val it = grid.timeIterator;
        var i = 0;
        while (it.hasNext) {
          ranTimes += it.next();
          if (it.hasNext) it.next();
        }
        method = new RannacherCentralBSM1FMethod(ranTimes.reverse.toArray[Double], null);
      } else if (methodType.equals("CN")) {
        method = new ThetaParabolic1DMethod();
      } else if (methodType.equals("EULER")) {
        method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
      } else if (methodType.equals("EULERG")) {
        method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
      } else if (methodType.equals("TRBDF2")) {
        method = new TRBDF2Parabolic1DMethod(payoff);
      } else if (methodType.equals("TRBDF2single")) {
        method = new TRBDF2SingleCentralBSM1FMethod()
      } else if (methodType.equals("LMG2G")) {
        method = new LMG2Parabolic1DMethod(payoff);
      } else if (methodType.equals("LMG2")) {
        method = new LMG2Parabolic1DMethod(payoff);
      } else if (methodType.equals("LMG3")) {
        method = new LRE3Parabolic1DMethod(payoff);
      } else if (methodType.equals("SCN")) {
        method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(), payoff);
      } else if (methodType.equals("C_SCN")) {
        val ranTimes = new ArrayBuffer[Double]();
        val it = grid.timeIterator;
        var i = 0;
        while (it.hasNext) {
          ranTimes += it.next();
          if (it.hasNext) it.next();
        }
        method = new SmoothCrankNicolsonCentralBSM1FMethod(ranTimes.toArray, null);
      }
      var pricer = new FDMSolver1D(spec, method, solver);
      if (timeSize(i) == 100) {
        pricer.listener = new Price2DListener()
      }

      pricer.smoother = smoother;
      var startTime = System.nanoTime();
      pricer.solve(payoff);
      if (pricer.listener != null) {
        val pw = new PrintWriter(new File("doc/lefloch_localextrapolation/price_european.txt"))
        pricer.listener.asInstanceOf[Price2DListener].print(pw)
        pw.close
      }
      val timeCN = (System.nanoTime() - startTime) * 1e-9;
      val priceCN = pricer.price(spot);
      val errorCN = math.abs(priceCN - priceRef);

      var changeCN = math.abs(priceCN - previousPriceCN);
      var ratioCN = previousChangeCN / changeCN;
      previousPriceCN = priceCN;
      previousChangeCN = changeCN;

      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(i), timeSize(i), priceCN, changeCN, ratioCN, timeCN);
      System.out.println(lineStr);

    }
  }

  class GilesPriceLine(m: String, s: String, xLen: Int, tLen: Int, p: Double, e: Double, t: Double, var maxError: Double, var maxGammaError: Double = Double.NaN) extends PriceLine(m: String, s: String, xLen: Int, tLen: Int, p: Double, e: Double, t: Double) {
    var ratio: Double = Double.NaN;
    var priceArray: Array[Double] = null;

    override def toString(): String = {
      var name = method;
      if (solver.length() > 0) {
        name = solver + "_" + name;
      }
      return "%d %s %.7f %.2e %.2e %.2e %.5f".format(timeSteps, name, price, maxGammaError, maxError, error, time);
    }

    def toLatexString(): String = {
      var name = method;
      if (solver.length() > 0) {
        name = solver + "\\_" + name;
      }
      return "%d & %d & %.2e & %.2e & %.2e & %.1f \\\\".format(spaceSteps, timeSteps, maxGammaError, maxError, time, ratio);
    }
  }

  def priceGilesLog(line: GilesPriceLine, boundaries: MeshBoundaries, strike: Double, payoff: FDPayoff, spot: Double, vol: Double, mu: Double, r: Double, priceRef: CubicPP) {
    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime, Math.log(boundaries.bottomSpace / spot), Math.log(boundaries.topSpace / spot))
    val grid: UniformMesh2D = new UniformMesh2D(
      line.spaceSteps,
      line.timeSteps,
      logBounds,
      Math.log(strike / spot));
    val spec = new ConstantLogBSM1FFDSpec(
      grid,
      vol,
      mu,
      r)
    grid.spaceTransform = new ExpTransformation(strike);
    val smoother = new SimpsonIntegralSmoother(Array(0.0));

    var solver: TridiagonalSolver = null;

    if (line.solver.equals("BS") || line.solver.equals("EO")) {
      solver = new ElliotOckendonTridiagonalSolver(payoff);
    } else if (line.solver.equals("SOR")) {
      solver = new SORTridiagonalSolver(payoff);
    } else {
      solver = new ThomasTridiagonalSolver();
    }
    var method: Parabolic1DMethod = null;
    if (line.method.equals("RAN")) {
      method = new RannacherCentralBSM1FMethod(Array(boundaries.lastTime));
    } else if (line.method.equals("CN")) {
      method = new ThetaParabolic1DMethod();
    } else if (line.method.equals("TRBDF2")) {
      method = new TRBDF2Parabolic1DMethod();
    } else if (line.method.equals("EULER")) {
      method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT);
    }
    var pricer = new FDMSolver1D(spec, method, solver);

    pricer.smoother = smoother;
    var startTime = System.nanoTime();
    pricer.solve(payoff);
    line.time = (System.nanoTime() - startTime) * 1e-9;
    line.price = pricer.price(0.0);
    val reference = priceRef.value(spot);
    line.error = Math.abs(line.price - reference);
    val min = 70;
    val max = 130;
    var i = 0;
    var priceArray = new Array[Double](grid.spaceVector.length);
    var x = grid.spaceTransform.transform(grid.spaceVector);
    while (i < grid.spaceVector.length) {
      priceArray(i) = pricer.price(grid.spaceVector(i));
      i += 1
    }
    line.priceArray = priceArray;
    if (line.solver.equals("RE")) {
      val reLine = new GilesPriceLine(line.method, "", line.spaceSteps, line.timeSteps / 2, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
      priceGilesLog(reLine, boundaries, strike, payoff, spot, vol, mu, r, priceRef);
      line.time = (System.nanoTime() - startTime) * 1e-9;
      val extrapolatedPrice = 2 * line.price - reLine.price;
      line.price = extrapolatedPrice;
      line.error = Math.abs(extrapolatedPrice - reference);
      i = 0;
      while (i < grid.spaceVector.length) {
        priceArray(i) = 2 * priceArray(i) - reLine.priceArray(i);
        i += 1
      }
      line.maxError = computeMaxPriceError(x, priceArray, priceRef, min, max);
      line.maxGammaError = computeMaxGammaError(x, priceArray, priceRef, min, max);
    } else {
      line.maxError = computeMaxPriceError(x, priceArray, priceRef, min, max);
      line.maxGammaError = computeMaxGammaError(x, priceArray, priceRef, min, max);
    }
  }

  test("giles-log-american") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;

    val lambda = Array(1.0, 2.0, 4.0, 8.0);
    val spaceSize = Array(40, 80, 160, 320, 640, 1280, 2560);

    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3);
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);
    val lines = new HashMap[String, ArrayBuffer[PriceLine]];
    //read File, interpolate with spline.
    val priceFunction = makePriceFromFile("giles_reference_10000.csv")
    //TODO compare max price in 60,140 and max gamma abs error
    val reference = priceFunction.value(spot);
    System.out.println("reference=" + reference);
    System.out.println("Spacesteps & Timesteps &  Max Gamma Error & Max Price Error & Time\\\\");
    System.out.println("Rannacher\\\\");
    for (k <- 0 until lambda.length) {
      System.out.println("lambda=" + lambda(k) + "\\\\");
      var previousLineRAN: GilesPriceLine = null;
      for (i <- 0 until spaceSize.length) {
        val spaceSteps = spaceSize(i);
        val timeSteps = Math.round(spaceSize(i) / lambda(k)).toInt
        var line: GilesPriceLine = null;
        line = new GilesPriceLine("RAN", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        priceGilesLog(line, boundaries, strike, payoff, spot, vol, mu, r, priceFunction);
        if (previousLineRAN != null) {
          line.ratio = previousLineRAN.maxError / line.maxError;
        }
        previousLineRAN = line;
        System.out.println(line.toLatexString());
      }
    }

    System.out.println("TR-BDF2\\\\");
    for (k <- 0 until lambda.length) {
      System.out.println("lambda=" + lambda(k) + "\\\\");
      var previousLineTRBDF2: GilesPriceLine = null;
      for (i <- 0 until spaceSize.length) {
        val spaceSteps = spaceSize(i);
        val timeSteps = Math.round(spaceSize(i) / lambda(k)).toInt
        var line: GilesPriceLine = null;

        line = new GilesPriceLine("TRBDF2", "BS", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        priceGilesLog(line, boundaries, strike, payoff, spot, vol, mu, r, priceFunction);
        if (previousLineTRBDF2 != null) {
          line.ratio = previousLineTRBDF2.maxError / line.maxError;
        }
        previousLineTRBDF2 = line;
        System.out.println(line.toLatexString());
      }
    }
    System.out.println("Rannacher with Richardson extrapolation\\\\");

    for (k <- 0 until lambda.length) {
      System.out.println("lambda=" + lambda(k) + "\\\\");
      var previousLineRAN: GilesPriceLine = null;
      for (i <- 0 until spaceSize.length) {
        val spaceSteps = spaceSize(i);
        val timeSteps = Math.round(spaceSize(i) / lambda(k)).toInt
        var line: GilesPriceLine = null;
        line = new GilesPriceLine("RAN", "RE", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        priceGilesLog(line, boundaries, strike, payoff, spot, vol, mu, r, priceFunction);
        if (previousLineRAN != null) {
          line.ratio = previousLineRAN.maxError / line.maxError;
        }
        previousLineRAN = line;
        System.out.println(line.toLatexString());
      }
    }
    System.out.println("TR-BDF2 with Richardson extrapolation\\\\");
    for (k <- 0 until lambda.length) {
      System.out.println("lambda=" + lambda(k) + "\\\\");
      var previousLineTRBDF2: GilesPriceLine = null;
      for (i <- 0 until spaceSize.length) {
        val spaceSteps = spaceSize(i);
        val timeSteps = Math.round(spaceSize(i) / lambda(k)).toInt
        var line: GilesPriceLine = null;

        line = new GilesPriceLine("TRBDF2", "RE", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        priceGilesLog(line, boundaries, strike, payoff, spot, vol, mu, r, priceFunction);
        if (previousLineTRBDF2 != null) {
          line.ratio = previousLineTRBDF2.maxError / line.maxError;
        }
        previousLineTRBDF2 = line;
        System.out.println(line.toLatexString());
      }
    }

    System.out.println("EULER with Richardson extrapolation\\\\");
    for (k <- 0 until lambda.length) {
      System.out.println("lambda=" + lambda(k) + "\\\\");
      var previousLineTRBDF2: GilesPriceLine = null;
      for (i <- 0 until spaceSize.length) {
        val spaceSteps = spaceSize(i);
        val timeSteps = Math.round(spaceSize(i) / lambda(k)).toInt
        var line: GilesPriceLine = null;

        line = new GilesPriceLine("EULER", "RE", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        priceGilesLog(line, boundaries, strike, payoff, spot, vol, mu, r, priceFunction);
        if (previousLineTRBDF2 != null) {
          line.ratio = previousLineTRBDF2.maxError / line.maxError;
        }
        previousLineTRBDF2 = line;
        System.out.println(line.toLatexString());
      }
    }
  }

  def computeMaxPriceError(x: Array[Double], price: Array[Double], reference: CubicPP, min: Double, max: Double): Double = {
    var i = 0;
    var maxError = 0.0;
    while (i < x.length && x(i) < min) i += 1;
    while (i < x.length && x(i) < max) {
      val error = Math.abs(price(i) - reference.value(x(i)))
      if (error > maxError) maxError = error;
      i += 1
    }
    return maxError;
  }

  def computeMaxGammaError(x: Array[Double], price: Array[Double], reference: CubicPP, min: Double, max: Double): Double = {
    var i = 0;
    var maxError = 0.0;
    var pp = CubicSpline.makeCubicSpline(x, price);
    while (i < x.length && x(i) < min) i += 1;
    while (i < x.length && x(i) < max) {
      val error = Math.abs(pp.secondDerivative(x(i)) - reference.secondDerivative(x(i)))
      if (error > maxError) maxError = error;
      i += 1
    }
    return maxError;
  }

  def makePriceFromFile(fileName: String): CubicPP = {
    val source = scala.io.Source.fromFile(fileName)
    val lines = source.getLines();
    val x = new ArrayBuffer[Double]();
    val y = new ArrayBuffer[Double]();
    lines.foreach(line => {
      val columns = line.split("\t");
      x += columns(0).toDouble;
      y += columns(1).toDouble;
    })
    source.close();
    return CubicSpline.makeCubicSpline(x.toArray[Double], y.toArray[Double]);
  }

  test("giles-american") {
    //TODO compute max error (in +-1STDEV).
    //TODO add relative convergence.
    val tte = 1.0;

    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    //read File, interpolate with spline.
    val priceFunction = makePriceFromFile("giles_reference_10000.csv")
    //TODO compare max price in 60,140 and max gamma abs error
    val reference = priceFunction.value(spot);
    System.out.println("reference=" + reference);
    val lambda = Array(1.0, 4.0, 8.0);
    val spaceSize = Array(40, 80, 160, 320, 640, 1280, 2560);

    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3);
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);
    val lines = new HashMap[String, ArrayBuffer[PriceLine]];

    System.out.println("ref error=" + (6.0903702250 + 7.57722903011449 - reference));

    for (k <- 0 until lambda.length) {
      for (i <- 0 until spaceSize.length) {
        val spaceSteps = spaceSize(i);
        val timeSteps = Math.round(spaceSize(i) / lambda(k)).toInt
        var line = new PriceLine("CN", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(lambda + " " + line);
        line = new PriceLine("RAN", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(lambda + " " + line);
        line = new PriceLine("TRBDF2", "", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(lambda + " " + line);
        line = new PriceLine("CN", "EO", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(lambda + " " + line);
        line = new PriceLine("RAN", "EO", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(lambda + " " + line);
        line = new PriceLine("TRBDF2", "EO", spaceSteps, timeSteps, Double.NaN, Double.NaN, Double.NaN);
        priceOSullivan(line, boundaries, strike, payoff, spot, vol, mu, r, reference);
        System.out.println(lambda + " " + line);
      }
    }
  }

  test("giles-european-smooth") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    val priceAnal = BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
      strike,
      spot,
      BlackScholesVanillaEuropean.variance(vol, tte),
      BlackScholesVanillaEuropean.df(mu, tte),
      BlackScholesVanillaEuropean.df(r, tte));

    val lambda = Array(1.0, 4.0);
    val spaceSize = Array(40, 80, 160, 320, 640, 1280, 2560);

    System.out.println("lambda spacesteps timesteps dx dt CNprice CNtime CN TRBDF2price TRBDF2time TRBDF2");

    for (k <- 0 until lambda.length) {
      for (i <- 0 until spaceSize.length) {
        val timeSize = Math.round(spaceSize(i) / lambda(k)).toInt
        val grid: UniformMesh2D = new UniformMesh2D(
          spaceSize(i),
          timeSize,
          MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3),
          strike);

        val spec = new ConstantBSM1FFDSpec(
          grid,
          vol,
          mu,
          r)

        val solver = new ThomasTridiagonalSolver();

        val payoff = new VanillaFDPayoff(false, strike, tte);
        val smoother = new SimpsonIntegralSmoother(Array(strike));

        var method: Parabolic1DMethod = new ThetaParabolic1DMethod();

        var pricer = new FDMSolver1D(spec, method, solver);

        pricer.smoother = smoother;
        var startTime = System.nanoTime();
        pricer.solve(payoff);
        val timeCN = (System.nanoTime() - startTime) * 1e-9;
        val priceCN = pricer.price(spot);
        val errorCN = Math.abs(priceCN - priceAnal);

        method = new TRBDF2Parabolic1DMethod(payoff);
        pricer = new FDMSolver1D(spec, method, solver);

        pricer.smoother = smoother;
        startTime = System.nanoTime();
        pricer.solve(payoff);
        val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
        val priceTRBDF2 = pricer.price(spot);
        val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);

        val dx = grid.spaceVector(1) - grid.spaceVector(0);
        val dt = 0;
        System.out.println(lambda(k) +
          " " + grid.spaceSize +
          " " + grid.timeIterator.length +
          " " + dx + " " + dt +
          " " + priceCN + " " + timeCN + " " + errorCN +
          " " + priceTRBDF2 + " " + timeTRBDF2 + " " + errorTRBDF2);
      }
    }

  }


  test("giles-european-projection") {
    val tte = 1.0;
    var spot = 100.0;
    var strike = 95.0;
    val vol = 0.3;
    val mu = 0.0;
    val r = 0.0;
    val priceAnal =
      (BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
        strike + 1e-5,
        spot,
        BlackScholesVanillaEuropean.variance(vol, tte),
        BlackScholesVanillaEuropean.df(mu, tte),
        BlackScholesVanillaEuropean.df(r, tte)) - BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
        strike - 1e-5,
        spot,
        BlackScholesVanillaEuropean.variance(vol, tte),
        BlackScholesVanillaEuropean.df(mu, tte),
        BlackScholesVanillaEuropean.df(r, tte))) / (2 * 1e-5);

    //      BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
    //      strike,
    //      spot,
    //      BlackScholesVanillaEuropean.variance(vol, tte),
    //      BlackScholesVanillaEuropean.df(mu, tte),
    //      BlackScholesVanillaEuropean.df(r, tte));

    val lambda = Array(1.0, 4.0);
    val spaceSize = Array(40, 80, 160, 320, 640, 1280, 2560);
    println("analytic price=" + priceAnal)
    System.out.println("lambda spacesteps timesteps dx dt CNprice CNtime CN TRBDF2price TRBDF2time TRBDF2");

    for (k <- 0 until lambda.length) {
      for (i <- 0 until spaceSize.length) {
        val timeSize = Math.round(spaceSize(i) / lambda(k)).toInt
        var grid: Mesh2D = new UniformMesh2D(
          spaceSize(i),
          timeSize,
          MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3),
          strike);

        var spec = new ConstantBSM1FFDSpec(
          grid,
          vol,
          mu,
          r)

        val solver = new ThomasTridiagonalSolver();
        var payoff = new DigitalFDPayoff(false, strike, 1.0, tte)
        //  new VanillaFDPayoff(false, strike, tte);
        val smoothers = Array[FDPayoffSmoother](new NoFDPayoffSmoother(), new SimpsonIntegralSmoother(Array(strike)), new ProjectionSmoother(Array(strike)))
        for (smoother <- smoothers) {

          var method: Parabolic1DMethod = new ThetaParabolic1DMethod();
          var pricer = new FDMSolver1D(spec, method, solver);

          pricer.smoother = smoother;
          var startTime = System.nanoTime();
          pricer.solve(payoff);
          var spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
          val timeCN = (System.nanoTime() - startTime) * 1e-9;
          val priceCN = spline.value(spot);
          val errorCN = Math.abs(priceCN - priceAnal);

          method = new TRBDF2Parabolic1DMethod(payoff);
          pricer = new FDMSolver1D(spec, method, solver);

          pricer.smoother = smoother;
          startTime = System.nanoTime();
          pricer.solve(payoff);
          spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
          val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
          val priceTRBDF2 = spline.value(spot);
          val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);

          val dx = grid.spaceVector(2) - grid.spaceVector(1);
          val dt = 0;
          println(smoother.getClass().getSimpleName() + " " + lambda(k) +
            " " + grid.spaceSize +
            " " + grid.timeIterator.length +
            " " + dx + " " + dt +
            " " + priceCN + " " + timeCN + " " + errorCN +
            " " + priceTRBDF2 + " " + timeTRBDF2 + " " + errorTRBDF2);
        }
        val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3)
        grid = new StaticAdaptiveMesh2D(new UniformMesh1D(spaceSize(i), boundaries.spaceBoundaries, strike, true),
          new UniformMesh1D(timeSize, boundaries.timeBoundaries))
        spec = new ConstantBSM1FFDSpec(
          grid,
          vol,
          mu,
          r)

        var method: Parabolic1DMethod = new ThetaParabolic1DMethod();
        var pricer = new FDMSolver1D(spec, method, solver);

        pricer.smoother = null;
        var startTime = System.nanoTime();
        pricer.solve(payoff);
        var spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
        val timeCN = (System.nanoTime() - startTime) * 1e-9;
        val priceCN = spline.value(spot);
        val errorCN = Math.abs(priceCN - priceAnal);

        method = new TRBDF2Parabolic1DMethod(payoff);
        pricer = new FDMSolver1D(spec, method, solver);

        pricer.smoother = null;
        startTime = System.nanoTime();
        pricer.solve(payoff);
        spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
        val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
        val priceTRBDF2 = spline.value(spot);
        val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);

        val dx = grid.spaceVector(2) - grid.spaceVector(1);
        val dt = 0;
        println("middle" + " " + lambda(k) +
          " " + grid.spaceSize +
          " " + grid.timeIterator.length +
          " " + dx + " " + dt +
          " " + priceCN + " " + timeCN + " " + errorCN +
          " " + priceTRBDF2 + " " + timeTRBDF2 + " " + errorTRBDF2);
      }
    }

  }

  test("giles-european-log-projection") {
    val tte = 1.0;
    var spot = 100.0;
    var strike = 95.0;
    val vol = 0.3;
    val mu = 0.0;
    val r = 0.0;
    val priceAnal =
    //      (BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
    //        strike+1e-5,
    //        spot,
    //        BlackScholesVanillaEuropean.variance(vol, tte),
    //        BlackScholesVanillaEuropean.df(mu, tte),
    //        BlackScholesVanillaEuropean.df(r, tte)) - BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
    //        strike-1e-5,
    //        spot,
    //        BlackScholesVanillaEuropean.variance(vol, tte),
    //        BlackScholesVanillaEuropean.df(mu, tte),
    //        BlackScholesVanillaEuropean.df(r, tte)))/(2*1e-5);

      BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
        strike,
        spot,
        BlackScholesVanillaEuropean.variance(vol, tte),
        BlackScholesVanillaEuropean.df(mu, tte),
        BlackScholesVanillaEuropean.df(r, tte));

    val lambda = Array(1.0, 4.0);
    val spaceSize = Array(40, 80, 160, 320, 640, 1280, 2560);
    println("analytic price=" + priceAnal)
    System.out.println("lambda spacesteps timesteps dx dt CNprice CNtime CN TRBDF2price TRBDF2time TRBDF2");

    for (k <- 0 until lambda.length) {
      for (i <- 0 until spaceSize.length) {
        val timeSize = Math.round(spaceSize(i) / lambda(k)).toInt
        val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3)
        val transform = new ExpTransformation(strike)
        var grid: Mesh2D = new StaticAdaptiveMesh2D(new TransformedUniformMesh1D(transform, spaceSize(i), new Mesh1DBoundaries(transform.inverseValue(boundaries.bottomSpace), transform.inverseValue(boundaries.topSpace)), new Point(transform.inverseValue(strike), false), true),
          new UniformMesh1D(timeSize, boundaries.timeBoundaries))

        var spec = new ConstantBSM1FFDSpec(
          grid,
          vol,
          mu,
          r)

        val solver = new ThomasTridiagonalSolver();
        var payoff = // new DigitalFDPayoff(false, strike, 1.0, tte)
          new VanillaFDPayoff(false, strike, tte);
        val smoothers = Array[FDPayoffSmoother](new NoFDPayoffSmoother(), new SimpsonIntegralSmoother(Array(strike)), new ProjectionSmoother(Array(strike)))
        for (smoother <- smoothers) {

          var method: Parabolic1DMethod = new ThetaParabolic1DMethod();
          var pricer = new FDMSolver1D(spec, method, solver);

          pricer.smoother = smoother;
          var startTime = System.nanoTime();
          pricer.solve(payoff);
          var spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
          val timeCN = (System.nanoTime() - startTime) * 1e-9;
          val priceCN = spline.value(spot);
          val errorCN = Math.abs(priceCN - priceAnal);

          method = new TRBDF2Parabolic1DMethod(payoff);
          pricer = new FDMSolver1D(spec, method, solver);

          pricer.smoother = smoother;
          startTime = System.nanoTime();
          pricer.solve(payoff);
          spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
          val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
          val priceTRBDF2 = spline.value(spot);
          val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);

          val dx = grid.spaceVector(2) - grid.spaceVector(1);
          val dt = 0;
          println(smoother.getClass().getSimpleName() + " " + lambda(k) +
            " " + grid.spaceSize +
            " " + grid.timeIterator.length +
            " " + dx + " " + dt +
            " " + priceCN + " " + timeCN + " " + errorCN +
            " " + priceTRBDF2 + " " + timeTRBDF2 + " " + errorTRBDF2);
        }
        grid = new StaticAdaptiveMesh2D(new TransformedUniformMesh1D(transform, spaceSize(i), new Mesh1DBoundaries(transform.inverseValue(boundaries.bottomSpace), transform.inverseValue(boundaries.topSpace)), new Point(transform.inverseValue(strike), true), true),
          new UniformMesh1D(timeSize, boundaries.timeBoundaries))
        spec = new ConstantBSM1FFDSpec(
          grid,
          vol,
          mu,
          r)

        var method: Parabolic1DMethod = new ThetaParabolic1DMethod();
        var pricer = new FDMSolver1D(spec, method, solver);

        pricer.smoother = null;
        var startTime = System.nanoTime();
        pricer.solve(payoff);
        var spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
        val timeCN = (System.nanoTime() - startTime) * 1e-9;
        val priceCN = spline.value(spot);
        val errorCN = Math.abs(priceCN - priceAnal);

        method = new TRBDF2Parabolic1DMethod(payoff);
        pricer = new FDMSolver1D(spec, method, solver);

        pricer.smoother = null;
        startTime = System.nanoTime();
        pricer.solve(payoff);
        spline = CubicSpline.makeCubicSpline(grid.spaceVector, pricer.price)
        val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
        val priceTRBDF2 = spline.value(spot);
        val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);

        val dx = grid.spaceVector(2) - grid.spaceVector(1);
        val dt = 0;
        println("middle" + " " + lambda(k) +
          " " + grid.spaceSize +
          " " + grid.timeIterator.length +
          " " + dx + " " + dt +
          " " + priceCN + " " + timeCN + " " + errorCN +
          " " + priceTRBDF2 + " " + timeTRBDF2 + " " + errorTRBDF2);
      }
    }

  }
  test("Digital Forsyth from uncertain") {
    var spaceSize = Array(60, 120, 240, 480, 960, 960 * 2, 960 * 4)
    var timeSize =
    //     Array(10, 20, 40, 80, 160, 320, 640)
      Array(25, 50, 100, 200, 400, 800, 1600) //, 1600*2)
    var spot = 100.0

    val tte = 0.25
    val r = 0.1
    val mu = 0.1
    val strike = 100
    val vol = 0.25

    val xMin = 0
    val xMax = 300
    System.out.println("spacesteps timesteps dx dt CNprice CNtime CN TRBDF2price TRBDF2time TRBDF2");

    var previousPriceCN, previousPriceTRBDF2, diffCN, diffTRBDF2: Double = 0.0
    for (i <- 0 until spaceSize.length) {
      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(i),
        timeSize(i),
        MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 5),
        strike);

      val spec = new ConstantBSM1FFDSpec(
        grid,
        vol,
        mu,
        r)

      val solver = new ThomasTridiagonalSolver();

      val payoff = new DigitalFDPayoff(true, strike, 1.0, tte);
      val smoother = new SplineSmoother(Array(strike)) //new SimpsonIntegralSmoother(Array(strike));

      var method: Parabolic1DMethod = new RannacherCentralBSM1FMethod(Array(tte))

      var pricer = new FDMSolver1D(spec, method, solver);
      val priceAnal = 0.4922364
      pricer.smoother = smoother;
      var startTime = System.nanoTime();
      pricer.solve(payoff);
      val timeCN = (System.nanoTime() - startTime) * 1e-9;
      val priceCN = pricer.price(spot);
      val errorCN = Math.abs(priceCN - priceAnal);
      val ratioCN = diffCN / (priceCN - previousPriceCN)
      diffCN = priceCN - previousPriceCN

      previousPriceCN = priceCN
      method = new TRBDF2Parabolic1DMethod(payoff);
      pricer = new FDMSolver1D(spec, method, solver);

      pricer.smoother = smoother;
      startTime = System.nanoTime();
      pricer.solve(payoff);
      val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
      val priceTRBDF2 = pricer.price(spot);
      val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);
      var ratioTRBDF2 = diffTRBDF2 / (priceTRBDF2 - previousPriceTRBDF2)
      diffTRBDF2 = priceTRBDF2 - previousPriceTRBDF2

      previousPriceTRBDF2 = priceTRBDF2
      val dx = grid.spaceVector(1) - grid.spaceVector(0);
      val dt = 0;
      System.out.println(
        grid.spaceSize +
          "\t" + grid.timeIterator.length +
          "\t" + dx +
          "\t" + priceCN + "\t" + timeCN + "\t" + ratioCN +
          "\t" + priceTRBDF2 + "\t" + timeTRBDF2 + "\t" + ratioTRBDF2);
    }
  }
  test("Digital Pooley") {
    val tte = 0.5;

    var spot = 40.0;
    var strike = 40.0;
    val vol = 0.3;
    val mu = 0.05;
    val r = 0.05;

    val spaceSize = Array(40, 80, 160, 320, 640, 1280)
    val timeSize = Array(25, 50, 100, 200, 400, 800)
    System.out.println("spacesteps timesteps dx dt CNprice CNtime CN TRBDF2price TRBDF2time TRBDF2");

    var previousPriceCN, previousPriceTRBDF2, diffCN, diffTRBDF2: Double = 0.0
    for (i <- 0 until spaceSize.length) {
      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(i),
        timeSize(i),
        new MeshBoundaries(0.0, tte, 0, 40 * 3.0),
        //MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 5),
        strike);

      val spec = new ConstantBSM1FFDSpec(
        grid,
        vol,
        mu,
        r)

      val solver = new ThomasTridiagonalSolver();

      val payoff = new DigitalFDPayoff(true, strike, 1.0, tte);
      val smoother = new SplineSmoother(Array(strike)) //new SimpsonIntegralSmoother(Array(strike));

      var method: Parabolic1DMethod = new ThetaParabolic1DMethod();

      var pricer = new FDMSolver1D(spec, method, solver);
      val priceAnal = 0.4922364
      pricer.smoother = smoother;
      var startTime = System.nanoTime();
      pricer.solve(payoff);
      val timeCN = (System.nanoTime() - startTime) * 1e-9;
      val priceCN = pricer.price(spot);
      val errorCN = Math.abs(priceCN - priceAnal);
      val ratioCN = diffCN / (priceCN - previousPriceCN)
      diffCN = priceCN - previousPriceCN

      previousPriceCN = priceCN
      method = new TRBDF2Parabolic1DMethod(payoff);
      pricer = new FDMSolver1D(spec, method, solver);

      pricer.smoother = smoother;
      startTime = System.nanoTime();
      pricer.solve(payoff);
      val timeTRBDF2 = (System.nanoTime() - startTime) * 1e-9;
      val priceTRBDF2 = pricer.price(spot);
      val errorTRBDF2 = Math.abs(priceTRBDF2 - priceAnal);
      var ratioTRBDF2 = diffTRBDF2 / (priceTRBDF2 - previousPriceTRBDF2)
      diffTRBDF2 = priceTRBDF2 - previousPriceTRBDF2

      previousPriceTRBDF2 = priceTRBDF2
      val dx = grid.spaceVector(1) - grid.spaceVector(0);
      val dt = 0;
      System.out.println(
        grid.spaceSize +
          " " + grid.timeIterator.length +
          " " + dx + " " + dt +
          " " + priceCN + " " + timeCN + " " + errorCN + " " + ratioCN +
          " " + priceTRBDF2 + " " + timeTRBDF2 + " " + errorTRBDF2 + " " + ratioTRBDF2);
    }

  }
}
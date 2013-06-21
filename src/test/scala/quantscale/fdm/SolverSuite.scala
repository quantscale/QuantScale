package test.quantscale.fdm
import java.util.Arrays
import scala.collection.mutable.ArrayBuffer
import org.scalatest.FunSuite
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.RannacherCentralBSM1FMethod
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.payoff.SimpsonIntegralSmoother
import quantscale.fdm.payoff.VanillaAmericanFDPayoff
import quantscale.fdm.payoff.VanillaFDPayoff
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.ElliotOckendonTridiagonalSolver
import quantscale.fdm.Epsilon
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.SORTridiagonalSolver
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.TridiagonalMatrix
import quantscale.fdm.TridiagonalSolver
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.fdm.mesh.UniformMesh2D
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod

class SolverSuite extends FunSuite {

  test("SOR-Convergence") {
    var tte = 1.0;
    var r = 0.01;
    var mu = 0.05;
    var strike = 100.0;
    var spot = 100.0;
    val vol = 0.4;
    var priceRef = 14.67882;
    val timeSize = Array( 200, 400, 800, 1600, 3200);
    var spaceSize = Array(9,17, 35, 69, 135, 269, 537, 1073, 1073 * 2 + 1);
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);

    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "BS");
    testForsyth(tte, r, mu, spot, vol, payoff, strike, priceRef, timeSize, spaceSize, "TRBDF2", "SOR");
    //TODO assertEquals(1e-5)
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
      val boundaries = //new GridBoundaries(0, tte, 0, 350);
            MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, 1.0, 3);

      val grid: UniformMesh2D = new UniformMesh2D(
        spaceSize(i),
        timeSize(i),
        boundaries,
        strike); //less will give higher european error (maybe because of order1 boundaries)

      val spec = new ConstantBSM1FFDSpec(
        grid,
        vol,
        mu,
        r)

      val smoother = new SimpsonIntegralSmoother(Array(strike));

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
        method = new RannacherCentralBSM1FMethod(Array(boundaries.lastTime));
      } else if (methodType.equals("C_RAN")) {
        val ranTimes = new ArrayBuffer[Double]();
        val it = grid.timeIterator;
        var i = 0;
        while (it.hasNext) {
          ranTimes += it.next();
          if (it.hasNext) it.next();
        }
        method = new RannacherCentralBSM1FMethod(ranTimes.reverse.toArray[Double]);
      } else if (methodType.equals("CN")) {
        method = new ThetaParabolic1DMethod();
      } else if (methodType.equals("TRBDF2")) {
        method = new TRBDF2Parabolic1DMethod(payoff);
      }
      var pricer = new FDMSolver1D(spec, method, solver);

      pricer.smoother = smoother;
      var startTime = System.nanoTime();
      pricer.solve(payoff);
      val timeCN = (System.nanoTime() - startTime) * 1e-9;
      val priceCN = pricer.price(spot);
      val errorCN = Math.abs(priceCN - priceRef);

      var changeCN = Math.abs(priceCN - previousPriceCN);
      var ratioCN = previousChangeCN / changeCN;
      previousPriceCN = priceCN;
      previousChangeCN = changeCN;

      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(i), timeSize(i), priceCN, changeCN, ratioCN, timeCN);
      System.out.println(lineStr);

    }
  }

  test("SOR-Solver") {
    val size = 10;
    val m = new TridiagonalMatrix(size);

    val x = new Array[Double](size);
    var i = 0;
    while (i < size) {
      m.lower(i) = Math.random;
      m.middle(i) = Math.random;
      m.upper(i) = m.lower(i) + 0.1; //Math.random;
      x(i) = Math.random;
      i += 1;
    }
    val y = new Array[Double](size);

    m.multiply(x, y);
    val dummyPayoff = new VanillaFDPayoff(false, 0.0, 0.0);
    val solver = new SORTridiagonalSolver(dummyPayoff);

    solver.init(size);
    val x0 = new Array[Double](size);
    solver.solve(m, y, x0);
    i = 0;
    System.out.println(Arrays.toString(x));
    System.out.println(Arrays.toString(x0));
    System.out.println(Epsilon.MACHINE_EPSILON_SQRT);
    while (i < size) {
      assert(Math.abs(x(i) - x0(i)) < Epsilon.MACHINE_EPSILON_SQRT, x(i) + " != " + x0(i) + " difference=" + Math.abs(x0(i) - x(i)));
      i += 1;
    }
  }
}
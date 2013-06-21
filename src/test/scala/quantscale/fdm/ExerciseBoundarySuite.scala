package test.quantscale.fdm
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.fdm.mesh.UniformMesh2D
import quantscale.fdm.ConstantLogBSM1FFDSpec
import java.io.PrintWriter
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.payoff.VanillaAmericanFDPayoff
import quantscale.fdm.ElliotOckendonTridiagonalSolver
import java.io.File
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.fdm.transform.ExpTransformation
import quantscale.fdm.payoff.BoundaryListener
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.method.ThetaParabolic1DMethod

@RunWith(classOf[JUnitRunner])
class ExerciseBoundarySuite extends FunSuite {
  test("american-atm-put") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 165.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    val spaceSize = 50;
    val timeSize = spaceSize;
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte)
    payoff.exerciseBoundary = new BoundaryListener()
    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3);

    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime, Math.log(boundaries.bottomSpace / strike), Math.log(boundaries.topSpace / strike))
    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      logBounds,
      0.0);

    val spec = new ConstantLogBSM1FFDSpec(
      grid,
      vol,
      mu,
      r);
    grid.spaceTransform = new ExpTransformation(strike);

    val f = new File("doc/lefloch_localextrapolation/exercise_american_atm.txt")
    val pw = new PrintWriter(f);
    pw.println("t\tS");

    val solver = new ElliotOckendonTridiagonalSolver(payoff);
    var method: Parabolic1DMethod = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
    var pricer = new FDMSolver1D(spec, method, solver);
    pricer.solve(payoff);
    for (i <- 0 until payoff.exerciseBoundary.t.size) {
      var line = payoff.exerciseBoundary.t(i) + "\t" + payoff.exerciseBoundary.boundary(i)
      println(line)
      pw.println(line)
    }
    pw.close()
  }
}
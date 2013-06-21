package test.quantscale.fdm

import java.io.File
import java.io.PrintWriter
import java.util.Arrays
import scala.collection.mutable.ArrayBuffer
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.analytic.BSMMeasure
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.analytic.CoxIngersollRossBond
import quantscale.analytic.GammaGreek
import quantscale.analytic.Price
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.ConstantCIRSpec
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.mesh.Mesh1DBoundaries
import quantscale.fdm.mesh.Mesh2D
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.fdm.mesh.Point
import quantscale.fdm.mesh.SMPoint
import quantscale.fdm.mesh.SinhMesh1D
import quantscale.fdm.mesh.SplineMesh1D
import quantscale.fdm.mesh.StaticAdaptiveMesh2D
import quantscale.fdm.mesh.UniformMesh1D
import quantscale.fdm.method.BackwardLinearOrder1Parabolic1DBoundaryFactory
import quantscale.fdm.method.BackwardPartialOrder2Parabolic1DBoundaryFactory
import quantscale.fdm.method.ForwardLinearOrder1Parabolic1DBoundaryFactory
import quantscale.fdm.method.ForwardPartialOrder2Parabolic1DBoundaryFactory
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.LRE3Parabolic1DMethod
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.RannacherCentralBSM1FMethod
import quantscale.fdm.method.SmoothCrankNicolsonCentralBSM1FMethod
import quantscale.fdm.method.TRBDF2CentralParabolic1DMethod
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.payoff.BondFDPayoff
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.payoff.VanillaAmericanFDPayoff
import quantscale.fdm.payoff.VanillaBondOptionFDPayoff
import quantscale.fdm.payoff.VanillaFDPayoff
import quantscale.fdm.method.ThetaParabolic1DMethod

@RunWith(classOf[JUnitRunner])
class BoundarySuite extends FunSuite {

  def createMethod(schemeName: String, payoff: FDPayoff, grid: Mesh2D): Parabolic1DMethod = {
    var method: Parabolic1DMethod = null
    val tte = grid.timeIterator.next
    schemeName match {
      case "LMG2-C1"      => method = new LMG2Parabolic1DMethod(payoff)
      case "LMG2" | "LMG2-T1" => 
        method = new LMG2Parabolic1DMethod(payoff)
        method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory
      case "LMG2-T2" =>
        method = new LMG2Parabolic1DMethod(payoff)
        method.lowerBoundary = ForwardPartialOrder2Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
      case "RAN"       => method = new RannacherCentralBSM1FMethod(Array(tte), payoff)
      case "CN"        => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
      case "TRBDF2-C1" => 
        method = new TRBDF2CentralParabolic1DMethod(payoff)
        method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory
      case "TRBDF2" | "TRBDF2-T1" =>
        method = new TRBDF2Parabolic1DMethod(payoff)
        method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory
      case "TRBDF2-T2" =>
        method = new TRBDF2Parabolic1DMethod(payoff)
        method.lowerBoundary = ForwardPartialOrder2Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
      case "EULER-C1" => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
      case "EULER" | "EULER-T1" =>
        //          method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory
      case "EULER-T2" =>
        method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        method.lowerBoundary = ForwardPartialOrder2Parabolic1DBoundaryFactory
        method.upperBoundary = BackwardPartialOrder2Parabolic1DBoundaryFactory
      case "LMG3" => method = new LRE3Parabolic1DMethod(payoff)
      case "SCN"  => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(tte), payoff)
      case "CRAN" =>
        val ranTimes = new ArrayBuffer[Double]();
        val it = grid.timeIterator;
        var i = 0;
        while (it.hasNext) {
          if (it.hasNext) it.next();
        }
        method = new RannacherCentralBSM1FMethod(ranTimes.toArray, payoff)
      case "CSCN" =>
        val ranTimes = new ArrayBuffer[Double]();
        val it = grid.timeIterator;
        var i = 0;
        while (it.hasNext) {
          ranTimes += it.next();
          if (it.hasNext) it.next();
        }
        method = new SmoothCrankNicolsonCentralBSM1FMethod(ranTimes.toArray, payoff)
    }
    return method
  }

  test("american-put") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.04;
    val r = 0.04;
    val spaceSize = 500;
    val timeSize = 5000;
    val isCall = false
    val payoff = new VanillaAmericanFDPayoff(isCall, strike, 0.0, tte);
    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 2);

    val spaceBoundaries = new Mesh1DBoundaries(boundaries.bottomSpace, boundaries.topSpace)
    val timeBoundaries = new Mesh1DBoundaries(0, tte)

    val spaceMesh = new SplineMesh1D(spaceSize, spaceBoundaries, Array(new SMPoint(strike, false, 0, 0.1))) //new LogUniformMesh1D(spaceSize, spaceBoundaries, strike)
    val farAwayPoint = spot * math.exp(4 * vol * math.sqrt(tte))
    val intermediatePoint = 0.5 * (farAwayPoint + spaceMesh.x(spaceMesh.size - 1))
    spaceMesh.insertPoints(Array(intermediatePoint, farAwayPoint))
    val grid: Mesh2D = new StaticAdaptiveMesh2D(
      spaceMesh,
      new UniformMesh1D(timeSize, timeBoundaries))

    val spec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r);
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("EULER", "EULER1", "EULER2") //,"TRBDF2")
    val f = new File("doc/lefloch_localextrapolation/boundary_american_put.txt")
    val pw = new PrintWriter(f);
    pw.println("S\tScheme\tPrice\tAPrice\tPriceError\tGamma\tAGamma\tGammaError");

    for (schemeName <- methodNames) {
      val solver = new ThomasTridiagonalSolver(payoff);
      var method = createMethod(schemeName, payoff, grid)
      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      var i = 0;
      var x = grid.spaceTransform.transform(grid.spaceVector);
      while (i < grid.spaceVector.length) {
        priceArray(i) = pricer.price(grid.spaceVector(i));
        i += 1
      }
      i = 1;
      val variance = BlackScholesVanillaEuropean.variance(vol, tte)
      val driftDf = BlackScholesVanillaEuropean.df(mu, tte)
      val discountDf = BlackScholesVanillaEuropean.df(r, tte)
      val priceAn = new Price()
      val gammaAn = new GammaGreek()
      val measures: Array[BSMMeasure] = Array(priceAn, gammaAn)
      while (i < grid.spaceVector.length - 1) {
        var gamma = 0.
        var delta = 0.
        var sdelta = 0.
        var sgamma = 0.
        delta = (priceArray(i + 1) - priceArray(i - 1)) / (x(i + 1) - x(i - 1))
        gamma = 2 * ((x(i + 1) - x(i)) * priceArray(i - 1) - (x(i + 1) - x(i - 1)) * priceArray(i) + (x(i) - x(i - 1)) * priceArray(i + 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)) * (x(i + 1) - x(i - 1)))

        BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, x(i), variance, driftDf, discountDf, measures)

        pw.println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + priceAn.value + "\t" + math.abs(priceArray(i) - priceAn.value) + "\t" + gamma + "\t" + gammaAn.value + "\t" + math.abs(gamma - gammaAn.value))
        println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + priceAn.value + "\t" + math.abs(priceArray(i) - priceAn.value) + "\t" + gamma + "\t" + gammaAn.value + "\t" + math.abs(gamma - gammaAn.value))
        i += 1
      }
    }
    pw.close()
  }

  test("european-put") {
    val tte = 0.5;

    var spot = 100.0;
    var strike = 120.0;
    val vol = 0.2;
    val mu = 0.15;
    val r = 0.05;
    val spaceSize = 200;
    val timeSize = 5000;
    val isCall = false
    val payoff = new VanillaFDPayoff(isCall, strike, tte);
    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 5);

    val spaceBoundaries = new Mesh1DBoundaries(boundaries.bottomSpace, boundaries.topSpace)
    val timeBoundaries = new Mesh1DBoundaries(0, tte)

    val spaceMesh = new SinhMesh1D(spaceSize, spaceBoundaries, false, Array(new Point(strike, false)), strike, SinhMesh1D.ALPHA_HIGHLY_NONUNIFORM)
    //new SplineMesh1D(spaceSize, spaceBoundaries, Array(new SMPoint(strike,false,0.2,1e-5)))//new LogUniformMesh1D(spaceSize, spaceBoundaries, strike)
    //    val farAwayPoint =spot * math.exp(4 * vol * math.sqrt(tte))
    //    val intermediatePoint = 0.5 * (farAwayPoint + spaceMesh.x(spaceMesh.size - 1))
    //    spaceMesh.insertPoints(Array(intermediatePoint,farAwayPoint))
    val grid: Mesh2D = new StaticAdaptiveMesh2D(
      spaceMesh,
      new UniformMesh1D(timeSize, timeBoundaries))

    val spec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r);
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("EULER", "EULER1", "EULER2") //,"TRBDF2")
    val f = new File("doc/lefloch_localextrapolation/boundary_european_put.txt")
    val pw = new PrintWriter(f);
    pw.println("S\tScheme\tPrice\tAPrice\tPriceError\tGamma\tAGamma\tGammaError");

    for (schemeName <- methodNames) {
      val solver = new ThomasTridiagonalSolver(payoff);
      var method = createMethod(schemeName, payoff, grid)
      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      var i = 0;
      var x = grid.spaceTransform.transform(grid.spaceVector);
      while (i < grid.spaceVector.length) {
        priceArray(i) = pricer.price(grid.spaceVector(i));
        i += 1
      }
      i = 1;
      val variance = BlackScholesVanillaEuropean.variance(vol, tte)
      val driftDf = BlackScholesVanillaEuropean.df(mu, tte)
      val discountDf = BlackScholesVanillaEuropean.df(r, tte)
      val priceAn = new Price()
      val gammaAn = new GammaGreek()
      val measures: Array[BSMMeasure] = Array(priceAn, gammaAn)
      while (i < grid.spaceVector.length - 1) {
        var gamma = 0.
        var delta = 0.
        var sdelta = 0.
        var sgamma = 0.
        delta = (priceArray(i + 1) - priceArray(i - 1)) / (x(i + 1) - x(i - 1))
        gamma = 2 * ((x(i + 1) - x(i)) * priceArray(i - 1) - (x(i + 1) - x(i - 1)) * priceArray(i) + (x(i) - x(i - 1)) * priceArray(i + 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)) * (x(i + 1) - x(i - 1)))

        BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, x(i), variance, driftDf, discountDf, measures)

        pw.println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + priceAn.value + "\t" + math.abs(priceArray(i) - priceAn.value) + "\t" + gamma + "\t" + gammaAn.value + "\t" + math.abs(gamma - gammaAn.value))
        println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + priceAn.value + "\t" + math.abs(priceArray(i) - priceAn.value) + "\t" + gamma + "\t" + gammaAn.value + "\t" + math.abs(gamma - gammaAn.value))
        i += 1
      }
    }
    pw.close()
  }

  test("cir") {
    val tte = 10.0;
    val aCIR = 0.04
    val bCIR = 0.8
    val sigmaCIR = 0.275
    val rStar = 0.05
    val rMax = .2
    val spaceSize = 200
    val timeSize = 5000

    val payoff = new BondFDPayoff(tte);

    val logBounds = new MeshBoundaries(0, tte, 0, rMax)
    val spaceTransform = new UniformMesh1D(spaceSize, logBounds.spaceBoundaries)
    spaceTransform.insertPoints(Array(0.5, 1.0))
    val timeTransform = new UniformMesh1D(timeSize, logBounds.timeBoundaries)

    val grid: Mesh2D = new StaticAdaptiveMesh2D(spaceTransform, timeTransform)

    val spec = new ConstantCIRSpec(grid, aCIR, bCIR, sigmaCIR)
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("EULER-C1", "EULER-T1", "EULER-T2") //"CN", "TRBDF2", "LMG2", "LMG2G", "EULER", "EULERG")
    val f = new File("doc/lefloch_localextrapolation/cir_bond.txt")
    val pw = new PrintWriter(f);
    pw.println("Scheme Spot Price AnPrice Error")
    for (schemeName <- methodNames) {
      val solver = new ThomasTridiagonalSolver(payoff);
      var method: Parabolic1DMethod = createMethod(schemeName, payoff, grid)

      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      priceArray = pricer.price
      for (i <- 0 until priceArray.length) {
        val spot = grid.spaceVector(i)
        val priceAn = CoxIngersollRossBond.price(bCIR, aCIR / bCIR, sigmaCIR, spot, tte)
        val error = priceArray(i) - priceAn
        println(schemeName + " " + spot + " " + priceArray(i) + " " + priceAn + " " + error)
        pw.println(schemeName + " " + spot + " " + priceArray(i) + " " + priceAn + " " + error)
      }
    }
    pw.close()
  }

  test("cir-call") {
    val S = 10.0
    val T = 5.0
    val isCall = true
    val strike = 0.7

    val aCIR = 0.04
    val bCIR = 0.8
    val sigmaCIR = 0.275
    val rStar = 0.05
    val rMax = .4
    val spaceSize = 201
    val timeSize = 1000

    val payoff = new VanillaBondOptionFDPayoff(isCall, strike, T, S);

    val logBounds = new MeshBoundaries(0, S, 0, rMax)
    val spaceTransform = new UniformMesh1D(spaceSize, logBounds.spaceBoundaries)
    spaceTransform.insertPoints(Array(0.5, 1.0))
    val timeTransform = new UniformMesh1D(timeSize, logBounds.timeBoundaries)
    timeTransform.insertPoints(Array(T))
    val grid: Mesh2D = new StaticAdaptiveMesh2D(spaceTransform, timeTransform)

    val spec = new ConstantCIRSpec(grid, aCIR, bCIR, sigmaCIR)
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("EULER-C1", "EULER-T1", "EULER-T2", "CN", "TRBDF2-T2", "LMG2-T2") //"CN", "TRBDF2", "LMG2", "LMG2G", "EULER", "EULERG")
    val f = new File("doc/lefloch_localextrapolation/cir_bondoption.txt")
    val pw = new PrintWriter(f);
    pw.println("Scheme Spot Price AnPrice Error")
    for (schemeName <- methodNames) {
      val solver = new ThomasTridiagonalSolver(payoff);
      var method: Parabolic1DMethod = createMethod(schemeName, payoff, grid)

      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      priceArray = pricer.price
      val startTime = System.nanoTime()
      for (i <- 0 until priceArray.length) {
        val spot = grid.spaceVector(i)
        val priceAn = CoxIngersollRossBond.priceVanillaOption(isCall, strike, bCIR, aCIR / bCIR, sigmaCIR, spot, T, S)
        val error = priceArray(i) - priceAn
        println(schemeName + " " + spot + " " + priceArray(i) + " " + priceAn + " " + error)
        pw.println(schemeName + " " + spot + " " + priceArray(i) + " " + priceAn + " " + error)
      }
      val endTime = System.nanoTime()
      //      println((endTime-startTime)*1e-9)
    }
    pw.close()
  }

}
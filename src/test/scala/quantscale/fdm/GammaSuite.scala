package test.quantscale.fdm

import java.io.File
import java.io.PrintWriter
import java.util.Arrays
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.scalatest.FunSuite
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.fdm.mesh.UniformMesh2D
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.LRE3Parabolic1DMethod
import quantscale.fdm.method.RannacherCentralBSM1FMethod
import quantscale.fdm.method.SmoothCrankNicolsonCentralBSM1FMethod
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.payoff.SimpsonIntegralSmoother
import quantscale.fdm.payoff.VanillaAmericanFDPayoff
import quantscale.fdm.payoff.VanillaBermudanFDPayoff
import quantscale.fdm.payoff.VanillaFDPayoff
import quantscale.fdm.payoff.VanillaFDPayoff
import quantscale.fdm.transform.ExpTransformation
import quantscale.fdm._
import quantscale.math.CubicSpline
import quantscale.math.NotAKnot
import quantscale.math.NotAKnot
import quantscale.analytic.GammaGreek
import quantscale.analytic.BSMMeasure
import quantscale.analytic.Price
import scala.collection.mutable.ArrayBuffer
import quantscale.fdm.mesh.Mesh2D
import quantscale.fdm.mesh.StaticAdaptiveMesh2D
import quantscale.fdm.mesh.LogUniformMesh1D
import quantscale.fdm.mesh.Mesh1DBoundaries
import quantscale.fdm.mesh.UniformMesh1D
import quantscale.fdm.listener.Price2DListener
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod

@RunWith(classOf[JUnitRunner])
class GammaSuite extends FunSuite {

  test("gamma-itm-log-american") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 165 //165.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    val spaceSize = 500;
    val timeSize = spaceSize / 32;
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);
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
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("CN", "TRBDF2", "RAN", "CRAN", "SCN", "CSCN", "LMG2", "LMG3")
    val f = new File("doc/lefloch_localextrapolation/gamma_american_itm.txt")
    val pw = new PrintWriter(f);
    pw.println("S\tScheme\tPrice\tGamma\tSGamma");

    for (schemeName <- methodNames) {
      val solver = new ElliotOckendonTridiagonalSolver(payoff);
      var method: Parabolic1DMethod = null
      schemeName match {
        case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
        case "RAN"    => method = new RannacherCentralBSM1FMethod(Array(tte), payoff)
        case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
        case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
        case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
        case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(tte), payoff)
        case "CRAN" =>
          val ranTimes = new ArrayBuffer[Double]();
          val it = grid.timeIterator;
          var i = 0;
          while (it.hasNext) {
            ranTimes += it.next();
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
      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      var i = 0;
      var x = grid.spaceTransform.transform(grid.spaceVector);
      priceArray = pricer.price
      //      while (i < grid.spaceVector.length) {
      //        priceArray(i) = pricer.price(grid.spaceVector(i));
      //        i += 1
      //      }
      var pp = CubicSpline.makeCubicSpline(x, priceArray);
      i = 0;
      while (i < grid.spaceVector.length) {
        var sgamma = pp.secondDerivative(x(i))
        var gamma = 0.
        if (i == 0 || i == grid.spaceVector.length - 1) {
        } else {
          //          delta = (priceArray(i + 1) - priceArray(i - 1)) / (x(i + 1) - x(i - 1))
          gamma = 2 * ((x(i + 1) - x(i)) * priceArray(i - 1) - (x(i + 1) - x(i - 1)) * priceArray(i) + (x(i) - x(i - 1)) * priceArray(i + 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)) * (x(i + 1) - x(i - 1)))
        }
        println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma)
        pw.println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma)
        i += 1
      }
    }
    pw.close()
  }

  test("gamma-itm-american") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 165 //165.0;
    val vol = 0.4;
    val mu = 0.02; //0.1 and -0.1 should pass
    val r = 0.0;
    val spaceSize = 500;
    val timeSize = spaceSize / 64;
    val payoff = new VanillaAmericanFDPayoff(false, strike, 0, tte);
    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, strike, vol, BlackScholesVanillaEuropean.df(mu, tte), 2);

    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime, boundaries.bottomSpace, boundaries.topSpace)
    val spaceBoundaries = new Mesh1DBoundaries(boundaries.bottomSpace, boundaries.topSpace)
    val timeBoundaries = new Mesh1DBoundaries(0, tte)

    val grid: Mesh2D = new StaticAdaptiveMesh2D(
      new LogUniformMesh1D(spaceSize, spaceBoundaries, strike),
      new UniformMesh1D(timeSize, timeBoundaries))

    //     new UniformMesh2D(
    //      spaceSize,
    //      timeSize,
    //      logBounds,
    //      strike);

    val spec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r);
    println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("EULER", "EULERG")
    //    Array("CN", "CRAN", "TRBDF2","EULER", "LMG2", "LMG3")
    val f = new File("doc/lefloch_localextrapolation/gamma_american_itm.txt")
    val pw = new PrintWriter(f)
    pw.println("S\tScheme\tPrice\tGamma\tSGamma");

    for (schemeName <- methodNames) {
      val solver = new ElliotOckendonTridiagonalSolver(payoff);
      var method: Parabolic1DMethod = null
      schemeName match {
        case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
        case "RAN"    => method = new RannacherCentralBSM1FMethod(Array(tte), payoff)
        case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
        case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
        case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "EULERG" => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
        case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(tte), payoff)
        case "CRAN" =>
          val ranTimes = new ArrayBuffer[Double]();
          val it = grid.timeIterator;
          var i = 0;
          while (it.hasNext) {
            ranTimes += it.next();
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
      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.listener = new Price2DListener
      pricer.solve(payoff);
      if (pricer.listener != null) {
        val pw = new PrintWriter(new File("doc/lefloch_localextrapolation/price_european.txt"))
        pricer.listener.asInstanceOf[Price2DListener].print(pw)
        pw.close
      }
      var priceArray = new Array[Double](grid.spaceVector.length);
      var i = 0;
      var x = grid.spaceTransform.transform(grid.spaceVector);
      priceArray = pricer.price
      //      while (i < grid.spaceVector.length) {
      //        priceArray(i) = pricer.price(grid.spaceVector(i));
      //        i += 1
      //      }
      val derivative = new Array[Double](x.length)
      CubicSpline.computeHarmonicFirstDerivativePCHIM(x, priceArray, derivative)
      var pp = CubicSpline.makeHermiteSpline(x, priceArray, derivative);
      i = 0;
      while (i < grid.spaceVector.length) {
        var sgamma = pp.secondDerivative(x(i))
        var gamma = 0.
        if (i == 0 || i == grid.spaceVector.length - 1) {
        } else {
          //          delta = (priceArray(i + 1) - priceArray(i - 1)) / (x(i + 1) - x(i - 1))
          gamma = 2 * ((x(i + 1) - x(i)) * priceArray(i - 1) - (x(i + 1) - x(i - 1)) * priceArray(i) + (x(i) - x(i - 1)) * priceArray(i + 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)) * (x(i + 1) - x(i - 1)))
        }
        println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma)
        pw.println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma)
        if (i >= grid.spaceVector.length - 4) {
          //if it breaks there is a boundary condition error
          //          assert(gamma < 1e-3, "gamma small enough "+gamma+ " for "+schemeName)
        }
        i += 1
      }
    }
    pw.close()
  }

  test("gamma-itm-european") {
    val tte = 1.0;

    var spot = 100.0;
    var strike = 165.0;
    val vol = 0.4;
    val mu = 0.02;
    val r = 0.0;
    val spaceSize = 500;
    val timeSize = spaceSize / 64; //64 => small distorsion with scn & ran
    val isCall = false
    val payoff = new VanillaFDPayoff(isCall, strike, tte);
    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 2);

    //    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime,boundaries.bottomSpace , boundaries.topSpace )
    //    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime, Math.log(boundaries.bottomSpace / strike), Math.log(boundaries.topSpace / strike))
    val spaceBoundaries = new Mesh1DBoundaries(boundaries.bottomSpace, boundaries.topSpace)
    val timeBoundaries = new Mesh1DBoundaries(0, tte)

    val spaceMesh = new LogUniformMesh1D(spaceSize, spaceBoundaries, strike)
    val farAwayPoint = spot*math.exp(4*vol*math.sqrt(tte))
    val intermediatePoint = 0.5*(farAwayPoint+spaceMesh.x(spaceMesh.size-1))
    spaceMesh.insertPoints(Array(farAwayPoint))
    val grid: Mesh2D = new StaticAdaptiveMesh2D(
      spaceMesh,
      new UniformMesh1D(timeSize, timeBoundaries))
    //    val grid: UniformMesh2D = new UniformMesh2D(
    //      spaceSize,
    //      timeSize,
    //      logBounds,
    //      0.);

    val spec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r);
    //    grid.spaceTransform = new ExpTransformation(strike);
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("EULER", "EULERG", "CN", "TRBDF2", "RAN", "SCN", "LMG2", "LMG3")
    val f = new File("doc/lefloch_localextrapolation/gamma_european_itm.txt")
    val pw = new PrintWriter(f);
    pw.println("S\tScheme\tPrice\tGamma\tSGamma\tAGamma");

    for (schemeName <- methodNames) {
      val solver = new ElliotOckendonTridiagonalSolver(payoff);
      var method: Parabolic1DMethod = null
      schemeName match {
        case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
        case "RAN"    => method = new RannacherCentralBSM1FMethod(Array(tte), payoff)
        case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
        case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
        case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "EULERG" => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
        case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(tte), payoff)
      }
      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      var i = 0;
      var x = grid.spaceTransform.transform(grid.spaceVector);
      while (i < grid.spaceVector.length) {
        priceArray(i) = pricer.price(grid.spaceVector(i));
        i += 1
      }
      var pp = CubicSpline.makeCubicSpline(x, priceArray)
      i = 0;
      while (i < grid.spaceVector.length) {
        var gamma = 0.
        var delta = 0.
        var sdelta = 0.
        var sgamma = 0.
        if (i == 0 || i == grid.spaceVector.length - 1) {
        } else {
          delta = (priceArray(i + 1) - priceArray(i - 1)) / (x(i + 1) - x(i - 1))
          gamma = 2 * ((x(i + 1) - x(i)) * priceArray(i - 1) - (x(i + 1) - x(i - 1)) * priceArray(i) + (x(i) - x(i - 1)) * priceArray(i + 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)) * (x(i + 1) - x(i - 1)))
          sdelta = pp.derivative(x(i))
          sgamma = pp.secondDerivative(x(i))
          //          assert(math.abs(sdelta - delta)<3e-2, "first derivative "+sdelta+" "+delta)
          //          assert(math.abs(sgamma-gamma)<1e-2,"second derivative "+sgamma+" "+gamma)
        }

        val gammaN = delta

        val variance = BlackScholesVanillaEuropean.variance(vol, tte)
        val driftDf = BlackScholesVanillaEuropean.df(mu, tte)
        val discountDf = BlackScholesVanillaEuropean.df(r, tte)
        val priceAn = new Price()
        val gammaAn = new GammaGreek()
        val measures: Array[BSMMeasure] = Array(priceAn, gammaAn)
        BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, x(i), variance, driftDf, discountDf, measures)

        println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma + "\t" + gammaAn.value + "\t" + priceAn.value)
        pw.println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma + "\t" + gammaAn.value)
        i += 1
      }
    }
    pw.close()
  }

  test("gamma-european-call") {
    val tte = 0.25;

    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = -0.1;
    val r = 0.1;
    val spaceSize = 269;
    val timeSize = 100
    val isCall = true
    val payoff = new VanillaFDPayoff(isCall, strike, tte);
    val boundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3);

    //    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime,boundaries.bottomSpace , boundaries.topSpace )
    val logBounds = new MeshBoundaries(boundaries.firstTime, boundaries.lastTime, boundaries.bottomSpace, boundaries.topSpace)
    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      logBounds,
      strike);

    val spec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r);
    //    grid.spaceTransform = new ExpTransformation(strike);
    System.out.println("spaceVector=" + Arrays.toString(grid.spaceVector));

    val methodNames = Array("CN", "TRBDF2", "RAN", "SCN", "LMG2", "LMG3")
    val f = new File("doc/lefloch_localextrapolation/gamma_european_call.txt")
    val pw = new PrintWriter(f);
    pw.println("S\tScheme\tPrice\tGamma\tSGamma\tAGamma");

    for (schemeName <- methodNames) {
      //ElliotOckendon fails this test TODO see why!
      val solver = new ThomasTridiagonalSolver(payoff);
      var method: Parabolic1DMethod = null
      schemeName match {
        case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
        case "RAN"    => method = new RannacherCentralBSM1FMethod(Array(tte), payoff)
        case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
        case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
        case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
        case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(tte), payoff)
      }
      var pricer = new FDMSolver1D(spec, method, solver);
      pricer.solve(payoff);
      var priceArray = new Array[Double](grid.spaceVector.length);
      var i = 0;
      var x = grid.spaceTransform.transform(grid.spaceVector);
      while (i < grid.spaceVector.length) {
        priceArray(i) = pricer.price(grid.spaceVector(i));
        i += 1
      }
      var pp = CubicSpline.makeCubicSpline(x, priceArray)
      i = 0;
      while (i < grid.spaceVector.length) {
        var gamma = 0.
        var delta = 0.
        var sdelta = 0.
        var sgamma = 0.
        if (i == 0 || i == grid.spaceVector.length - 1) {
        } else {
          delta = (priceArray(i + 1) - priceArray(i - 1)) / (x(i + 1) - x(i - 1))
          gamma = 2 * ((x(i + 1) - x(i)) * priceArray(i - 1) - (x(i + 1) - x(i - 1)) * priceArray(i) + (x(i) - x(i - 1)) * priceArray(i + 1)) / ((x(i + 1) - x(i)) * (x(i) - x(i - 1)) * (x(i + 1) - x(i - 1)))
          sdelta = pp.derivative(x(i))
          sgamma = pp.secondDerivative(x(i))
          //          assert(math.abs(sdelta - delta)<3e-2, "first derivative "+sdelta+" "+delta)
          //          assert(math.abs(sgamma-gamma)<1e-2,"second derivative "+sgamma+" "+gamma)
        }

        val gammaN = delta

        val variance = BlackScholesVanillaEuropean.variance(vol, tte)
        val driftDf = BlackScholesVanillaEuropean.df(mu, tte)
        val discountDf = BlackScholesVanillaEuropean.df(r, tte)
        val priceAn = new Price()
        val gammaAn = new GammaGreek()
        val measures: Array[BSMMeasure] = Array(priceAn, gammaAn)
        BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, x(i), variance, driftDf, discountDf, measures)

        println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma + "\t" + gammaAn.value + "\t" + priceAn.value)
        pw.println(x(i) + "\t" + schemeName + "\t" + priceArray(i) + "\t" + gamma + "\t" + sgamma + "\t" + gammaAn.value)
        i += 1
      }
    }
    pw.close()
  }

}


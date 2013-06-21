package test.quantscale.fdm

import scala.collection.mutable.ArrayBuffer
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.analytic.AnalyticBarrierPayoff
import quantscale.analytic.BlackScholesMertonBarrier
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.mesh.Mesh1DBoundaries
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.fdm.mesh.Point
import quantscale.fdm.mesh.SMPoint
import quantscale.fdm.mesh.SinhMesh1D
import quantscale.fdm.mesh.SplineMesh1D
import quantscale.fdm.mesh.StaticAdaptiveMesh2D
import quantscale.fdm.mesh.UniformMesh1D
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.LRE3Parabolic1DMethod
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.RannacherCentralBSM1FMethod
import quantscale.fdm.method.SmoothCrankNicolsonCentralBSM1FMethod
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.payoff.KOBarrierFDPayoff
import quantscale.fdm.payoff.KOBarrierFDPayoff
import quantscale.fdm.payoff.KOBarrierFDPayoff
import quantscale.fdm.payoff.VanillaFDPayoff
import quantscale.analytic.BarrierType._
import quantscale.fdm.method.ForwardLinearOrder1Parabolic1DBoundaryFactory
import quantscale.fdm.method.BackwardLinearOrder1Parabolic1DBoundaryFactory
import quantscale.fdm.method.RK2Parabolic1DMethod
import quantscale.fdm.method.PredictorCorrectorParabolic1DMethod

@RunWith(classOf[JUnitRunner])
class BarrierSuite extends FunSuite {

  test("gobet-down-out") {
    val spot = 100.0
    val tte = 0.2
    val nKODays = 5
    val koTimes = new ArrayBuffer[Double]
    var i = 0
    while (i < nKODays) {
      koTimes += tte * (i + 1) / nKODays
      i += 1
    }
    val r = 0.1
    val q = 0.0
    val vol = 0.3
    val strike = spot
    val downBarrier = 0.99 * spot
    val upBarrier = Double.MaxValue
    if (upBarrier != Double.MaxValue) {
      val ybar = BlackScholesMertonBarrier.computeYbar(math.log(spot), math.log(upBarrier), vol, math.sqrt(tte / nKODays))
      val barrierShift = upBarrier * BlackScholesMertonBarrier.computeRelativeShift(math.sqrt(tte / nKODays), 0.0, vol, ybar)
      //      val anPrice = BlackScholesMertonBarrier.priceUpAndOut2(true, strike, barrierShift, spot, vol * vol * tte, math.exp(-(r - q) * tte), math.exp(-r * tte), math.exp(-(r - q) * tte), math.exp(-r * tte))
      val anPrice = BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(true, strike, barrierShift, UP_AND_OUT,0.0), spot, vol * vol * tte, math.exp(-(r - q) * tte), math.exp(-(r - q) * tte), math.exp(-r * tte), math.exp(-r * tte))

      println("UP=" + anPrice)
    } else {
      val ybar = BlackScholesMertonBarrier.computeYbar(math.log(spot), math.log(downBarrier), vol, math.sqrt(tte / nKODays))
      val barrierShift = downBarrier / BlackScholesMertonBarrier.computeRelativeShift(math.sqrt(tte / nKODays), 0.0, vol, ybar)
      val anPrice = BlackScholesMertonBarrier.priceDownAndOut(true, strike, barrierShift, spot, vol * vol * tte, math.exp(-(r - q) * tte), math.exp(-r * tte), math.exp(-(r - q) * tte), math.exp(-r * tte))

      println("DOWN=" + anPrice)
    }
    val payoff = new KOBarrierFDPayoff(true, strike, upBarrier, downBarrier, tte, koTimes.toArray, false)

    //the following makes order 4 for every scheme, except Rann at beginning
    //    val spaceSize= Array(31,62,125,250,500,1000,2000,4000)
    //    val timeSize = Array(nKODays, nKODays*2, nKODays*4, nKODays*8, nKODays*16,nKODays*32,nKODays*64,nKODays*128)

    //constant space, varying time
    val spaceSize = 
      Array(1000, 1000, 1000, 1000, 1000, 1000, 1000)
//      Array(100, 100, 100, 100, 100, 100, 100)
    val timeSize = Array(nKODays, nKODays * 2, nKODays * 4, nKODays * 8, nKODays * 16, nKODays * 64)

    //constant time, varying space => always order 2 except cn
    //    val spaceSize=Array(31,62,125,250,500,1000,2000,4000)  
    //    val timeSize = Array(nKODays,nKODays,nKODays,nKODays,nKODays,nKODays,nKODays,nKODays)
    //    val spaceSize = Array(62, 62, 62,      125, 125,  125,   250, 250, 250,  500, 500,  500,  1000, 1000, 1000, 2000,2000,2000,4000,4000, 4000)
    //    val timeSize = Array(nKODays, nKODays*2, nKODays*4, nKODays, nKODays*2, nKODays*4, nKODays, nKODays*2, nKODays*4,nKODays, nKODays*2, nKODays*4,nKODays, nKODays*2, nKODays*4,nKODays, nKODays*2,nKODays*4,nKODays, nKODays*2,nKODays*4 )
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "RAN")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "LMG3")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "TRBDF2")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "SCN")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "EULER")
//    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "EXP")
//        priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "RK2")
//                priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "PC")
  }

  test("discrete-knock-out") {
    val spot = 2016.0
    val tte = 1.0
    val nKODays = 120
    val koTimes = new ArrayBuffer[Double]
    var i = 0
    while (i < nKODays) {
      koTimes += tte * (i + 1) / nKODays
      i += 1
    }
    val r = 0.02
    val q = 0
    val vol = 0.3
    val strike = spot
    val downBarrier = Double.MinValue //0.8 * spot
    val upBarrier = 1.1 * spot // Double.MaxValue

    val ybar = BlackScholesMertonBarrier.computeYbar(math.log(spot), math.log(upBarrier), vol, math.sqrt(tte / nKODays))
    val barrierShift = upBarrier * BlackScholesMertonBarrier.computeRelativeShift(math.sqrt(tte / nKODays), 0.0, vol, ybar)
    val anPrice = BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(true, strike, barrierShift, UP_AND_OUT), spot, vol * vol * tte, math.exp(-(r - q) * tte), math.exp(-(r - q) * tte), math.exp(-r * tte),math.exp(-r * tte))
    println("AN=" + anPrice)
    val payoff = new KOBarrierFDPayoff(true, strike, upBarrier, downBarrier, tte, koTimes.toArray, false)

    //the following makes order 4 for every scheme, except Rann at beginning
    //    val spaceSize= Array(31,62,125,250,500,1000,2000,4000)
    //    val timeSize = Array(nKODays, nKODays*2, nKODays*4, nKODays*8, nKODays*16,nKODays*32,nKODays*64,nKODays*128)

    //constant space, varying time
    val spaceSize = Array(1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000)
    val timeSize = Array(nKODays, nKODays * 2, nKODays * 4, nKODays * 8, nKODays * 16, nKODays * 32, nKODays * 64)

    //constant time, varying space => always order 2 except cn
    //    val spaceSize=Array(31,62,125,250,500,1000,2000,4000)  
    //    val timeSize = Array(nKODays,nKODays,nKODays,nKODays,nKODays,nKODays,nKODays,nKODays)
    //    val spaceSize = Array(62, 62, 62,      125, 125,  125,   250, 250, 250,  500, 500,  500,  1000, 1000, 1000, 2000,2000,2000,4000,4000, 4000)
    //    val timeSize = Array(nKODays, nKODays*2, nKODays*4, nKODays, nKODays*2, nKODays*4, nKODays, nKODays*2, nKODays*4,nKODays, nKODays*2, nKODays*4,nKODays, nKODays*2, nKODays*4,nKODays, nKODays*2,nKODays*4,nKODays, nKODays*2,nKODays*4 )
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "CN")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "RAN")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "SCN")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "TRBDF2")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "LMG2")
    priceKO(spot, vol, r - q, r, payoff, spaceSize, timeSize, "LMG3")

  }

  def priceKO(spot: Double, vol: Double, mu: Double, r: Double, payoff: KOBarrierFDPayoff, spaceSize: Array[Int], timeSize: Array[Int], schemeName: String) {
    val upBarrier = payoff.upLevel
    val downBarrier = payoff.downLevel
    val tte = payoff.expiryTime
    val meshBoundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, math.exp(-mu * tte), 4)
    //    val spaceMesh = new UniformMesh1D(16000, new Mesh1DBoundaries(meshBoundaries.bottomSpace, meshBoundaries.topSpace))
    //    val timeMesh = new UniformMesh1D(1000, new Mesh1DBoundaries(0, tte))

    val starPoint = if (Math.abs(upBarrier - spot) < Math.abs(downBarrier - spot)) upBarrier else downBarrier

    val specialPoints =
      if (upBarrier < Double.MaxValue)
        if (upBarrier > payoff.strike) Array(new Point(payoff.strike, false), new Point(upBarrier, true))
        else Array(new Point(upBarrier, true), new Point(payoff.strike, false))
      else if (downBarrier < payoff.strike) Array(new Point(downBarrier, true), new Point(payoff.strike, false))
      else Array(new Point(payoff.strike, false), new Point(downBarrier, true))
    println(schemeName)
    var previousPrice = 0.
    val timePoints = new Array[SMPoint](payoff.koTimes.length)
    for (i <- 0 until timePoints.length) {
      timePoints(i) = new SMPoint(payoff.koTimes(i))
    }
    var previousChange = 0.
    for (i <- 0 until timeSize.length) {
      var tSize = timeSize(i)
      var method: Parabolic1DMethod = null
      schemeName match {
        case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
        case "RAN"    => method = new RannacherCentralBSM1FMethod(payoff.koTimes)
        case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
        case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
        case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
        case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
        case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(payoff.koTimes, payoff)
        case "EXP" =>
          method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_EXPLICIT)
          tSize = tSize * spaceSize(i)
        case "RK2" => 
          method = new RK2Parabolic1DMethod(payoff)
        tSize = tSize * spaceSize(i)
        case "PC" =>
          method = new PredictorCorrectorParabolic1DMethod(payoff)
          tSize = tSize * spaceSize(i)
      }
      val spaceMesh = 
                new SinhMesh1D(spaceSize(i), new Mesh1DBoundaries(meshBoundaries.bottomSpace, meshBoundaries.topSpace), false, specialPoints, starPoint, SinhMesh1D.ALPHA_ALMOST_UNIFORM)
//        new SinhMesh1D(spaceSize(i), new Mesh1DBoundaries(meshBoundaries.bottomSpace, meshBoundaries.topSpace), false, specialPoints, starPoint, SinhMesh1D.ALPHA_HIGHLY_NONUNIFORM / 4)

      val timeMesh = new SplineMesh1D(tSize, new Mesh1DBoundaries(0, tte), timePoints)
//            val timeMesh = new UniformMesh1D(timeSize(i), new Mesh1DBoundaries(0, tte))
//            timeMesh.insertPoints(payoff.koTimes.toArray)
      val mesh = new StaticAdaptiveMesh2D(spaceMesh, timeMesh)
      val spec = new ConstantBSM1FFDSpec(
        mesh,
        vol,
        mu,
        r)
      

      method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
      method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory
      var solver = new ThomasTridiagonalSolver(payoff)
      var pricer = new FDMSolver1D(spec, method, solver)
      var rtTime = System.nanoTime();
      pricer.solve(payoff)
      val price = pricer.price(spot)
      val time = (System.nanoTime() - rtTime) * 1e-9
      val change = price - previousPrice
      val ratio = previousChange / change
      previousPrice = price
      previousChange = change

      val lineStr = "%s & %s & %.6f & %.2e & %.1f & %.2e".format(spaceSize(i), tSize, price, change, ratio, time);
      println(lineStr)
    }
  }
}
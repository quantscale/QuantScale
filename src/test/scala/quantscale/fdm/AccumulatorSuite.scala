package test.quantscale.fdm

import java.io.File
import java.io.PrintWriter
import scala.collection.mutable.ArrayBuffer
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import javax.time.calendar.LocalDate
import javax.time.calendar.LocalDateTime
import javax.time.calendar.Period
import quantscale.analytic.BarrierType._
import quantscale.analytic.BlackScholesMertonBarrier
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.RE2FDMSolver1D
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.listener.GammaListener
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.fdm.mesh.Point
import quantscale.fdm.mesh.SMPoint
import quantscale.fdm.mesh.SinhMesh1D
import quantscale.fdm.mesh.SplineMesh1D
import quantscale.fdm.mesh.StaticAdaptiveMesh2D
import quantscale.fdm.mesh.UniformMesh1D
import quantscale.fdm.method.BackwardLinearOrder1Parabolic1DBoundaryFactory
import quantscale.fdm.method.ForwardLinearOrder1Parabolic1DBoundaryFactory
import quantscale.fdm.method.LMG2Parabolic1DMethod
import quantscale.fdm.method.LMG3Parabolic1DMethod
import quantscale.fdm.method.LRE3Parabolic1DMethod
import quantscale.fdm.method.Parabolic1DMethod
import quantscale.fdm.method.RannacherCentralBSM1FMethod
import quantscale.fdm.method.SmoothCrankNicolsonCentralBSM1FMethod
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.method.ThetaParabolic1DMethod
import quantscale.fdm.payoff.AccumulatorFDPayoff
import quantscale.fdm.payoff.AccumulatorKODABisFDPayoff
import quantscale.fdm.payoff.AccumulatorKODAFDPayoff
import quantscale.fdm.payoff.AccumulatorTARNFDPayoff
import quantscale.fdm.payoff.BackwardSchedule
import quantscale.fdm.RE3FDMSolver1D
import quantscale.fdm.method.TGAParabolic1DMethod
import quantscale.fdm.method.TGASerialParabolic1DMethod
import quantscale.fdm.method.TGABisParabolic1DMethod

@RunWith(classOf[JUnitRunner])
class AccumulatorSuite extends FunSuite {

  private def diffACT365(begin: LocalDate, end: LocalDate): Double = {
    return Period.daysBetween(begin, end).getDays() / 365.0
  }

  private def diffACT365(begin: LocalDateTime, end: LocalDateTime): Double = {
    val timeDay = (end.getHourOfDay() - begin.getHourOfDay()) / 24.0 + (end.getMinuteOfHour() - begin.getMinuteOfHour()) / (24 * 60.0)
    return (Period.daysBetween(begin, end).getDays() + timeDay) / 365.0
  }

  test("accumulator-greeks") {
    val r = 0
    val mu = 0
    val vol = 0.3
    val buySell = 1
    val spot = 2016.0
    val forwardPrice = 0.9 * spot
    val accrualAbove = 1200.0
    val accrualBelow = 600.0
    val knockOutBarrier = 1.1 * spot

    val startDate = LocalDate.of(2010, 9, 1)
    val valueDate = startDate
    val accrualDelay = 1
    val paymentDelay = 30
    val nPeriods = 12
    val maturityDate = valueDate.plusDays(nPeriods * paymentDelay)

    val knockouts = new ArrayBuffer[Double]()
    val payments = new ArrayBuffer[Double]()
    var tempDate = startDate.plusDays(paymentDelay)
    for (i <- 0 until nPeriods) {
      payments += diffACT365(valueDate, tempDate)
      tempDate = tempDate.plusDays(paymentDelay)
    }

    payments += diffACT365(valueDate, maturityDate)
    val periodAccruals = new ArrayBuffer[Double]()
    tempDate = startDate.plusDays(accrualDelay);
    while (tempDate.isBefore(maturityDate) || tempDate.equals(maturityDate)) {
      val yearFrac = diffACT365(valueDate, tempDate)
      periodAccruals += yearFrac
      knockouts += yearFrac
      tempDate = tempDate.plusDays(accrualDelay);
    }

    val paymentsSchedule = new BackwardSchedule(payments.toArray)
    val periodAccrualsSchedule = new BackwardSchedule(periodAccruals.toArray)
    val knockoutsSchedule = new BackwardSchedule(knockouts.toArray)

    val payoff = new AccumulatorFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier)

    val fdSchemes = Array("SCN", "LMG3", "RAN", "EULER", "CN", "TRBDF2")
    val tSize = Array(360 + 1)
    val xSize = Array(2000)
    //    val fdSchemes = Array("TRBDF2", "SCN", "LMG2", "LMG3", "EULER", "RAN", "CN")

    val tte = diffACT365(valueDate, maturityDate)
    val meshBoundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, math.exp(-mu * tte), 4)
    val f = new File("doc/lefloch_localextrapolation/gamma_accumulator.txt")
    val pw = new PrintWriter(f);
    pw.println("Scheme S Price Delta Gamma")
    for (timeSize <- tSize) {
      for (scheme <- fdSchemes) {
        var method: Parabolic1DMethod = null

        for (spaceSize <- xSize) {
          val spaceMesh = new SplineMesh1D(
            spaceSize, meshBoundaries.spaceBoundaries,
            Array(new SMPoint(knockOutBarrier, true, spot * 0.05, 0.1)))
          //        val spaceMesh = new UniformMesh1D(spaceSize, meshBoundaries.spaceBoundaries)
          val timeMesh = new UniformMesh1D(timeSize, meshBoundaries.timeBoundaries)

          //todo should be done by solver aware of payoff? or mesh.initPayoff? 
          timeMesh.insertPoints(payments.toArray)
          timeMesh.insertPoints(periodAccruals.toArray)
          timeMesh.insertPoints(knockouts.toArray)

          val mesh = new StaticAdaptiveMesh2D(spaceMesh, timeMesh)
          val spec = new ConstantBSM1FFDSpec(mesh, vol, mu, r)
          scheme match {
            case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
            case "RAN"    => method = new RannacherCentralBSM1FMethod(payoff.knockouts.times, payoff)
            case "CN"     => method = new ThetaParabolic1DMethod()
            case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
            case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
            case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
            case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(payoff.knockouts.times, payoff)
          }
          //     method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
          method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
          method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory

          var solver = new ThomasTridiagonalSolver(payoff)
          var pricer = new FDMSolver1D(spec, method, solver)
          val myListener = new GammaListener(Array(0.0))
          pricer.listener = myListener
          var rtTime = System.nanoTime();
          pricer.solve(payoff)
          val time = (System.nanoTime() - rtTime) * 1e-9
          var i = 1
          while (i < spaceSize - 1) {
            val spot = spaceMesh.x(i)
            val price = pricer.price(i)
            val delta = myListener.delta(i)
            val gamma = myListener.gamma(i)
            pw.println(scheme + " " + spot + " " + price + " " + delta + " " + gamma)
            i += 1
          }
        }
      }
    }
    pw.close()
  }

  test("accumulator-nikkei", SpecificTest) {
    val r = 0
    val mu = 0
    val vol = 0.3
    val buySell = 1
    val spot = 2016.0
    val forwardPrice = 0.9 * spot
    val accrualAbove = 1200.0
    val accrualBelow = 600.0
    val knockOutBarrier = 1.1 * spot

    val startDate = LocalDateTime.of(2010, 9, 1, 0, 0, 0)
    val valueDate = startDate
    val accrualDelay = 1
    val paymentDelay = 30
    val nPeriods = 12
    val isKoAfterAccrual = false
    val maturityDate = valueDate.plusDays(nPeriods * paymentDelay)

    val knockouts = new ArrayBuffer[Double]()
    val payments = new ArrayBuffer[Double]()
    var tempDate = startDate.plusDays(paymentDelay)
    for (i <- 0 until nPeriods) {
      payments += diffACT365(valueDate, tempDate)
      tempDate = tempDate.plusDays(paymentDelay)
    }

    payments += diffACT365(valueDate, maturityDate)
    val periodAccruals = new ArrayBuffer[Double]()
    tempDate = startDate.plusHours(math.round(accrualDelay * 24).toInt);
    while (tempDate.isBefore(maturityDate) || tempDate.equals(maturityDate)) {
      val yearFrac = diffACT365(valueDate, tempDate)
      periodAccruals += yearFrac
      knockouts += yearFrac
      tempDate = tempDate.plusHours(math.round(accrualDelay * 24).toInt);
    }

    val paymentsSchedule = new BackwardSchedule(payments.toArray)
    val periodAccrualsSchedule = new BackwardSchedule(periodAccruals.toArray)
    val knockoutsSchedule = new BackwardSchedule(knockouts.toArray)

    val anPrice = BlackScholesMertonBarrier.priceAccumulator(buySell, isKoAfterAccrual, periodAccruals, knockouts, payments, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, spot, vol, r, mu)
    println("AN=" + anPrice)
    val tte = diffACT365(valueDate, maturityDate)
    //
    //    val sqrtKoFrequency = math.sqrt(tte/knockouts.length)
    //    val koVariance = tte*vol*vol
    //    val ybar= BlackScholesMertonBarrier.computeYbar(math.log(spot), math.log(knockOutBarrier), vol, sqrtKoFrequency)
    //    val shift = BlackScholesMertonBarrier.computeRelativeShift(sqrtKoFrequency, 0.0, vol, ybar)
    //    val koBeforePrice = accrualAbove*BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(true, forwardPrice, knockOutBarrier*shift,UP_AND_IN), spot, koVariance, 1.0, 1.0,1.0,1.0)-
    //    accrualBelow*BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(false, forwardPrice, knockOutBarrier*shift, UP_AND_IN), spot, koVariance, 1.0, 1.0,1.0,1.0)
    //    println("AN+="+(anPrice+koBeforePrice))
    val payoff = new AccumulatorKODAFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKoAfterAccrual)

    val tSize = Array(360 + 1, 360 * 2 + 1, 360 * 4 + 1) //Array(36 * 64 + 1, 36 * 128 + 1)
    val xSize = Array(801, 3201)
    val fdSchemes = Array("TGAb","TGA")//Array("TRBDF2", "SCN", "LMG2", "RE2", "TGAb", "TGA","TGAS","LMG3" , "LRE3", "RE3")
    //    val tSize = Array(360 + 1, 360 * 2 + 1, 360 * 4 + 1, 360 * 8 + 1, 360 * 16 + 1)
    //    val xSize = Array(125, 250, 500, 1000, 2000, 4000, 8000)
    //    val fdSchemes = Array("TRBDF2", "SCN", "LMG2", "LMG3", "EULER", "RAN", "CN")

    val meshBoundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, math.exp(-mu * tte), 4)

    for (timeSize <- tSize) {
      for (scheme <- fdSchemes) {
        var method: Parabolic1DMethod = null

        for (spaceSize <- xSize) {
          //          val spaceMesh = new SplineMesh1D(
          //            spaceSize, meshBoundaries.spaceBoundaries,
          //            Array(new SMPoint(knockOutBarrier, true, spot * 0.05, 0.1)))

          val spaceMesh = new SinhMesh1D(spaceSize, meshBoundaries.spaceBoundaries, false, Array(new Point(knockOutBarrier, true)), knockOutBarrier, SinhMesh1D.ALPHA_HIGHLY_NONUNIFORM)
          //        val spaceMesh = new UniformMesh1D(spaceSize, meshBoundaries.spaceBoundaries)
          val timeMesh = new UniformMesh1D(timeSize, meshBoundaries.timeBoundaries)

          //todo should be done by solver aware of payoff? or mesh.initPayoff? 
          timeMesh.insertPoints(payments.toArray)
          timeMesh.insertPoints(periodAccruals.toArray)
          timeMesh.insertPoints(knockouts.toArray)

          val mesh = new StaticAdaptiveMesh2D(spaceMesh, timeMesh)
          val spec = new ConstantBSM1FFDSpec(mesh, vol, mu, r)
          scheme match {
            case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
            case "RAN"    => method = new RannacherCentralBSM1FMethod(payoff.knockouts.times, payoff)
            case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
            case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
            case "EULER" | "RE2" | "RE3"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
            case "LMG3"   => method = new LMG3Parabolic1DMethod(payoff)
            case "LRE3"   => method = new LRE3Parabolic1DMethod(payoff)
            case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(payoff.knockouts.times, payoff)
            case "TGA" => method = new TGAParabolic1DMethod(payoff,2-math.sqrt(2))
            case "TGAb" => method = new TGABisParabolic1DMethod(payoff,1.999999-math.sqrt(2))
            case "TGAS" => method = new TGASerialParabolic1DMethod(payoff,2-math.sqrt(2))
          }
          //     method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
          method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
          method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory

          
          var solver = new ThomasTridiagonalSolver(payoff)
          var pricer =new FDMSolver1D(spec, method, solver) 
          scheme match {
            case "RE2" => pricer = new RE2FDMSolver1D(pricer)
            case "RE3" => pricer = new RE3FDMSolver1D(pricer)
            case _ =>
          }
          var rtTime = System.nanoTime();
          pricer.solve(payoff)
          val price = pricer.price(spot)
          val time = (System.nanoTime() - rtTime) * 1e-9

          println(scheme + " " + timeSize + " " + timeMesh.size + " " + spaceSize + " " + price + " " + time)
        }
      }
    }
    assert(true)
  }

  test("accumulator-nikkei-with-rate", SpecificTest) {
    val r = 0.1
    val mu = 0.1
    val vol = 0.3
    val buySell = 1
    val spot = 2016.0
    val forwardPrice = 0.9 * spot
    val accrualAbove = 1200.0
    val accrualBelow = 600.0
    val knockOutBarrier = 1.1 * spot

    val startDate = LocalDateTime.of(2010, 9, 1, 0, 0, 0)
    val valueDate = startDate
    val accrualDelay = 1
    val paymentDelay = 30
    val nPeriods = 12
    val isKoAfterAccrual = true
    val maturityDate = valueDate.plusDays(nPeriods * paymentDelay)

    val knockouts = new ArrayBuffer[Double]()
    val payments = new ArrayBuffer[Double]()
    var tempDate = startDate.plusDays(paymentDelay)
    for (i <- 0 until nPeriods) {
      payments += diffACT365(valueDate, tempDate)
      tempDate = tempDate.plusDays(paymentDelay)
    }

    payments += diffACT365(valueDate, maturityDate)
    val periodAccruals = new ArrayBuffer[Double]()
    tempDate = startDate.plusHours(math.round(accrualDelay * 24).toInt);
    while (tempDate.isBefore(maturityDate) || tempDate.equals(maturityDate)) {
      val yearFrac = diffACT365(valueDate, tempDate)
      periodAccruals += yearFrac
      knockouts += yearFrac
      tempDate = tempDate.plusHours(math.round(accrualDelay * 24).toInt);
    }

    val paymentsSchedule = new BackwardSchedule(payments.toArray)
    val periodAccrualsSchedule = new BackwardSchedule(periodAccruals.toArray)
    val knockoutsSchedule = new BackwardSchedule(knockouts.toArray)

    val anPrice = BlackScholesMertonBarrier.priceAccumulator(buySell, isKoAfterAccrual, periodAccruals, knockouts, payments, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, spot, vol, r, mu)
    println("AN=" + anPrice)
    val tte = diffACT365(valueDate, maturityDate)
    //
    //    val sqrtKoFrequency = math.sqrt(tte/knockouts.length)
    //    val koVariance = tte*vol*vol
    //    val ybar= BlackScholesMertonBarrier.computeYbar(math.log(spot), math.log(knockOutBarrier), vol, sqrtKoFrequency)
    //    val shift = BlackScholesMertonBarrier.computeRelativeShift(sqrtKoFrequency, 0.0, vol, ybar)
    //    val koBeforePrice = accrualAbove*BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(true, forwardPrice, knockOutBarrier*shift,UP_AND_IN), spot, koVariance, 1.0, 1.0,1.0,1.0)-
    //    accrualBelow*BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(false, forwardPrice, knockOutBarrier*shift, UP_AND_IN), spot, koVariance, 1.0, 1.0,1.0,1.0)
    //    println("AN+="+(anPrice+koBeforePrice))
    val payoffRepl = new AccumulatorKODAFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKoAfterAccrual)
    val payoffDummy = new AccumulatorKODABisFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKoAfterAccrual, Double.MaxValue)
//  val payoffTest = new AccumulatorKODATerFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKoAfterAccrual)

    val payoffs = Array(payoffRepl, payoffDummy)
    val tSize = Array(360 + 1, 360 * 2 + 1, 360 * 4 + 1) //Array(36 * 64 + 1, 36 * 128 + 1)
    val xSize = Array(801, 3201)
    val fdSchemes = Array("TRBDF2", "SCN", "LMG3")
    //    val tSize = Array(360 + 1, 360 * 2 + 1, 360 * 4 + 1, 360 * 8 + 1, 360 * 16 + 1)
    //    val xSize = Array(125, 250, 500, 1000, 2000, 4000, 8000)
    //    val fdSchemes = Array("TRBDF2", "SCN", "LMG2", "LMG3", "EULER", "RAN", "CN")

    val meshBoundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, math.exp(-mu * tte), 4)

    for (timeSize <- tSize) {
      for (scheme <- fdSchemes) {
        var method: Parabolic1DMethod = null

        for (spaceSize <- xSize) {
          for (payoff <- payoffs) {
            //          val spaceMesh = new SplineMesh1D(
            //            spaceSize, meshBoundaries.spaceBoundaries,
            //            Array(new SMPoint(knockOutBarrier, true, spot * 0.05, 0.1)))

            val spaceMesh = new SinhMesh1D(spaceSize, meshBoundaries.spaceBoundaries, false, Array(new Point(knockOutBarrier, true)), knockOutBarrier, SinhMesh1D.ALPHA_HIGHLY_NONUNIFORM)
            //        val spaceMesh = new UniformMesh1D(spaceSize, meshBoundaries.spaceBoundaries)
            val timeMesh = new UniformMesh1D(timeSize, meshBoundaries.timeBoundaries)

            //todo should be done by solver aware of payoff? or mesh.initPayoff? 
            timeMesh.insertPoints(payments.toArray)
            timeMesh.insertPoints(periodAccruals.toArray)
            timeMesh.insertPoints(knockouts.toArray)

            val mesh = new StaticAdaptiveMesh2D(spaceMesh, timeMesh)
            val spec = new ConstantBSM1FFDSpec(mesh, vol, mu, r)
            scheme match {
              case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
              case "RAN"    => method = new RannacherCentralBSM1FMethod(Array(), payoff)
              case "CN"     => method = new ThetaParabolic1DMethod()
              case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
              case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
              case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
              case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(), payoff)
            }
            //     method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
            method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
            method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory

            var solver = new ThomasTridiagonalSolver(payoff)
            var pricer = new FDMSolver1D(spec, method, solver)
            var rtTime = System.nanoTime();
            pricer.solve(payoff)
            val price = pricer.price(spot)
            val time = (System.nanoTime() - rtTime) * 1e-9

            println(payoff.getClass().getSimpleName() + " " + scheme + " " + timeSize + " " + timeMesh.size + " " + spaceSize + " " + price + " " + time)
          }
        }
      }
    }
    assert(true)
  }

  test("accumulator-tarn-nikkei", SpecificTest) {
    val r = 0
    val mu = 0
    val vol = 0.3
    val buySell = 1
    val spot = 2016.0
    val forwardPrice = 0.9 * spot
    val accrualAbove = 1200.0
    val accrualBelow = 600.0
    val knockOutBarrier = 1.1 * spot
    val tarnLevel = 300000 //300000 12 expecting 360 4229786 720 4225446
    //60000 12 expecting 720 8004804
    //30000 2 expecting 5229465
    //10000 1 expecting 483839.98
    //1000 1 //expecting 241919.97294875103 
    //1000*30.0 1 //expecting 5249439

    val startDate = LocalDateTime.of(2010, 9, 1, 0, 0, 0)
    val valueDate = startDate
    val accrualDelay = 1
    val paymentDelay = 30
    val nPeriods = 12
    val isKoAfterAccrual = true
    val maturityDate = valueDate.plusDays(nPeriods * paymentDelay)

    val knockouts = new ArrayBuffer[Double]()
    val payments = new ArrayBuffer[Double]()
    var tempDate = startDate.plusDays(paymentDelay)
    for (i <- 0 until nPeriods) {
      payments += diffACT365(valueDate, tempDate)
      tempDate = tempDate.plusDays(paymentDelay)
    }

    payments += diffACT365(valueDate, maturityDate)
    val periodAccruals = new ArrayBuffer[Double]()
    tempDate = startDate.plusHours(math.round(accrualDelay * 24).toInt);
    while (tempDate.isBefore(maturityDate) || tempDate.equals(maturityDate)) {
      val yearFrac = diffACT365(valueDate, tempDate)
      periodAccruals += yearFrac
      knockouts += yearFrac
      tempDate = tempDate.plusHours(math.round(accrualDelay * 24).toInt);
    }

    val paymentsSchedule = new BackwardSchedule(payments.toArray)
    val periodAccrualsSchedule = new BackwardSchedule(periodAccruals.toArray)
    val knockoutsSchedule = new BackwardSchedule(knockouts.toArray)

    val anPrice = BlackScholesMertonBarrier.priceAccumulator(buySell, isKoAfterAccrual, periodAccruals, knockouts, payments, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, spot, vol, r, mu)
    println("AN=" + anPrice)
    val tte = diffACT365(valueDate, maturityDate)
    //
    //    val sqrtKoFrequency = math.sqrt(tte/knockouts.length)
    //    val koVariance = tte*vol*vol
    //    val ybar= BlackScholesMertonBarrier.computeYbar(math.log(spot), math.log(knockOutBarrier), vol, sqrtKoFrequency)
    //    val shift = BlackScholesMertonBarrier.computeRelativeShift(sqrtKoFrequency, 0.0, vol, ybar)
    //    val koBeforePrice = accrualAbove*BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(true, forwardPrice, knockOutBarrier*shift,UP_AND_IN), spot, koVariance, 1.0, 1.0,1.0,1.0)-
    //    accrualBelow*BlackScholesMertonBarrier.price(new AnalyticBarrierPayoff(false, forwardPrice, knockOutBarrier*shift, UP_AND_IN), spot, koVariance, 1.0, 1.0,1.0,1.0)
    //    println("AN+="+(anPrice+koBeforePrice))
    val payoffTARNKODA =
      new AccumulatorTARNFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKoAfterAccrual, tarnLevel)
    val payoffTARN = new AccumulatorKODABisFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, isKoAfterAccrual, tarnLevel)
    val payoffs = Array( payoffTARNKODA)
    val tSize = Array(paymentDelay * nPeriods + 1, paymentDelay * nPeriods * 2 + 1, paymentDelay * nPeriods * 4 + 1) //Array(36 * 64 + 1, 36 * 128 + 1)
    val xSize = Array(801, 3201)
    val fdSchemes = Array("EULER", "CN", "SCN", "LMG2", "RE2", "LMG3", "LRE3", "RE3")
    //    val tSize = Array(360 + 1, 360 * 2 + 1, 360 * 4 + 1, 360 * 8 + 1, 360 * 16 + 1)
    //    val xSize = Array(125, 250, 500, 1000, 2000, 4000, 8000)
    //    val fdSchemes = Array("TRBDF2", "SCN", "LMG2", "LMG3", "EULER", "RAN", "CN")

    val meshBoundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, math.exp(-mu * tte), 4)

    for (timeSize <- tSize) {
      for (scheme <- fdSchemes) {
        for (payoff <- payoffs) {
          var method: Parabolic1DMethod = null

          for (spaceSize <- xSize) {
            //          val spaceMesh = new SplineMesh1D(spaceSize, meshBoundaries.spaceBoundaries,
            //            Array(new SMPoint(knockOutBarrier, true, spot * 0.05, 0.1)))

            val spaceMesh = new SinhMesh1D(spaceSize, meshBoundaries.spaceBoundaries, false, Array(new Point(knockOutBarrier, true)), knockOutBarrier, SinhMesh1D.ALPHA_HIGHLY_NONUNIFORM)

            //                  val spaceMesh = new UniformMesh1D(spaceSize, meshBoundaries.spaceBoundaries)
            val timeMesh = new UniformMesh1D(timeSize, meshBoundaries.timeBoundaries)

            //todo should be done by solver aware of payoff? or mesh.initPayoff? 
            timeMesh.insertPoints(payments.toArray)
            timeMesh.insertPoints(periodAccruals.toArray)
            timeMesh.insertPoints(knockouts.toArray)

            val mesh = new StaticAdaptiveMesh2D(spaceMesh, timeMesh)
            val spec = new ConstantBSM1FFDSpec(mesh, vol, mu, r)
             scheme match {
            case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
            case "RAN"    => method = new RannacherCentralBSM1FMethod(Array(), payoff)
            case "CN"     => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_CRANK_NICOLSON)
            case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
            case "EULER" | "RE2" | "RE3"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
            case "LMG3"   => method = new LMG3Parabolic1DMethod(payoff)
            case "LRE3"   => method = new LRE3Parabolic1DMethod(payoff)
            case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(Array(), payoff)
          }
            //     method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
            method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
            method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory

            var solver = new ThomasTridiagonalSolver(payoff)
            var pricer = new FDMSolver1D(spec, method, solver)
                   scheme match {
            case "RE2" => pricer = new RE2FDMSolver1D(pricer)
            case "RE3" => pricer = new RE3FDMSolver1D(pricer)
            case _ =>
          }
            var rtTime = System.nanoTime();
            pricer.solve(payoff)
            val price = pricer.price(spot)
            val time = (System.nanoTime() - rtTime) * 1e-9

            println(payoff.getClass().getSimpleName() + " " + scheme + " " + timeSize + " " + timeMesh.size + " " + spaceSize + " " + price + " " + time)
          }
        }
      }
    }
    assert(true)
  }

  test("accumulator-anal", SpecificTest) {
    val r = 0.00
    val mu = r
    val vol = 0.2
    val buySell = 1
    val spot = 100.0
    val forwardPrice = 0.9 * spot
    val accrualAbove = 1.0
    val accrualBelow = 2.0
    val knockOutBarrier = 1.05 * spot

    val startDate = LocalDateTime.of(2010, 9, 1, 0, 0, 0)
    val valueDate = startDate
    val accrualDelay = 1
    val paymentDelay = 21
    val nPeriods = 12
    val maturityDate = valueDate.plusDays(nPeriods * paymentDelay)

    val knockouts = new ArrayBuffer[Double]()
    val payments = new ArrayBuffer[Double]()
    var tempDate = startDate.plusDays(paymentDelay)
    for (i <- 0 until nPeriods) {
      payments += diffACT365(valueDate, tempDate) * 365.0 / 252.0
      tempDate = tempDate.plusDays(paymentDelay)
    }

    payments += diffACT365(valueDate, maturityDate) * 365.0 / 252.0
    val periodAccruals = new ArrayBuffer[Double]()
    tempDate = startDate.plusHours(math.round(accrualDelay * 24).toInt);
    while (tempDate.isBefore(maturityDate) || tempDate.equals(maturityDate)) {
      val yearFrac = diffACT365(valueDate, tempDate) * 365.0 / 252.0
      periodAccruals += yearFrac
      knockouts += yearFrac
      tempDate = tempDate.plusHours(math.round(accrualDelay * 24).toInt);
    }

    val paymentsSchedule = new BackwardSchedule(payments.toArray)
    val periodAccrualsSchedule = new BackwardSchedule(periodAccruals.toArray)
    val knockoutsSchedule = new BackwardSchedule(knockouts.toArray)

    val anPrice = BlackScholesMertonBarrier.priceAccumulator(buySell, false, periodAccruals, knockouts, payments, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, spot, vol, r, mu)
    println("AN=" + anPrice)
    val payoff = new AccumulatorKODAFDPayoff(buySell, periodAccrualsSchedule, knockoutsSchedule, paymentsSchedule, forwardPrice, accrualAbove, accrualBelow, knockOutBarrier, false)

    val tSize = Array(252 + 1, 252 * 2 + 1, 252 * 4 + 1) //Array(36 * 64 + 1, 36 * 128 + 1)
    val xSize = Array(801, 3201)
    val fdSchemes = Array("TRBDF2", "SCN", "LMG3")
    //    val tSize = Array(360 + 1, 360 * 2 + 1, 360 * 4 + 1, 360 * 8 + 1, 360 * 16 + 1)
    //    val xSize = Array(125, 250, 500, 1000, 2000, 4000, 8000)
    //    val fdSchemes = Array("TRBDF2", "SCN", "LMG2", "LMG3", "EULER", "RAN", "CN")

    val tte = diffACT365(valueDate, maturityDate) * 365.0 / 252.0
    val meshBoundaries = MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, math.exp(-mu * tte), 4)

    for (timeSize <- tSize) {
      for (scheme <- fdSchemes) {
        var method: Parabolic1DMethod = null

        for (spaceSize <- xSize) {
          //          val spaceMesh = new SplineMesh1D(
          //            spaceSize, meshBoundaries.spaceBoundaries,
          //            Array(new SMPoint(knockOutBarrier, true, spot * 0.05, 0.1)))

          val spaceMesh = new SinhMesh1D(spaceSize, meshBoundaries.spaceBoundaries, false, Array(new Point(knockOutBarrier, true)), knockOutBarrier, SinhMesh1D.ALPHA_HIGHLY_NONUNIFORM)
          //        val spaceMesh = new UniformMesh1D(spaceSize, meshBoundaries.spaceBoundaries)
          val timeMesh = new UniformMesh1D(timeSize, meshBoundaries.timeBoundaries)

          //todo should be done by solver aware of payoff? or mesh.initPayoff? 
          timeMesh.insertPoints(payments.toArray)
          timeMesh.insertPoints(periodAccruals.toArray)
          timeMesh.insertPoints(knockouts.toArray)

          val mesh = new StaticAdaptiveMesh2D(spaceMesh, timeMesh)
          val spec = new ConstantBSM1FFDSpec(mesh, vol, mu, r)
          scheme match {
            case "LMG2"   => method = new LMG2Parabolic1DMethod(payoff)
            case "RAN"    => method = new RannacherCentralBSM1FMethod(payoff.knockouts.times, payoff)
            case "CN"     => method = new ThetaParabolic1DMethod()
            case "TRBDF2" => method = new TRBDF2Parabolic1DMethod(payoff)
            case "EULER"  => method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
            case "LMG3"   => method = new LRE3Parabolic1DMethod(payoff)
            case "SCN"    => method = new SmoothCrankNicolsonCentralBSM1FMethod(payoff.knockouts.times, payoff)
          }
          //     method = new ThetaParabolic1DMethod(ThetaParabolic1DMethod.THETA_IMPLICIT)
          method.lowerBoundary = ForwardLinearOrder1Parabolic1DBoundaryFactory
          method.upperBoundary = BackwardLinearOrder1Parabolic1DBoundaryFactory

          var solver = new ThomasTridiagonalSolver(payoff)
          var pricer = new FDMSolver1D(spec, method, solver)
          var rtTime = System.nanoTime();
          pricer.solve(payoff)
          val price = pricer.price(spot)
          val time = (System.nanoTime() - rtTime) * 1e-9

          println(scheme + " " + timeSize + " " + timeMesh.size + " " + spaceSize + " " + price + " " + time)
        }
      }
    }
    assert(true)
  }

}
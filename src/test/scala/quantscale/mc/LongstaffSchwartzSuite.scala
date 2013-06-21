package test.quantscale.mc

import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.fdm.mesh.UniformMesh2D
import quantscale.fdm.ThomasTridiagonalSolver
import quantscale.fdm.payoff.SimpsonIntegralSmoother
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.fdm.method.TRBDF2Parabolic1DMethod
import quantscale.fdm.FDMSolver1D
import quantscale.fdm.ConstantBSM1FFDSpec
import quantscale.fdm.payoff.VanillaBermudanFDPayoff
import quantscale.fdm.mesh.MeshBoundaries
import quantscale.mc.ConstantBSM1FMCSpec
import quantscale.mc.BSMPathGenerator
import quantscale.mc.ConstantDiscounter
import quantscale.mc.BSMPathPricer
import quantscale.mc.MonteCarloSimulation
import quantscale.random.SobolFactory
import quantscale.mc.payoff.VanillaBermudanMCPayoff
import quantscale.mc.BSMMCSpec
import quantscale.mc.BrownianBridgeFactory
import quantscale.mc.LongstaffSchwartzPathPricer
import quantscale.mc.LongstaffSchwartzSimulation
import quantscale.random.PseudoRandomSequenceFactory
import quantscale.random.Well1024a
import quantscale.mc.GlassermanYuSimulation
import quantscale.mc.payoff.VanillaBermudanGYMCPayoff
import quantscale.mc.payoff.VanillaBasketBermudanGYMCPayoff
import quantscale.mc.payoff.VanillaBasketBermudanMCPayoff
import quantscale.mc.BSM1FMCSpec
import quantscale.random.MersenneTwister19937x64
import quantscale.mc.BrownianTransformFactory
import quantscale.random.LecuyerMRG32k3a
import quantscale.random.RandomSequenceFactory
import quantscale.random.Well19937

@RunWith(classOf[JUnitRunner])
class LongstaffSchwartzSuite extends FunSuite {
  test("bsm-bermudan") {
    val isCall = false
    val tte = 1.0;
    var spaceSize: Int = 800;
    var timeSize: Int = 100;
    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    val exerciseTimes = Array[Double](tte / 4.0, tte / 2.0, tte * 0.75, tte)

    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 3),
      strike);

    val fdspec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r)

    val fdpayoff = new VanillaBermudanFDPayoff(isCall, strike, exerciseTimes);
    val smoother = new SimpsonIntegralSmoother(Array(strike));
    val solver = new ThomasTridiagonalSolver(fdpayoff);
    val method = new TRBDF2Parabolic1DMethod(fdpayoff);
    val pricer = new FDMSolver1D(fdspec, method, solver);
    pricer.smoother = smoother;

    pricer.solve(fdpayoff);
    val priceAnal = BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
      strike,
      spot,
      BlackScholesVanillaEuropean.variance(vol, tte),
      BlackScholesVanillaEuropean.df(mu, tte),
      BlackScholesVanillaEuropean.df(r, tte));

    val price = pricer.price(spot);
    val refPrice = 13.386303
    println("priceAnal=" + priceAnal + " priceRef=" + refPrice + " price=" + price);

    //    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - priceAnal));

    val payoff = new VanillaBermudanMCPayoff(isCall, strike, exerciseTimes, exerciseTimes)
    val discounter = new ConstantDiscounter(r)
    val spec = new ConstantBSM1FMCSpec(spot, mu, vol)
    val rngFactory =
      SobolFactory
    //              new PseudoRandomSequenceFactory(
    //        new MersenneTwister19937x64())
    //                new Well1024a())
    //        new LecuyerMRG32k3a())
    val pathGen = new BSMPathGenerator(new BSMMCSpec(spec), rngFactory, null)
    var pathPricer = new LongstaffSchwartzPathPricer(payoff, discounter)
    val simNow = new LongstaffSchwartzSimulation(pathGen, pathPricer)
    val priceMc = simNow.evaluateForwardBackward(1024 * 16 + 1, 1024 * 32 * 4 + 1)
    println(priceMc)
    val payoffGY = new VanillaBermudanGYMCPayoff(isCall, strike, exerciseTimes, exerciseTimes)
    pathGen.pathListener = payoffGY.martingalePower.martingalePower
    pathPricer = new LongstaffSchwartzPathPricer(payoffGY, discounter)
    val simGY = new GlassermanYuSimulation(pathGen, pathPricer)
    val priceGY = simGY.evaluateForwardBackward(1024 * 16 + 1, 1024 * 32 * 4 + 1)
    println(priceGY)
    assert(Math.abs(priceMc - price) < 1e-2, "Monte-Carlo price " + priceMc + " too far off reference " + price)
  }

  test("bsm-bermudan-table") {
    val isCall = false
    val tte = 1.0;
    var spaceSize: Int = 4000;
    var timeSize: Int = 800;
    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    val exerciseTimes = Array[Double](tte / 4.0, tte / 2.0, tte * 0.75, tte)
    val nTrainingList = Array(512-1,1024-1, 1024 * 4 -1, 1024 * 16 -1, 1024 * 32 -1, 1024 * 64-1, 1024 * 128-1 )
//    val nSimList = Array(1024 * 1024 + 1)
    val nSimList = Array(1024 * 16-1, 1024 * 32-1 , 1024 * 64 -1, 1024 * 128-1 , 1024 * 256-1 , 1024 * 1024-1 , 1024*2048-1)
    val rngList = Array("Sobol", "Sobol-BB", "MT19937x64", "Well19937", "MRG32k3a")
    val grid: UniformMesh2D = new UniformMesh2D(
      spaceSize,
      timeSize,
      MeshBoundaries.makeBoundariesWithStdDev(tte, spot, vol, BlackScholesVanillaEuropean.df(mu, tte), 4),
      strike);

    val fdspec = new ConstantBSM1FFDSpec(
      grid,
      vol,
      mu,
      r)

    val fdpayoff = new VanillaBermudanFDPayoff(isCall, strike, exerciseTimes);
    val smoother = new SimpsonIntegralSmoother(Array(strike));
    val solver = new ThomasTridiagonalSolver(fdpayoff);
    val method = new TRBDF2Parabolic1DMethod(fdpayoff);
    val pricer = new FDMSolver1D(fdspec, method, solver);
    pricer.smoother = smoother;

    pricer.solve(fdpayoff);
    val priceAnal = BlackScholesVanillaEuropean.priceEuropeanVanilla(false,
      strike,
      spot,
      BlackScholesVanillaEuropean.variance(vol, tte),
      BlackScholesVanillaEuropean.df(mu, tte),
      BlackScholesVanillaEuropean.df(r, tte));

    val price = pricer.price(spot);
    val refPrice = 13.386303
    println("priceAnal=" + priceAnal + " priceRef=" + refPrice + " price=" + price);

    //    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - priceAnal));
    println("Method\tTrainingPaths\tSimulationPaths\tRNG\tTrainingValue\tValue\tTrainingError\tError")
    for (nTraining <- nTrainingList) {
      for (nSim <- nSimList) {
        for (rng <- rngList) {
          val payoff = new VanillaBermudanMCPayoff(isCall, strike, exerciseTimes, exerciseTimes)
          val discounter = new ConstantDiscounter(r)
          val spec = new ConstantBSM1FMCSpec(spot, mu, vol)
          val (rngFactory, transformFactory) = createRNGFactory(rng)
          val pathGen = new BSMPathGenerator(new BSMMCSpec(spec), rngFactory, transformFactory)
          var pathPricer = new LongstaffSchwartzPathPricer(payoff, discounter)
          var simNow = new LongstaffSchwartzSimulation(pathGen, pathPricer)
          var priceMc = simNow.evaluateForwardBackward(nTraining, nSim)
          println("LS\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simNow.trainingMean + "\t" + priceMc+ "\t"+(simNow.trainingMean-price)+"\t"+ (priceMc-price))
          val payoffGY = new VanillaBermudanGYMCPayoff(isCall, strike, exerciseTimes, exerciseTimes)
          pathGen.pathListener = payoffGY.martingalePower.martingalePower

//          pathPricer = new LongstaffSchwartzPathPricer(payoffGY, discounter)
//           simNow = new LongstaffSchwartzSimulation(pathGen, pathPricer)
//          priceMc = simNow.evaluateForwardBackward(nTraining, nSim)
//          println("LS-M\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simNow.trainingMean + "\t" + priceMc+ "\t"+(simNow.trainingMean-price)+"\t"+ (priceMc-price))
//          
          pathPricer = new LongstaffSchwartzPathPricer(payoffGY, discounter)
          var simGY = new GlassermanYuSimulation(pathGen, pathPricer)
          simGY.isITMRegression = true
          var priceGY = simGY.evaluateForwardBackward(nTraining, nSim)
          println("GY-Low-In\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.lowMean+"\t"+(simGY.trainingMean-price)+"\t"+ (simGY.lowMean-price))
          println("GY-Up-In\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.upMean+"\t"+ (simGY.trainingMean-price)+"\t"+ (simGY.upMean-price))
                    simGY = new GlassermanYuSimulation(pathGen, pathPricer)
          simGY.isITMRegression = false
          priceGY = simGY.evaluateForwardBackward(nTraining, nSim)
          println("GY-Low\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.lowMean+ "\t"+(simGY.trainingMean-price)+"\t"+ (simGY.lowMean-price))
          println("GY-Up\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.upMean+ "\t"+ (simGY.trainingMean-price)+"\t"+ (simGY.upMean-price))
        }
      }
    }
  }

  def createRNGFactory(rng: String): (RandomSequenceFactory, BrownianTransformFactory) = {
    rng match {
      case "Sobol"      => return (SobolFactory, null)
      case "Sobol-BB"   => return (SobolFactory, BrownianBridgeFactory)
      case "MT19937x64" => return (new PseudoRandomSequenceFactory(new MersenneTwister19937x64()), null)
      case "Well1024a"  => return (new PseudoRandomSequenceFactory(new Well1024a()), null)
      case "MRG32k3a"   => return (new PseudoRandomSequenceFactory(new LecuyerMRG32k3a()), null)
      case "Well19937" => return (new PseudoRandomSequenceFactory(new Well19937()), null)
    }
  }

  test("bsm-basket-bermudan") {
    val isCall = false
    val tte = 1.0;
    var spot = Array(100.0, 100.0);
    var strike = 1.0;
    val vol = Array(0.4, 0.2);
    val mu = Array(0.05, 0.05);
    val r = 0.05
    val weights = Array(0.4, 0.6)
    val correlation = Array(Array(1.0, 0.5), Array(0.5, 1.0))
    val exerciseTimes = Array[Double](tte / 4.0, tte / 2.0, tte * 0.75, tte)

    val refPrice = 13.386303

    val payoff = new VanillaBasketBermudanMCPayoff(isCall, strike, weights, spot, exerciseTimes, exerciseTimes)
    val discounter = new ConstantDiscounter(r)
    val specs: Array[BSM1FMCSpec] = Array(new ConstantBSM1FMCSpec(spot(0), mu(0), vol(0)), new ConstantBSM1FMCSpec(spot(1), mu(1), vol(1)))

    val rngFactory =
      SobolFactory
    //              new PseudoRandomSequenceFactory(
    //        new MersenneTwister19937x64())
    //                new Well1024a())
    //        new LecuyerMRG32k3a())
    val pathGen = new BSMPathGenerator(new BSMMCSpec(specs, correlation), rngFactory, BrownianBridgeFactory)
    var pathPricer = new LongstaffSchwartzPathPricer(payoff, discounter)
    val simNow = new LongstaffSchwartzSimulation(pathGen, pathPricer)
    val priceMc = simNow.evaluateForwardBackward(1024 * 1024*4 - 1, 1024 * 1024 *8 - 1)
    println(priceMc+" "+(simNow.trainingMean*(1024 * 1024*4 - 1)+priceMc*(1024 * 1024 *8 - 1))/(1024 * 1024 *4 - 1+1024 * 1024*8 - 1))
    val payoffGY = new VanillaBasketBermudanGYMCPayoff(isCall, strike, weights, spot, exerciseTimes, exerciseTimes)
    pathGen.pathListener = payoffGY.martingalePower.martingalePower
    pathPricer = new LongstaffSchwartzPathPricer(payoffGY, discounter)
    val simGY = new GlassermanYuSimulation(pathGen, pathPricer)
    val priceGY = simGY.evaluateForwardBackward(1024 * 1024*4-1, 1024 * 1024*8 - 1)
    println(priceGY+" "+(simGY.trainingMean*(1024 * 1024*4 - 1)+simGY.lowMean*(1024 * 1024 *8 - 1))/(1024 * 1024 *4 - 1+1024 * 1024*8 - 1))
    assert(Math.abs(priceMc - priceGY) < 1e-2, "Monte-Carlo price " + priceMc + " too far off reference " + priceGY)
  }
  
  test("bsm-basket-bermudan-table") {
    val isCall = false
    val tte = 1.0;
    var spaceSize: Int = 4000;
    var timeSize: Int = 800;
    var spot = Array(100.0, 100.0);
    var strike = 1.0;
    val vol = Array(0.4, 0.2);
    val mu = Array(0.05, 0.05);
    val r = 0.05
    val weights = Array(0.4, 0.6)
    val correlation = Array(Array(1.0, 0.5), Array(0.5, 1.0))
    val exerciseTimes = Array[Double](tte / 4.0, tte / 2.0, tte * 0.75, tte)
    val nTrainingList = Array(512-1,1024-1, 1024 * 4 -1, 1024 * 16 -1, 1024 * 32 -1, 1024 * 64-1, 1024 * 128-1 )
//    val nSimList = Array(1024 * 1024 + 1)
    val nSimList = Array(1024 * 16-1, 1024 * 32-1 , 1024 * 64 -1, 1024 * 128-1 , 1024 * 256-1 , 1024 * 1024-1 , 1024*2048-1)
    val rngList = Array("Sobol", "Sobol-BB", "MT19937x64", "Well19937", "MRG32k3a")

      val specs: Array[BSM1FMCSpec] = Array(new ConstantBSM1FMCSpec(spot(0), mu(0), vol(0)), new ConstantBSM1FMCSpec(spot(1), mu(1), vol(1)))


     val payoff = new VanillaBasketBermudanMCPayoff(isCall, strike, weights, spot, exerciseTimes, exerciseTimes)
    val price = 0.07550

    //    assert(Math.abs(price - refPrice) < 1e-1, price + " too far from " + refPrice + " diff=" + Math.abs(price - priceAnal));
    println("Method\tTrainingPaths\tSimulationPaths\tRNG\tTrainingValue\tValue\tTrainingError\tError")
    for (nTraining <- nTrainingList) {
      for (nSim <- nSimList) {
        for (rng <- rngList) {
          val payoff = new VanillaBasketBermudanMCPayoff(isCall, strike, weights, spot, exerciseTimes, exerciseTimes)
          val discounter = new ConstantDiscounter(r)
          val (rngFactory, transformFactory) = createRNGFactory(rng)
          val pathGen = new BSMPathGenerator(new BSMMCSpec(specs, correlation), rngFactory, transformFactory)
          var pathPricer = new LongstaffSchwartzPathPricer(payoff, discounter)
          var simNow = new LongstaffSchwartzSimulation(pathGen, pathPricer)
          var priceMc = simNow.evaluateForwardBackward(nTraining, nSim)
          println("LS\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simNow.trainingMean + "\t" + priceMc+ "\t"+(simNow.trainingMean-price)+"\t"+ (priceMc-price))
          val payoffGY = new VanillaBasketBermudanGYMCPayoff(isCall, strike, weights, spot, exerciseTimes, exerciseTimes)
          pathGen.pathListener = payoffGY.martingalePower.martingalePower
          pathPricer = new LongstaffSchwartzPathPricer(payoffGY, discounter)
          var simGY = new GlassermanYuSimulation(pathGen, pathPricer)
          simGY.isITMRegression = true
          var priceGY = simGY.evaluateForwardBackward(nTraining, nSim)
          println("GY-Low-In\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.lowMean+"\t"+(simGY.trainingMean-price)+"\t"+ (simGY.lowMean-price))
          println("GY-Up-In\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.upMean+"\t"+ (simGY.trainingMean-price)+"\t"+ (simGY.upMean-price))
                    simGY = new GlassermanYuSimulation(pathGen, pathPricer)
          simGY.isITMRegression = false
          priceGY = simGY.evaluateForwardBackward(nTraining, nSim)
          println("GY-Low\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.lowMean+ "\t"+(simGY.trainingMean-price)+"\t"+ (simGY.lowMean-price))
          println("GY-Up\t" + nTraining + "\t" + nSim + "\t" + rng + "\t" + simGY.trainingMean + "\t" + simGY.upMean+ "\t"+ (simGY.trainingMean-price)+"\t"+ (simGY.upMean-price))
        }
      }
    }
  }

}
package quantscale.mc

import org.scalatest.FunSuite
import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import quantscale.analytic._
import quantscale.mc.payoff.{VanillaBasketMCPayoff, VanillaMCPayoff}
import quantscale.random.SobolFactory
import quantscale.analytic.Price
import quantscale.analytic.DeltaGreek

@RunWith(classOf[JUnitRunner])
class AdjointMCSuite extends FunSuite {
  test("bsm-european") {
    val tte = 2.0;
    var spot = 100.0;
    var strike = 90.0;
    val vol = 0.3;
    val mu = 0.05;
    val r = 0.05;
    val isCall = true
    val measures: Array[BSMMeasure] = Array(new Price(), new DeltaGreek())
    BlackScholesVanillaEuropean.priceEuropeanVanillaAdjoint(isCall,
      strike,
      spot,
      new Variance(vol, tte),
      new DiscountFactor(mu, tte),
      new DiscountFactor(r, tte),
      measures)

    val weights = Array[Double](1.0)

    val payoff = new VanillaBasketMCPayoff(isCall, strike, weights, tte, tte)
    val discounter = new ConstantDiscounter(r)
    val rngFactory = SobolFactory
    for (i <- 0 until 10) {
      var startTime = System.nanoTime()
      var spec = new ConstantBSM1FMCSpec(spot, mu, vol)
      var pathGen = new BSMPathGenerator(new BSMMCSpec(spec), rngFactory, BrownianBridgeFactory)
      var pathPricer = new BSMPathPricer(payoff, discounter)
      var sim = new MonteCarloSimulation(pathGen, pathPricer)
      val numSimulations: Int = 1024 * 64
      val (priceMc, deltaMc) = sim.evaluateAdjointForward(numSimulations)
      var endTime = System.nanoTime()
      val priceAnal = measures(0).value
      println(priceAnal)
      println(priceMc)
      assert(Math.abs(priceMc - priceAnal) < 1e-2, "Monte Carlo price " + priceMc + " too far of reference price " + priceAnal)
      val deltaAnal = measures(1).value
      println(deltaAnal)
      println(deltaMc(0))
      println("adjoint delta error = " + (deltaMc(0) - deltaAnal) + " in " + (endTime - startTime) * 1e-9)
      startTime = System.nanoTime()
      spec = new ConstantBSM1FMCSpec(spot, mu, vol)
      pathGen = new BSMPathGenerator(new BSMMCSpec(spec), rngFactory, BrownianBridgeFactory)
      pathPricer = new BSMPathPricer(payoff, discounter)
      sim = new MonteCarloSimulation(pathGen, pathPricer)
      val priceMcRef = sim.evaluateForward(numSimulations)
      val eps = 1e-3
      val spotUp = spot + eps
      spec = new ConstantBSM1FMCSpec(spotUp, mu, vol)
      pathGen = new BSMPathGenerator(new BSMMCSpec(spec), rngFactory, BrownianBridgeFactory)
      pathPricer = new BSMPathPricer(payoff, discounter)
      sim = new MonteCarloSimulation(pathGen, pathPricer)
      val priceMcUp = sim.evaluateForward(numSimulations)
      val deltaFD = (priceMcUp - priceMcRef) / eps
      endTime = System.nanoTime()
      println(deltaFD)
      println("fd delta error = " + (deltaFD - deltaAnal) + " in " + (endTime - startTime) * 1e-9)
    }
  }

  test("bsm-basket") {
    val tte = 1.0;
    var spot = 100.0;
    var strike = 100.0;
    val vol = Array(0.3,0.2,0.25,0.4);
    val mu = Array(0,0.01,0.02,0.05);
    var correlation = Array(Array(1.0, 0.1, 0.1,0.1),
                            Array(0.1,1.0,0.1,0.1),
                            Array(0.1,0.1,1.0,0.1),
                            Array(0.1,0.1,0.1,1.0))

    val r = 0.01;
    val isCall = true

    val weights = Array[Double](0.3,0.2,0.25,0.25)

    val payoff = new VanillaBasketMCPayoff(isCall, strike, weights, tte, tte)
    val discounter = new ConstantDiscounter(r)
    val rngFactory = SobolFactory
    for (i <- 0 until 10) {
      var iAsset = 0
      var startTime = System.nanoTime()
      var spec = Array.ofDim[BSM1FMCSpec](weights.length)
      while (iAsset < weights.length) {
       spec(iAsset) = new ConstantBSM1FMCSpec(spot, mu(iAsset), vol(iAsset))
        iAsset +=1
      }
      var pathGen = new BSMPathGenerator(new BSMMCSpec(spec, correlation), rngFactory, BrownianBridgeFactory)
      var pathPricer = new BSMPathPricer(payoff, discounter)
      var sim = new MonteCarloSimulation(pathGen, pathPricer)
      val numSimulations: Int = 1024 * 64
      val (priceMc, deltaMc) = sim.evaluateAdjointForward(numSimulations)
      var endTime = System.nanoTime()
      println("adjoint "+priceMc +" "+deltaMc.mkString("[",",","]")+ " in " + (endTime - startTime) * 1e-9)
      startTime = System.nanoTime()
      iAsset = 0
      pathGen = new BSMPathGenerator(new BSMMCSpec(spec, correlation), rngFactory, BrownianBridgeFactory)
      pathPricer = new BSMPathPricer(payoff, discounter)
      sim = new MonteCarloSimulation(pathGen, pathPricer)
      val priceMcRef= sim.evaluateForward(numSimulations)
      val eps = 1e-3
      var deltaFD = Array.ofDim[Double](weights.length)
      while (iAsset < weights.length) {
        val oldSpec = spec(iAsset)
           spec(iAsset) = new ConstantBSM1FMCSpec(spot+eps, mu(iAsset), vol(iAsset))

       pathGen = new BSMPathGenerator(new BSMMCSpec(spec, correlation), rngFactory, BrownianBridgeFactory)
       pathPricer = new BSMPathPricer(payoff, discounter)
       sim = new MonteCarloSimulation(pathGen, pathPricer)
       val priceMcUp= sim.evaluateForward(numSimulations)
        deltaFD(iAsset) = (priceMcUp-priceMcRef)/eps;
        spec(iAsset) = oldSpec
        iAsset += 1
      }
       endTime = System.nanoTime()
      println("fd "+priceMcRef +" "+deltaFD.mkString("[",",","]")+ " in " + (endTime - startTime) * 1e-9)

    }
  }
}

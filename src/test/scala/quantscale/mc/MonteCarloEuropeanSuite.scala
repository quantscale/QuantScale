package test.quantscale.mc

import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.analytic.BlackScholesVanillaEuropean
import quantscale.mc.BSM1FMCSpec
import quantscale.mc.BSMMCSpec
import quantscale.mc.BSMPathGenerator
import quantscale.mc.BSMPathPricer
import quantscale.mc.BrownianBridgeFactory
import quantscale.mc.ConstantBSM1FMCSpec
import quantscale.mc.ConstantDiscounter
import quantscale.mc.MonteCarloSimulation
import quantscale.mc.payoff.VanillaMCPayoff
import quantscale.random.LecuyerMRG32k3a
import quantscale.random.MersenneTwister19937x64
import quantscale.random.PseudoRandomSequenceFactory
import quantscale.random.Sobol
import quantscale.random.SobolFactory
import quantscale.random.SobolFactory
import quantscale.random.Well1024a
import java.util.Arrays

@RunWith(classOf[JUnitRunner])
class MonteCarloEuropeanSuite extends FunSuite {

  test("lecuyer-uniform") {
    var average = 0.0
    val rng = new LecuyerMRG32k3a
    val n = 100000
    for (i <- 0 to n) {
      average += rng.nextDouble
    }
    average /= n
    assert(math.abs(average - 0.5) < 1e-3, "average is " + average)
  }

  test("sobol-first-numbers") {
    var ref = Array(Array(0.5, 0.5, 0.5),
      Array(0.75, 0.25, 0.25),
      Array(0.25, 0.75, 0.75),
      Array(0.375, 0.375, 0.625),
      Array(0.875, 0.875, 0.125),
      Array(0.625, 0.125, 0.875),
      Array(0.125, 0.625, 0.375),
      Array(0.1875, 0.3125, 0.9375),
      Array(0.6875, 0.8125, 0.4375))

    var rng = SobolFactory.makeRandomSequenceGenerator(3)
    var i = 0
    var points = Array.ofDim[Double](3)
    while (i < 9) {
      rng.nextSequence(points)
      println(Arrays.toString(points))
      for (j <- 0 until 3) {
        assert(points(j) == ref(i)(j), "points for sequence " + i + " not equals: expected " + ref(i)(j) + " but was " + points(j))
      }

      i += 1
    }
  }

  test("bsm-european") {
    val tte = 1.0;
    var spot = 100.0;
    var strike = 100.0;
    val vol = 0.4;
    val mu = 0.05;
    val r = 0.05;
    val isCall = true
    val priceAnal = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall,
      strike,
      spot,
      BlackScholesVanillaEuropean.variance(vol, tte),
      BlackScholesVanillaEuropean.df(mu, tte),
      BlackScholesVanillaEuropean.df(r, tte));

    val payoff = new VanillaMCPayoff(isCall, strike, tte, tte)
    val discounter = new ConstantDiscounter(r)
    val spec = new ConstantBSM1FMCSpec(spot, mu, vol)
    val rngFactory =
      SobolFactory
    //      new PseudoRandomSequenceFactory(
    //        new MersenneTwister19937x64())
    //        new Well1024a())
    //        new LecuyerMRG32k3a())
    val pathGen = new BSMPathGenerator(new BSMMCSpec(spec), rngFactory, BrownianBridgeFactory)
    val pathPricer = new BSMPathPricer(payoff, discounter)
    val sim = new MonteCarloSimulation(pathGen, pathPricer)
    val priceMc = sim.evaluateForward(1024 * 32)
    println(priceAnal)
    println(priceMc)
    assert(Math.abs(priceMc - priceAnal) < 1e-2, "Monte Carlo price " + priceMc + " too far of reference price " + priceAnal)
  }
}
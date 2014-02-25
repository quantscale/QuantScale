package quantscale.analytic

import org.scalatest.FunSuite

class BarrierAdjointSuite extends FunSuite {
  test("adjoint-vs-fd") {
    val isCall = true
    val strike = 99.0
    val barrier = 121.0
    var spot = 100.0
    val tte=0.7
    val vol = 0.31
    val growthRate = 0.01
    val discountRate = 0.01
    val values = BlackScholesMertonBarrier.priceUpAndOutAdjoint(isCall, strike, barrier, spot, tte, vol, growthRate, discountRate)
    val eps = 1e-7
    val vm0 = (BlackScholesMertonBarrier.priceUpAndOutAdjoint(isCall, strike, barrier, spot+eps, tte, vol, growthRate, discountRate)(0)-values(0))/eps
    val vm1 = (BlackScholesMertonBarrier.priceUpAndOutAdjoint(isCall, strike, barrier, spot, tte, vol+eps, growthRate, discountRate)(0)-values(0))/eps
    val vm2 = (BlackScholesMertonBarrier.priceUpAndOutAdjoint(isCall, strike, barrier, spot, tte, vol, growthRate+eps, discountRate)(0)-values(0))/eps
    val vm3 = (BlackScholesMertonBarrier.priceUpAndOutAdjoint(isCall, strike, barrier, spot, tte, vol, growthRate, discountRate+eps)(0)-values(0))/eps
    val vm4 = (BlackScholesMertonBarrier.priceUpAndOutAdjoint(isCall, strike, barrier, spot, tte+eps, vol, growthRate, discountRate)(0)-values(0))/eps
    println(values(1), vm0)
    println(values(2), vm1)
    println(values(3), vm2)
    println(values(4), vm3)
    println(values(5), vm4)
  }
}

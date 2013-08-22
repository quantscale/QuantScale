package quantscale.fdm.sabr

import quantscale.analytic.SABRModelSpec
import quantscale.fdm.TridiagonalMatrix

class HaganRichardsonEulerNSABRTransformedDensitySolver(spec: SABRModelSpec, forward: Double, T: Double, size: Int = 100, timeSteps: Int = 10, nDeviations: Double = 3.0)
  extends HaganSABRTransformedDensitySolver(spec, forward, T, size, timeSteps, nDeviations) {

  var order = 5
  private var EmCache: Array[Array[Double]] = null

  override def solve() {
    EmCache = Array.ofDim(order, size)
    buildEmCache()
    val EmInit = Em_.clone()
    var io = 0
    var Po: Array[Double] = null
    var PLo, PRo = 0.0
    while (io < order) {
      Array.copy(EmInit, 0, Em_, 0, size)
      solve(io + 1)
      if (io == 0) {
        Po = P0_
        PLo = PL_
        PRo = PR_
      } else {
        val pow = math.pow(2, io)
        var j = size - 2
        while (j > 0) {
          Po(j) = (pow * P0_(j) - Po(j)) / (pow - 1)
          j -= 1
        }
        PLo = (pow * PL_ - PLo) / (pow - 1)
        PRo = (pow * PR_ - PRo) / (pow - 1)
      }
      io += 1
    }

    P0_ = Po
    PL_ = PLo
    PR_ = PRo
  }


  private[sabr] def buildEmCache() {
    buildEm()
    var io = order - 1
    val dto = dt / math.pow(2, io) //EmCache(i) = cache for timestep 2^i
    val EmCacheo = Array.ofDim[Double](size)
    var j = size - 2
    val efactor = spec.rho * spec.nu * spec.alpha
    while (j > 0) {
      EmCacheo(j) = math.exp(efactor * Gammam_(j) * dto)
      j -= 1
    }
    EmCacheo(0) = EmCacheo(1)
    EmCacheo(size - 1) = EmCacheo(size - 2)
    EmCache(io) = EmCacheo
    io -= 1
    while (io >= 0) {
      EmCache(io) = Array.ofDim[Double](size)
      var j = size - 1
      val EmCacheio = EmCache(io)
      val EmCacheio1 = EmCache(io + 1)

      while (j >= 0) {
        EmCacheio(j) = EmCacheio1(j) * EmCacheio1(j)
        j -= 1
      }
      io -= 1
    }
  }

  private[sabr] def advanceEm(order: Int, M: Array[Double]) {
    var j = 0
    val Mcache = EmCache(order - 1)

    while (j < size) {
      M(j) = M(j) * Mcache(j)
      j += 1
    }
  }

  def solve(order: Int) {
    val divisor = math.pow(2, order - 1)
    val dto = dt / divisor
    tri1 = new TridiagonalMatrix(size)
    P0_ = computeP()
    P1_ = Array.ofDim(size)
    PL_ = 0.0
    PR_ = 0.0
    var t = T
    var tIndex = 0
    while (tIndex < timeSteps * divisor) {
      t -= dto
      advanceEm(order, Em_)
      computeSystem(dto, Em_, tri1)
      P0_(0) = 0
      P0_(size - 1) = 0
      solver.solve(tri1, P0_, P1_)
      PL_ += dto * computedPLdt(Em_, P1_)
      PR_ += dto * computedPRdt(Em_, P1_)
      //   printSumQF("RE",t,Q1, QL_, QR_)
      val Qtmp = P0_
      P0_ = P1_
      P1_ = Qtmp
      tIndex += 1
    }
  }


}
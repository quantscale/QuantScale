package quantscale.analytic

import quantscale.fdm.payoff.AccumulatorFDPayoff
import scala.collection.mutable.ArrayBuffer
import quantscale.fdm.Epsilon
import java.sql.Time
import quantscale.math.CumulativeBivariateNormalDistribution
import quantscale.math.Function1D

object BarrierType extends Enumeration {
  type BarrierType = Value
  val UP_AND_OUT, DOWN_AND_OUT, UP_AND_IN, DOWN_AND_IN = Value
}
import quantscale.analytic.BarrierType._

class AnalyticBarrierPayoff(val isCall: Boolean, val strike: Double, val barrier: Double, val barrierType: BarrierType, val rebate: Double = 0.0) {

}

object BlackScholesMertonBarrier {
  def priceOutPartialEnd(isCall: Boolean, isUp: Boolean,
                         K: Double,
                         B: Double,
                         S: Double,
                         vol1: Double,
                         vol2b: Double, //black vol
                         r: Double,
                         q: Double,
                         t1: Double,
                         T2: Double,
                         bvd: CumulativeBivariateNormalDistribution,
                         cnd: Function1D): Double = {
    val p = new Price
    priceOutPartialEnd(isCall, isUp, K, B, S, vol1, vol2b, r, q, t1, T2, Array[BSMMeasure](p), bvd, cnd)
    return p.value
  }

  private def dM(cnd: quantscale.math.Function1D, us: Int, e1: Double, g1: Double, rho: Double, de1: Double, dg1: Double): Double = {
    return us * (de1 * NormalDistribution.value(e1) * cnd.value((us * g1 - rho * us * e1) / math.sqrt(1 - rho * rho)) +
      dg1 * NormalDistribution.value(g1) * cnd.value((us * e1 - rho * us * g1) / math.sqrt(1 - rho * rho)))
  }

  private def d2M(cnd: quantscale.math.Function1D, us: Int, e1: Double, g1: Double, rho: Double, de1: Double, dg1: Double, d2e1: Double, d2g1: Double): Double = {
    return us * ((d2e1 - de1 * de1 * e1) * NormalDistribution.value(e1) * cnd.value((us * g1 - rho * us * e1) / math.sqrt(1 - rho * rho)) +
      de1 * NormalDistribution.value(e1) * (us * dg1 - rho * us * de1) / math.sqrt(1 - rho * rho) * NormalDistribution.value((us * e1 - rho * us * g1) / math.sqrt(1 - rho * rho)) +
      (d2g1 - dg1 * dg1 * g1) * NormalDistribution.value(g1) * cnd.value((us * e1 - rho * us * g1) / math.sqrt(1 - rho * rho)) +
      dg1 * NormalDistribution.value(g1) * (us * de1 - rho * us * dg1) / math.sqrt(1 - rho * rho) * NormalDistribution.value((us * e1 - rho * us * g1) / math.sqrt(1 - rho * rho)))
  }

  def priceOutPartialEnd(isCall: Boolean, isUp: Boolean,
                         K: Double,
                         B: Double,
                         S: Double,
                         vol1: Double,
                         vol2b: Double, //black v, 
                         r: Double,
                         q: Double,
                         t1: Double,
                         T2: Double,
                         measures: Array[BSMMeasure],
                         bvd: CumulativeBivariateNormalDistribution,
                         cnd: Function1D = CodyCND): Array[Double] = {
    val cs = if (isCall) 1 else -1
    val us = if (isUp) 1 else -1
    val FT2 = S * math.exp((r - q) * T2)
    val dfT2 = math.exp(-r * T2)
    val vol2 = math.sqrt((vol2b * vol2b * T2 - vol1 * vol1 * t1) / (T2 - t1)) //forward vol
    val varT2 = vol2 * vol2 * (T2 - t1) + vol1 * vol1 * t1
    val sqrtVarT2 = math.sqrt(varT2)
    val vart1 = vol1 * vol1 * t1
    val sqrtVart1 = math.sqrt(vart1)
    val eterm = ((r - q) * vol1 * vol1 / (vol2 * vol2) - (r - q)) * t1
    val d1 = (math.log(S / K) + (r - q) * T2 + 0.5 * varT2) / sqrtVarT2
    val d2 = d1 - sqrtVarT2
    val f1 = d1 + 2*( math.log(B / S) + eterm) / sqrtVarT2
    val f2 = f1 - sqrtVarT2
    val e1 = (math.log(S / B) + (r - q) * t1 + 0.5 * vart1) / sqrtVart1
    val e2 = e1 - sqrtVart1
    val e3 = e1 + 2 * (math.log(B / S) + eterm) / sqrtVart1
    val e4 = e3 - sqrtVart1
    val g1 = (math.log(S / B) + (r - q) * T2 + 0.5 * varT2) / sqrtVarT2
    val g2 = g1 - sqrtVarT2
    val g3 = g1 + 2 * (math.log(B / S) + eterm) / sqrtVarT2
    val g4 = g3 - sqrtVarT2
    val mu = ((r - q) - vol2 * vol2 * 0.5) / (vol2 * vol2)
    val rho = math.sqrt(vart1 / varT2)
    val exp = math.exp(eterm)
    val Z0 = bvd.value(-us * e1, -us * g1, rho, cnd)
    val P1 = math.pow(B / S * exp, 2 * (mu + 1))
    val M1 = bvd.value(us * e3, -us * g3, -rho, cnd)
    val Z2 = bvd.value(-us * e2, -us * g2, rho, cnd)
    val P3 = math.pow(B / S * exp, 2 * mu)
    val M3 = bvd.value(us * e4, -us * g4, -rho, cnd)
    val Z4 = bvd.value(-us * e1, -us * d1, rho, cnd)
    val M5 = bvd.value(-us * f1, us * e3, -rho, cnd)
    val Z6 = bvd.value(-us * e2, -us * d2, rho, cnd)
    val M7 = bvd.value(-us * f2, us * e4, -rho, cnd)

    val cA = us * FT2 * dfT2 * (Z0 - P1 * M1) - us * K * dfT2 * (Z2 - P3 * M3)
    val cB = -cs * FT2 * dfT2 * 0.5 * (1 + cs * us) * (Z4 - P1 * M5) + cs * 0.5 * (1 + cs * us) * K * dfT2 * (Z6 - P3 * M7)

    for (measure <- measures) {
      measure match {
        case Price() => measure.value = cA + cB
        case DeltaGreek() => {
          val de1 = 1.0 / (sqrtVart1 * S)
          val dg1 = 1.0 / (sqrtVarT2 * S)
          val de3 = -de1
          val dg3 = -dg1
          val de2 = de1;
          val dg2 = dg1;
          val de4 = de3;
          val dg4 = dg3;
          val dd1 = dg1;
          val df1 = -dg1
          val dd2 = dd1;
          val df2 = df1
          val dZ0 = dM(cnd, us, -e1, -g1, rho, -de1, -dg1)
          val d2e1 = -de1 / S
          val d2g1 = -dg1 / S
          val dM1 = dM(cnd, us, e3, -g3, -rho, de3, -dg3)
          val dZ2 = dM(cnd, us, -e2, -g2, rho, -de2, -dg2)
          val dM3 = dM(cnd, us, e4, -g4, -rho, de4, -dg4)
          val dZ4 = dM(cnd, us, -e1, -d1, rho, -de1, -dd1)
          val dM5 = dM(cnd, us, -f1, e3, -rho, -df1, de3)
          val dZ6 = dM(cnd, us, -e2, -d2, rho, -de2, -dd2)
          val dM7 = dM(cnd, us, -f2, e4, -rho, -df2, de4)

          val dP1 = -(2 * (mu + 1)) * P1 / S
          val dP3 = -2 * mu * P3 / S

          //          println("an details " + dZ0 + " " + dP1 + " " + dM1 + " " + dZ2 + " " + dP3 + " " + dM3 + " " + dZ4 + " " + dM5 + " " + dZ6 + " " + dM7)
          val deltaA = us * FT2 / S * dfT2 * (Z0 - P1 * M1) +
            us * FT2 * dfT2 * (dZ0 - dP1 * M1 - P1 * dM1) -
            us * K * dfT2 * (dZ2 - dP3 * M3 - P3 * dM3)
          val deltaB = -cs * FT2 / S * dfT2 * 0.5 * (1 + cs * us) * (Z4 - P1 * M5) -
            cs * FT2 * dfT2 * 0.5 * (1 + cs * us) * (dZ4 - dP1 * M5 - P1 * dM5) +
            cs * 0.5 * (1 + cs * us) * K * dfT2 * (dZ6 - dP3 * M7 - P3 * dM7)

          measure.value = deltaA + deltaB
        }
        case GammaGreek() => {
          val de1 = 1.0 / (sqrtVart1 * S)
          val dg1 = 1.0 / (sqrtVarT2 * S)
          val de3 = -de1
          val dg3 = -dg1
          val de2 = de1;
          val dg2 = dg1;
          val de4 = de3;
          val dg4 = dg3;
          val dd1 = dg1;
          val df1 = -dg1
          val dd2 = dd1;
          val df2 = df1
          val dZ0 = dM(cnd, us, -e1, -g1, rho, -de1, -dg1)
          val d2e1 = -de1 / S
          val d2g1 = -dg1 / S
          val d2e3 = -d2e1
          val d2g3 = -d2g1
          val d2e2 = d2e1
          val d2g2 = d2g1
          val d2e4 = d2e3
          val d2g4 = d2g3
          val d2d1 = d2g1
          val d2f1 = -d2g1
          val d2d2 = d2d1
          val d2f2 = d2f1
          val d2Z0 = d2M(cnd, us, -e1, -g1, rho, -de1, -dg1, -d2e1, -d2g1)
          val dM1 = dM(cnd, us, e3, -g3, -rho, de3, -dg3)
          val d2M1 = d2M(cnd, us, e3, -g3, -rho, de3, -dg3, d2e3, -d2g3)
          val dZ2 = dM(cnd, us, -e2, -g2, rho, -de2, -dg2)
          val d2Z2 = d2M(cnd, us, -e2, -g2, rho, -de2, -dg2, -d2e2, -d2g2)
          val dM3 = dM(cnd, us, e4, -g4, -rho, de4, -dg4)
          val d2M3 = d2M(cnd, us, e4, -g4, -rho, de4, -dg4, d2e4, -d2g4)
          val dZ4 = dM(cnd, us, -e1, -d1, rho, -de1, -dd1)
          val d2Z4 = d2M(cnd, us, -e1, -d1, rho, -de1, -dd1, -d2e1, -d2d1)
          val dM5 = dM(cnd, us, -f1, e3, -rho, -df1, de3)
          val d2M5 = d2M(cnd, us, -f1, e3, -rho, -df1, de3, -d2f1, d2e3)
          val dZ6 = dM(cnd, us, -e2, -d2, rho, -de2, -dd2)
          val d2Z6 = d2M(cnd, us, -e2, -d2, rho, -de2, -dd2, -d2e2, -d2d2)
          val dM7 = dM(cnd, us, -f2, e4, -rho, -df2, de4)
          val d2M7 = d2M(cnd, us, -f2, e4, -rho, -df2, de4, -d2f2, d2e4)

          val dP1 = -(2 * (mu + 1)) * P1 / S
          val d2P1 = -(2 * (mu + 1) + 1) * dP1 / S
          val dP3 = -2 * mu * P3 / S
          val d2P3 = -(2 * mu + 1) * dP3 / S

          //          println("an details " + dZ0 + " " + dP1 + " " + dM1 + " " + dZ2 + " " + dP3 + " " + dM3 + " " + dZ4 + " " + dM5 + " " + dZ6 + " " + dM7)
          //          println("an details " + d2Z0 + " " + d2P1 + " " + d2M1 + " " + d2Z2 + " " + d2P3 + " " + d2M3 + " " + d2Z4 + " " + d2M5 + " " + d2Z6 + " " + d2M7)
          val deltaB = -cs * FT2 / S * dfT2 * 0.5 * (1 + cs * us) * (Z4 - P1 * M5) -
            cs * FT2 * dfT2 * 0.5 * (1 + cs * us) * (dZ4 - dP1 * M5 - P1 * dM5) +
            cs * 0.5 * (1 + cs * us) * K * dfT2 * (dZ6 - dP3 * M7 - P3 * dM7)

          val gammaA = 2 * us * FT2 / S * dfT2 * (dZ0 - dP1 * M1 - P1 * dM1) +
            us * FT2 * dfT2 * (d2Z0 - d2P1 * M1 - 2 * dP1 * dM1 - P1 * d2M1) -
            us * K * dfT2 * (d2Z2 - d2P3 * M3 - 2 * dP3 * dM3 - P3 * d2M3)

          val gammaB = -2 * cs * FT2 / S * dfT2 * 0.5 * (1 + cs * us) * (dZ4 - dP1 * M5 - P1 * dM5) -
            cs * FT2 * dfT2 * 0.5 * (1 + cs * us) * (d2Z4 - d2P1 * M5 - 2 * dP1 * dM5 - P1 * d2M5) +
            cs * 0.5 * (1 + cs * us) * K * dfT2 * (d2Z6 - d2P3 * M7 - 2 * dP3 * dM7 - P3 * d2M7)
          measure.value = gammaA + gammaB
        }
      }
    }
    return Array(Z0, P1, M1, Z2, P3, M3, Z4, M5, Z6, M7)
  }

  def priceDownAndOut(isCall: Boolean, strike: Double, barrier: Double, spot: Double, variance: Double, driftDf: Double, discountDf: Double, driftDfDelivery: Double, discountDfDelivery: Double): Double = {
    val bs0 = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance, driftDf, discountDf)
    val bs1 = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, barrier * barrier / spot, variance, driftDf, discountDf)
    return bs0 - math.pow(barrier / spot, -1.0 - 2 * math.log(driftDf) / variance) * bs1
  }

  def priceUpAndOut2(isCall: Boolean, strike: Double, barrier: Double, spot: Double, variance: Double, driftDf: Double, discountDf: Double, driftDfDelivery: Double, discountDfDelivery: Double): Double = {
    val bs0 = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, spot, variance, driftDf, discountDf)
    val bsH = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, barrier, spot, variance, driftDf, discountDf)
    val bs1 = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, strike, barrier * barrier / spot, variance, driftDf, discountDf)
    val bs2 = BlackScholesVanillaEuropean.priceEuropeanVanilla(isCall, barrier, barrier * barrier / spot, variance, driftDf, discountDf)
    val sqrtVar = math.sqrt(variance)
    val dsh = (math.log(spot / barrier / driftDf) + variance * 0.5) / sqrtVar
    val dhs = (math.log(barrier / spot / driftDf) + variance * 0.5) / sqrtVar
    val ndsh = CumulativeNormalDistribution.value(dsh)
    val ndhs = CumulativeNormalDistribution.value(dhs)
    return bs0 - bsH - (barrier - strike) * discountDf * ndsh -
      math.pow(barrier / spot, -1.0 - 2 * math.log(driftDf) / variance) * (bs1 - bs2 - (barrier - strike) * discountDf * ndhs)
  }
  def priceUpAndOut(isCall: Boolean, strike: Double, barrier: Double, spot: Double, variance: Double, driftDf: Double, discountDf: Double, driftDfDelivery: Double, discountDfDelivery: Double): Double = {
    val sqrtVar = math.sqrt(variance)
    val x = (math.log(spot / strike / driftDf) + variance * 0.5) / sqrtVar
    val x1 = (math.log(spot / barrier / driftDf) + variance * 0.5) / sqrtVar
    val y = (math.log(barrier * barrier / (spot * strike) / driftDf) + 0.5 * variance) / sqrtVar
    val y1 = (math.log(barrier / spot / driftDf) + 0.5 * variance) / sqrtVar
    val lambda = -0.5 - (math.log(driftDf)) / variance

    val sign = if (isCall) 1 else -1
    val nx = CumulativeNormalDistribution.value(sign * x)
    val nxd = CumulativeNormalDistribution.value(sign * (x - sqrtVar))
    val nx1 = CumulativeNormalDistribution.value(x1)
    val nx1d = CumulativeNormalDistribution.value(x1 - sqrtVar)
    val ny = CumulativeNormalDistribution.value(-y)
    val ny1 = CumulativeNormalDistribution.value(-y1)
    val nyd = CumulativeNormalDistribution.value(-y + sqrtVar)
    val ny1d = CumulativeNormalDistribution.value(-y1 + sqrtVar)

    val undiscPrice = if (isCall) spot * driftDfDelivery * (nx - nx1 +
      math.pow(barrier / spot, 2 * lambda) * (ny - ny1)) -
      strike * (nxd - nx1d + math.pow(barrier / spot, 2 * lambda - 2) * (nyd - ny1d))
    else strike * (nxd - math.pow(barrier / spot, 2 * lambda - 2) * (nyd)) -
      spot * driftDfDelivery * (nx - math.pow(barrier / spot, 2 * lambda) * ny)
    return discountDfDelivery * undiscPrice
  }

  def priceAccumulator(buySell: Int, isKoAfterAccrual: Boolean,
                       periodAccrual: ArrayBuffer[Double],
                       knockouts: ArrayBuffer[Double],
                       payments: ArrayBuffer[Double],
                       strike: Double,
                       accrualAbove: Double,
                       accrualBelow: Double,
                       barrier: Double,
                       spot: Double,
                       vol: Double,
                       discountRate: Double,
                       driftRate: Double): Double = {

    //        val koFrequency = knockouts(knockouts.length-1) / knockouts.length
    //    val sqrtKoFrequency = math.sqrt(koFrequency)
    //    val ybar = computeYbar(math.log(spot), math.log(barrier), vol, sqrtKoFrequency)
    //    val shift = computeRelativeShift(sqrtKoFrequency, 0.0, vol, ybar)

    var sum = 0.0

    var i = 0
    val lastTp = payments(payments.length - 1)
    for (acc <- periodAccrual) {
      i += 1
      if (acc < Epsilon.MACHINE_EPSILON_SQRT) {
        //TODO barrier is hit
        sum += accrualAbove * math.max((spot - strike), 0) - accrualBelow * math.max(strike - spot, 0)
      } else {
        val tp = getPayment(payments, acc)
        val driftDf = math.exp(-driftRate * acc)
        val discountDf = math.exp(-discountRate * acc)
        val driftDfDelivery = math.exp(-driftRate * tp)
        val discountDfDelivery = math.exp(-discountRate * tp)
        val variance = vol * vol * acc
        //                val koFrequency = tp/ getNumberOfPayments(payments, tp)
        //                val sqrtKoFrequency = math.sqrt(koFrequency)
        //                val ybar = computeYbar(math.log(spot), math.log(barrier), vol, math.sqrt(tp))
        //                val shift = computeRelativeShift(sqrtKoFrequency, 0.0, vol)
        val koFrequency = acc / i
        val sqrtKoFrequency = math.sqrt(koFrequency)
        val ybar = computeYbar(math.log(spot), math.log(barrier), vol, sqrtKoFrequency)
        val shift = computeRelativeShift(sqrtKoFrequency, 0.0, vol, ybar)

        sum += accrualAbove * price(new AnalyticBarrierPayoff(true, strike, barrier * shift, UP_AND_OUT), spot, variance, driftDf, driftDfDelivery, discountDfDelivery, discountDfDelivery) -
          accrualBelow * price(new AnalyticBarrierPayoff(false, strike, barrier * shift, UP_AND_OUT), spot, variance, driftDf, driftDfDelivery, discountDfDelivery, discountDfDelivery)

        if (isKoAfterAccrual) {
          if (acc == lastTp) {
            sum += accrualAbove * price(new AnalyticBarrierPayoff(true, strike, barrier * shift, UP_AND_IN), spot, variance, driftDf, driftDfDelivery, discountDfDelivery, discountDfDelivery) -
              accrualBelow * price(new AnalyticBarrierPayoff(false, strike, barrier * shift, UP_AND_IN), spot, variance, driftDf, driftDfDelivery, discountDfDelivery, discountDfDelivery)
          }
        }
      }
    }

    return sum
  }

  private def getPayment(payments: ArrayBuffer[Double], d: Double): Double = {
    var i = 0
    while (i < payments.length) {
      if (payments(i) > d - Epsilon.MACHINE_EPSILON_SQRT) return payments(i)
      i += 1
    }
    throw new RuntimeException("Date not found " + d + " last date was " + payments(payments.length - 1))
  }

  def priceCallMinus2PutsKO(strike: Double, barrier: Double, spot: Double, variance: Double, driftDf: Double, discountDf: Double, driftDfPayment: Double, discountDfPayment: Double): Double = {
    val sqrtVar = math.sqrt(variance)
    val x = (math.log(spot / strike / driftDf) + variance * 0.5) / sqrtVar
    val x1 = (math.log(spot / barrier / driftDf) + variance * 0.5) / sqrtVar
    val y = (math.log(barrier * barrier / (spot * strike) / driftDf) + 0.5 * variance) / sqrtVar
    val y1 = (math.log(barrier / spot / driftDf) + 0.5 * variance) / sqrtVar
    val lambda = 0.5 - math.log(driftDf) / variance

    val nx = CumulativeNormalDistribution.value(x)
    val nxd = CumulativeNormalDistribution.value(x - sqrtVar)
    val nx1 = CumulativeNormalDistribution.value(x1)
    val nx1d = CumulativeNormalDistribution.value(x1 - sqrtVar)
    val ny = CumulativeNormalDistribution.value(-y)
    val ny1 = CumulativeNormalDistribution.value(-y1)
    val nyd = CumulativeNormalDistribution.value(-y + sqrtVar)
    val ny1d = CumulativeNormalDistribution.value(-y1 + sqrtVar)

    val undiscPrice = spot * driftDfPayment * (2 - nx - nx1 -
      math.pow(barrier / spot, 2 * lambda) * (ny + ny1)) -
      strike * (2 - nxd - nx1d + math.pow(barrier / spot, 2 * lambda - 2) * (nyd + ny1d))
    return discountDfPayment * undiscPrice
  }

  def computeRelativeShift(sqrtBarrierFrequency: Double, sqrtOriginalFrequency: Double, v: Double, ybar: Double = 0.5826): Double = {
    return math.exp(ybar * v * (sqrtBarrierFrequency - sqrtOriginalFrequency));
  }

  def computeYbar(logx0: Double, loglevel: Double, v0: Double, sqrtdt: Double): Double = {
    val u = math.abs(logx0 - loglevel) / (v0 * sqrtdt);
    return 0.5826 + 0.1245 * math.exp(-2.7 * math.pow(u, 1.2));
  }

  def price(payoff: AnalyticBarrierPayoff, spot: Double, variance: Double, driftDf: Double, driftDfDelivery: Double, discountDfRebate: Double, discountDfDelivery: Double): Double = {
    return new QuantlibBarrierAnalyticEngine(payoff, spot, variance, driftDf, driftDfDelivery, discountDfRebate, discountDfDelivery).calculate()
  }

  class QuantlibBarrierAnalyticEngine(val payoff: AnalyticBarrierPayoff, val underlying: Double, val variance: Double, val driftDf: Double, val driftDfDelivery: Double, val discountDfRebate: Double, val discountDfDelivery: Double) {
    private val sqrtVar = math.sqrt(variance)

    private val mu = -math.log(driftDf) / (variance) - 0.5;

    private val muSigma = (1 + mu) * sqrtVar;

    def calculate(): Double = {
      if (payoff.strike < 0.0) throw new RuntimeException("strike must be positive")
      if (underlying <= 0) throw new RuntimeException("negative or null underlying given")

      //        QL_REQUIRE(!triggered(spot), "barrier touched");

      var value = 0.
      if (payoff.isCall) {
        payoff.barrierType match {
          case DOWN_AND_IN =>
            if (payoff.strike >= payoff.barrier)
              value = C(1, 1) + E(1)
            else
              value = A(1) - B(1) + D(1, 1) + E(1);
          case UP_AND_IN =>
            if (payoff.strike >= payoff.barrier)
              value = A(1) + E(-1);
            else
              value = B(1) - C(-1, 1) + D(-1, 1) + E(-1);
          case DOWN_AND_OUT =>
            if (payoff.strike >= payoff.barrier)
              value = A(1) - C(1, 1) + F(1);
            else
              value = B(1) - D(1, 1) + F(1);
          case UP_AND_OUT =>
            if (payoff.strike >= payoff.barrier)
              value = F(-1);
            else
              value = A(1) - B(1) + C(-1, 1) - D(-1, 1) + F(-1);

        }
      } else {
        payoff.barrierType match {
          case DOWN_AND_IN =>
            if (payoff.strike >= payoff.barrier)
              value = B(-1) - C(1, -1) + D(1, -1) + E(1);
            else
              value = A(-1) + E(1);
          case UP_AND_IN =>
            if (payoff.strike >= payoff.barrier)
              value = A(-1) - B(-1) + D(-1, -1) + E(-1);
            else
              value = C(-1, -1) + E(-1);
          case DOWN_AND_OUT =>
            if (payoff.strike >= payoff.barrier)
              value = A(-1) - B(-1) + C(1, -1) - D(1, -1) + F(1);
            else
              value = F(1);
          case UP_AND_OUT =>
            if (payoff.strike >= payoff.barrier)
              value = B(-1) - D(-1, -1) + F(-1);
            else
              value = A(-1) - C(-1, -1) + F(-1);
        }
      }
      return value
    }

    private def A(phi: Int): Double = {
      val x1 =
        math.log(underlying / payoff.strike) / sqrtVar + muSigma;
      val N1 = CumulativeNormalDistribution.value(phi * x1);
      val N2 = CumulativeNormalDistribution.value(phi * (x1 - sqrtVar));
      return phi * discountDfDelivery * (underlying / driftDfDelivery * N1 - payoff.strike * N2);
    }

    private def B(phi: Int): Double = {
      val x2 =
        math.log(underlying / payoff.barrier) / sqrtVar + muSigma;
      val N1 = CumulativeNormalDistribution.value(phi * x2);
      val N2 = CumulativeNormalDistribution.value(phi * (x2 - sqrtVar));
      return phi * discountDfDelivery * (underlying / driftDfDelivery * N1 - payoff.strike * N2);
    }

    private def C(eta: Int, phi: Int): Double = {
      val HS = payoff.barrier / underlying;
      val powHS0 = math.pow(HS, 2 * mu);
      val powHS1 = powHS0 * HS * HS;
      val y1 = math.log(payoff.barrier * HS / payoff.strike) / sqrtVar + muSigma;
      val N1 = CumulativeNormalDistribution.value(eta * y1);
      val N2 = CumulativeNormalDistribution.value(eta * (y1 - sqrtVar));
      return phi * discountDfDelivery * (underlying / driftDfDelivery * powHS1 * N1 - payoff.strike * powHS0 * N2);
    }

    private def D(eta: Int, phi: Int): Double = {
      val HS = payoff.barrier / underlying;
      val powHS0 = math.pow(HS, 2 * mu);
      val powHS1 = powHS0 * HS * HS;
      val y2 = math.log(payoff.barrier / underlying) / sqrtVar + muSigma;
      val N1 = CumulativeNormalDistribution.value(eta * y2);
      val N2 = CumulativeNormalDistribution.value(eta * (y2 - sqrtVar));
      return phi * discountDfDelivery * (underlying / driftDfDelivery * powHS1 * N1 - payoff.strike * powHS0 * N2);
    }

    private def E(eta: Int): Double = {
      if (payoff.rebate > 0) {
        val powHS0 = math.pow(payoff.barrier / underlying, 2 * mu);
        val x2 =
          math.log(underlying / payoff.barrier) / sqrtVar + muSigma;
        val y2 =
          math.log(payoff.barrier / underlying) / sqrtVar + muSigma;
        val N1 = CumulativeNormalDistribution.value(eta * (x2 - sqrtVar));
        val N2 = CumulativeNormalDistribution.value(eta * (y2 - sqrtVar));
        return payoff.rebate * discountDfDelivery * (N1 - powHS0 * N2);
      } else {
        return 0.0;
      }
    }

    private def F(eta: Int): Double = {
      if (payoff.rebate > 0) {
        val m = mu
        val lambda = math.sqrt(m * m - 2.0 * math.log(discountDfRebate) / variance) //here df to expiry, not payment!
        val HS = payoff.barrier / underlying;
        val powHSplus = math.pow(HS, m + lambda);
        val powHSminus = math.pow(HS, m - lambda);

        val z = math.log(payoff.barrier / underlying) / sqrtVar + lambda * sqrtVar;

        val N1 = CumulativeNormalDistribution.value(eta * z);
        val N2 = CumulativeNormalDistribution.value(eta * (z - 2.0 * lambda * sqrtVar));
        return payoff.rebate * (powHSplus * N1 + powHSminus * N2);
      } else {
        return 0.0;
      }
    }
  }

}
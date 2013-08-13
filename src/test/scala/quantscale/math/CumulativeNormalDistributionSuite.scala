package test.quantscale.math

import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.analytic._
import quantscale.math.CumulativeBivariateNormalDistribution
import quantscale.math.GenzWestLeFlochBVD
import quantscale.math.GenzWestBVD
import quantscale.analytic.Price
import quantscale.analytic.DeltaGreek
import quantscale.analytic.GammaGreek

@RunWith(classOf[JUnitRunner])
class CumulativeNormalDistributionSuite extends FunSuite {

  test("bivariate-derivative") {
    val t1 = 0.21 //(30.0*3)/365//28/11 - 6/9
    val T2 = 0.30 //(30.0*3+15)/365 //21/12-6/9
    val q = 0.233172168 / 100
    val r = 0.695744403 / 100
    val vol2b = 9.44 / 100
    val vol1 = vol2b * 1.194
    val B = 1.2360
    val vol2 = math.sqrt((vol2b * vol2b * T2 - vol1 * vol1 * t1) / (T2 - t1)) //forward vol
    val varT2 = vol2 * vol2 * (T2 - t1) + vol1 * vol1 * t1
    val sqrtVarT2 = math.sqrt(varT2)
    val vart1 = vol1 * vol1 * t1
    val sqrtVart1 = math.sqrt(vart1)
    val rho = math.sqrt(vart1 / varT2)
    for (i <- 0 until 100) {
      val S0 = 1.26 + i * 0.01
      val eps = 1e-5
      var S = S0 + eps
      var e1 = (math.log(S / B) + (r - q) * t1 + 0.5 * vart1) / sqrtVart1
      var g1 = (math.log(S / B) + (r - q) * T2 + 0.5 * varT2) / sqrtVarT2
      val bvd_up = GenzWestLeFlochBVD.value(e1, g1, rho)
      S = S0 - eps
      e1 = (math.log(S / B) + (r - q) * t1 + 0.5 * vart1) / sqrtVart1
      g1 = (math.log(S / B) + (r - q) * T2 + 0.5 * varT2) / sqrtVarT2
      val bvd_down = GenzWestLeFlochBVD.value(e1, g1, rho)
      val der_num = (bvd_up - bvd_down) / (2 * eps)
      println("NUM " + der_num)
      S = S0
      e1 = (math.log(S / B) + (r - q) * t1 + 0.5 * vart1) / sqrtVart1
      g1 = (math.log(S / B) + (r - q) * T2 + 0.5 * varT2) / sqrtVarT2
      val dmda = math.exp(-0.5 * e1 * e1) / math.sqrt(2 * math.Pi) * CodyCND.value((g1 - rho * e1) / math.sqrt(1 - rho * rho))
      val dmdb = math.exp(-0.5 * g1 * g1) / math.sqrt(2 * math.Pi) * CodyCND.value((e1 - rho * g1) / math.sqrt(1 - rho * rho))
      val der_an = 1.0 / (sqrtVart1 * S) * dmda + 1.0 / (sqrtVarT2 * S) * dmdb
      println("ANA " + der_an + " ratio=" + der_num / der_an)
      assert(math.abs(der_num - der_an) < eps)
    }
  }

  test("bivariate-symmetry") {
    val x = 0.1301020461100371
    val y = 0.08933866025092878
    val rho = -0.9956254315755497
    val z1 = GenzWestBVD.value(x, y, rho)
    val z2 = GenzWestBVD.value(y, x, rho)
    assert(math.abs(z1 - z2) < 1e-15)
    println(z1)
    println(z2)
  }
  test("partialbarrier-gamma-bad") {
    val isCall = false
    val isUp = false
    val K = 1.3500
    val B = 1.2360
    val S0 = 1.26 //1.2619
    val vol2 = 9.44 / 100
    val vol1 = vol2 * 1.194
    println(vol1)
    val q = 0.233172168 / 100
    val r = 0.695744403 / 100
    val t1 = 0.21 //(30.0*3)/365//28/11 - 6/9
    val T2 = 0.30 //(30.0*3+15)/365 //21/12-6/9
    println("S\tPrice\tDelta\tGamma\tDeltaAnalytic\tGammaAnalytic\tDeltaError\tGammaError")
    val price = new Price()
    val deltaAn = new DeltaGreek()
    val gammaAn = new GammaGreek()
    val measures = Array[BSMMeasure](price, deltaAn, gammaAn)
    val priceUp = new Price()
    val measuresUp = Array[BSMMeasure](priceUp)
    val priceDown = new Price()
    val measuresDown = Array[BSMMeasure](priceDown)
    val eps = 1e-5
    for (i <- 0 until 100) {
      var S = S0 + 1e-3 * i
      val details = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S, vol1, vol2, r, q, t1, T2, measures, GenzWestBVD, CodyCND)

      val c = price.value
      val detailsUp = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S * (1 + eps), vol1, vol2, r, q, t1, T2, measuresUp, GenzWestBVD, CodyCND)
      val detailsDown = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S * (1 - eps), vol1, vol2, r, q, t1, T2, measuresDown, GenzWestBVD, CodyCND)
      val delta = (priceUp.value - priceDown.value) / (2 * S * eps)
      val gamma = (priceUp.value + priceDown.value - 2 * c) / (S * S * eps * eps)
      var j = 0
      //      while (j < details.length) {
      ////        val deriv = (detailsUp(j) -detailsDown(j))/(2*S*1e-5)
      ////        print(deriv+" ")
      //        val d2 = (detailsUp(j) + detailsDown(j)-2*details(j))/(S*S*eps*eps)
      //        j+=1
      //      }
      //      println
      println(S + "\t" + c + "\t" + delta + "\t" + gamma + "\t" + deltaAn.value + "\t" + gammaAn.value + "\t" + (deltaAn.value - delta) + "\t" + (gammaAn.value - gamma))
    }
  }

  test("partialbarrier-gamma-ok") {
    val isCall = false
    val isUp = false
    val K = 1.3500
    val B = 1.2360
    val S0 = 1.26 //1.2619
    val vol2 = 9.44 / 100
    val vol1 = vol2 * 1.194
    println(vol1)
    val q = 0.233172168 / 100
    val r = 0.695744403 / 100
    val t1 = 0.21 //(30.0*3)/365//28/11 - 6/9
    val T2 = 0.30 //(30.0*3+15)/365 //21/12-6/9
    println("S\tPrice\tDelta\tGamma\tDeltaAnalytic\tGammaAnalytic\tDeltaError\tGammaError")
    val price = new Price()
    val deltaAn = new DeltaGreek()
    val gammaAn = new GammaGreek()
    val measures = Array[BSMMeasure](price, deltaAn, gammaAn)
    val priceUp = new Price()
    val measuresUp = Array[BSMMeasure](priceUp)
    val priceDown = new Price()
    val measuresDown = Array[BSMMeasure](priceDown)
    val eps = 1e-5
    val cnd = CodyCND //HartCND
    val bvd = GenzWestLeFlochBVD
    for (i <- 0 until 100) {
      var S = S0 + 1e-3 * i
      val details = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S, vol1, vol2, r, q, t1, T2, measures, bvd, cnd)

      val c = price.value
      val detailsUp = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S * (1 + eps), vol1, vol2, r, q, t1, T2, measuresUp, bvd, cnd)
      val detailsDown = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S * (1 - eps), vol1, vol2, r, q, t1, T2, measuresDown, bvd, cnd)
      val delta = (priceUp.value - priceDown.value) / (2 * S * eps)
      val gamma = (priceUp.value + priceDown.value - 2 * c) / (S * S * eps * eps)
      var j = 0
      //      while (j < details.length) {
      ////        val deriv = (detailsUp(j) -detailsDown(j))/(2*S*1e-5)
      ////        print(deriv+" ")
      //        val d2 = (detailsUp(j) + detailsDown(j)-2*details(j))/(S*S*eps*eps)
      //        print(d2+" ")
      //        j+=1
      //      }
      //      println
      println(S + "\t" + c + "\t" + delta + "\t" + gamma + "\t" + deltaAn.value + "\t" + gammaAn.value + "\t" + (deltaAn.value - delta) + "\t" + (gammaAn.value - gamma))
    }
  }

  test("partialbarrier-upoutcall-gamma-ok") {
    val isCall = true
    val isUp = true
    val K = 1.2600
    val B = 1.3750
    val S0 = 1.35 //1.2619
    val vol2 = 9.44 / 100
    val vol1 = vol2 * 1.194
    println(vol1)
    val r = 0.695744403 / 100
    val q = 0.233172168 / 100
    val t1 = 0.21 //(30.0*3)/365//28/11 - 6/9
    val T2 = 0.30 //(30.0*3+15)/365 //21/12-6/9
    println("S\tPrice\tDelta\tGamma\tDeltaAnalytic\tGammaAnalytic\tDeltaError\tGammaError")
    val price = new Price()
    val deltaAn = new DeltaGreek()
    val gammaAn = new GammaGreek()
    val measures = Array[BSMMeasure](price, deltaAn, gammaAn)
    val priceUp = new Price()
    val measuresUp = Array[BSMMeasure](priceUp)
    val priceDown = new Price()
    val measuresDown = Array[BSMMeasure](priceDown)
    val eps = 1e-5
    val cnd = CodyCND //HartCND
    val bvd = GenzWestLeFlochBVD
    for (i <- 0 until 100) {
      var S = S0 + 1e-3 * i
      val details = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S, vol1, vol2, r, q, t1, T2, measures, bvd, cnd)

      val c = price.value
      val detailsUp = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S * (1 + eps), vol1, vol2, r, q, t1, T2, measuresUp, bvd, cnd)
      val detailsDown = BlackScholesMertonBarrier.priceOutPartialEnd(isCall, isUp, K, B, S * (1 - eps), vol1, vol2, r, q, t1, T2, measuresDown, bvd, cnd)
      val delta = (priceUp.value - priceDown.value) / (2 * S * eps)
      val gamma = (priceUp.value + priceDown.value - 2 * c) / (S * S * eps * eps)
      var j = 0
      //      while (j < details.length) {
      ////        val deriv = (detailsUp(j) -detailsDown(j))/(2*S*1e-5)
      ////        print(deriv+" ")
      //        val d2 = (detailsUp(j) + detailsDown(j)-2*details(j))/(S*S*eps*eps)
      //        print(d2+" ")
      //        j+=1
      //      }
      //      println
      println(S + "\t" + c + "\t" + delta + "\t" + gamma + "\t" + deltaAn.value + "\t" + gammaAn.value + "\t" + (deltaAn.value - delta) + "\t" + (gammaAn.value - gamma))
    }
  }

  test("bivariate-monotonicity") {
    var x = -6.9410
    var y = 6.9410
    var corr = -0.999

    val iter = 100
    println("x\tAlgorithm\tValue")
    for (i <- 0 until iter) {
      var v0 = GenzWestBVD.value(x, y, corr, SchonfelderCND)
      var v1 = GenzWestLeFlochBVD.value(x, y, corr, SchonfelderCND)
      println(x + "\tGenzWest\t" + v0)
      println(x + "\tGenzWestFixed\t" + v1)
      x = x + 1e-6
    }

  }

  test("performance") {
    val low = -37 // -12.0
    val high = 8.0 // 5.0
    val maxIter = 10000
    for (k<- 0 until 10) {

      var start = System.nanoTime()
      var sum = 0.0
      for (i <- 0 until maxIter) {
        val x = low + (high - low) * i / (maxIter - 1)
         sum += CodyCND.value(x)
      }
      var elapsed = System.nanoTime()-start
      println("Cody  "+elapsed*1e-9+" "+sum)
      start = System.nanoTime()
      sum = 0.0
      for (i <- 0 until maxIter) {
        val x = low + (high - low) * i / (maxIter - 1)
        sum += CodyInlinedCND.value(x)
      }
      elapsed = System.nanoTime()-start
      println("CodyI "+elapsed*1e-9+" "+sum)
      start = System.nanoTime()
      sum = 0.0
      for (i <- 0 until maxIter) {
        val x = low + (high - low) * i / (maxIter - 1)
        sum += SchonfelderCND.value(x)
      }
      elapsed = System.nanoTime()-start
      println("Schon "+elapsed*1e-9+" "+sum)
      start = System.nanoTime()
      sum = 0.0
      for (i <- 0 until maxIter) {
        val x = low + (high - low) * i / (maxIter - 1)
        sum += OouraCND.value(x)
      }
      elapsed = System.nanoTime()-start
      println("Ooura "+elapsed*1e-9+" "+sum)
      start = System.nanoTime()
      sum = 0.0
      for (i <- 0 until maxIter) {
        val x = low + (high - low) * i / (maxIter - 1)
        sum += JohnsonCND.value(x)
      }
      elapsed = System.nanoTime()-start
      println("Johns "+elapsed*1e-9+" "+sum)

    }
  }
  test("absolute-middle-precision") {
    val low = -37 // -12.0
    val high = 8.0 // 5.0
    val maxIter = 1000
    val accuracy = 1e-15

    for (i <- 0 until maxIter) {
      val x = low + (high - low) * i / (maxIter - 1)
      val cndCody = CodyCND.value(x)
      val cndHart = HartCND.value(x)
      val cndSchon = SchonfelderCND.value(x)
      val cndHill = HillAS66CND.value(x)
      val cndOo = OouraCND.value(x)
      val cndJ = JohnsonCND.value(x)
      assert(math.abs(cndCody - cndHart) < accuracy, "cody vs hart " + cndCody + " " + cndHart + " at " + x)
      assert(math.abs(cndSchon - cndHart) < accuracy, "Schonfelder vs hart " + cndSchon + " " + cndHart + " at " + x)
      assert(math.abs(cndHill - cndHart) < 1e-11, "Hill vs hart " + cndHill + " " + cndHart + " at " + x)
      assert(math.abs(cndOo - cndHart) < accuracy, "Ooura vs hart" + cndOo + " " + cndHart + " at" + x)
      assert(math.abs(cndJ - cndHart) < accuracy, "Johnson vs hart " + cndJ + " " + cndHart + " at " + x)
    }
  }

  test("relative-high-precision") {
    val low = -37 // -12.0
    val high = 8.0 // 5.0
    val maxIter = (37+8)*4+1//100
    val accuracy = 1e-15
    val eps = 1e-6
//    println("x\tApproximation\tValue\tRelativeError\tDiff")
println("x\tCody\tHart\tHartQL\tSchon\tHill\tOoura\tJohnson")
    
        for (i <- 0 until maxIter) {
          val x = low + (high - low) * i / (maxIter - 1)
//    val x = -16.6
    val cndCody = CodyCND.value(x)
    val cndHart = HartCND.value(x)
    val cndHartQL = HartCND.improve(x, cndHart)
    val cndSchon = SchonfelderCND.value(x)
    val cndHill = HillAS66CND.value(x)
    val cndOo = OouraCND.value(x)
    val cndJ = JohnsonCND.value(x)
    val cndCodyOpp = 1.0 - CodyCND.value(-x)
    val cndHartOpp = 1.0 - HartCND.value(-x)
    val cndSchonOpp = 1.0 - SchonfelderCND.value(-x)
    val cndHillOpp = 1.0 - HillAS66CND.value(-x)
    val y = x + eps
    val dCody = CodyCND.value(y) - CodyCND.value(x)
    val dHart = HartCND.value(y) - HartCND.value(x)
    val dSchon = SchonfelderCND.value(y) - SchonfelderCND.value(x)
    val dHill = HillAS66CND.value(y) - HillAS66CND.value(x)
    val dOo = OouraCND.value(y) - OouraCND.value(x)
    val dCodyOpp = CodyCND.value(-x) - CodyCND.value(-y)
    val dHartOpp = HartCND.value(-x) - HartCND.value(-y)
    val dSchonOpp = SchonfelderCND.value(-x) - SchonfelderCND.value(-y)
    val dHillOpp = HillAS66CND.value(-x) - HillAS66CND.value(-y)
    println(x + "\t"+cndCody+"\t"+cndHart+"\t"+cndHartQL+"\t"+cndSchon+"\t"+cndHill+"\t"+cndOo+"\t"+cndJ)
    
    //    println(x + "\tCody\t" + cndCody + "\t" + 0.0 + "\t" + dCody)
//    println(x + "\tHart\t" + cndHart + "\t" + math.abs((cndCody - cndHart) / cndCody) + "\t" + dHart)
//    println(x + "\tHartQL\t" + cndHartQL + "\t" + math.abs((cndCody - cndHartQL) / cndCody) + "\t" + dHart)
//    println(x + "\tSchonfelder\t" + cndSchon + "\t" + math.abs((cndSchon - cndCody) / cndCody) + "\t" + dSchon)
//    println(x + "\tHill\t" + cndHill + "\t" + math.abs((cndHill - cndCody) / cndCody) + "\t" + dHill)
//    println(x + "\tOoura\t" + cndOo + "\t" + math.abs((cndOo - cndCody) / cndCody) + "\t" + dOo)
//    println(x + "\t1-Cody\t" + cndCodyOpp + "\t" + math.abs((cndCody - cndCodyOpp) / cndCody) + "\t" + dCodyOpp)
//    println(x + "\t1-Hart\t" + cndHartOpp + "\t" + math.abs((cndCody - cndHartOpp) / cndCody) + "\t" + dHartOpp)
//    println(x + "\t1-Schonfelder\t" + cndSchonOpp + "\t" + math.abs((cndSchonOpp - cndCody) / cndCody) + "\t" + dSchonOpp)
//    println(x + "\t1-Hill\t" + cndHillOpp + "\t" + math.abs((cndHillOpp - cndCody) / cndCody) + "\t" + dHillOpp)
//          assert(math.abs((cndCody-cndHart)/cndHart)<accuracy,"cody vs hart "+cndCody+" "+cndHart+" at "+x)
//          assert(math.abs((cndSchon-cndHart)/cndHart)<accuracy,"Schonfelder vs hart "+cndSchon+" "+cndHart+" at "+x)
//          assert(math.abs((cndHill-cndHart)/cndHart)<1e-11,"Hill vs hart "+cndHill+" "+cndHart+" at "+x)
//                    assert(math.abs((cndOo-cndHart)/cndHart)<accuracy,"Ooura vs hart "+cndOo+" "+cndHart+" at "+x)
        }
  }
  
   test("relative-high-precision-exact") {
    val low = -37 // -12.0
    val high = 8.0 // 5.0
    val maxIter = (37+8)*4+1//100
    val accuracy = 1e-15
    println("x\tMethod\tValue")

        for (i <- 0 until maxIter) {
          val x = low + (high - low) * i / (maxIter - 1)
//    val x = -16.6
    val cndCody = CodyCND.value(x)
          println(x + "\tCody\t"+cndCody)
    val cndHart = HartCND.value(x)
          println(x + "\tHart\t"+cndHart)
    val cndHartQL = HartCND.improve(x, cndHart)
          println(x + "\tHartQL\t"+cndHartQL)
    val cndSchon = SchonfelderCND.value(x)
          println(x + "\tSchonfelder\t"+cndSchon)
    val cndHill = HillAS66CND.value(x)
          println(x + "\tHill\t"+cndHill)
    val cndOo = OouraCND.value(x)
          println(x + "\tOoura\t"+cndOo)
    val cndJ = JohnsonCND.value(x)
          println(x + "\tJohnson\t"+cndJ)

        }
        }

     test("relative-high-precision-exact-200") {
       val low = -37 // -12.0
       val high = 8.0 // 5.0
       val maxIter = 200//100
       val accuracy = 1e-15
       println("x\tMethod\tValue")

       for (i <- 0 until maxIter) {
         val x = low + (high - low) * i / (maxIter - 1)
         //    val x = -16.6
//         val cndCody = CodyCND.value(x)
//         println(x + "\tCody\t"+cndCody)
//         val cndHart = HartCND.value(x)
//         println(x + "\tHart\t"+cndHart)
//         val cndHartQL = HartCND.improve(x, cndHart)
//         println(x + "\tHartQL\t"+cndHartQL)
//         val cndSchon = SchonfelderCND.value(x)
//         println(x + "\tSchonfelder\t"+cndSchon)
//         val cndHill = HillAS66CND.value(x)
//         println(x + "\tHill\t"+cndHill)
//         val cndOo = OouraCND.value(x)
//         println(x + "\tOoura\t"+cndOo)
         val cndJ = JohnsonCND.value(x)
         println(x + "\tJohnson\t"+cndJ)


       }
  }

      test("relative-high-precision-rand") {
    val low = -37 // -12.0
    val high = 8.0 // 5.0
    val maxIter = 200
    val accuracy = 1e-15
 println("x\tMethod\tValue")
//println("x\tCody\tHart\tHartQL\tSchon\tHill\tOoura")
    
        for (i <- 0 until maxIter) {
          val x = low + (high - low) * i / (maxIter - 1)
//    val x = -16.6
    val cndCody = CodyCND.value(x)
    val cndHart = HartCND.value(x)
    val cndHartQL = HartCND.improve(x, cndHart)
    val cndSchon = SchonfelderCND.value(x)
    val cndHill = HillAS66CND.value(x)
    val cndOo = OouraCND.value(x)
    val cndSun = SunCND.value(x)
    val cndCodyOpp = 1.0 - CodyCND.value(-x)
    val cndHartOpp = 1.0 - HartCND.value(-x)
    val cndSchonOpp = 1.0 - SchonfelderCND.value(-x)
    val cndHillOpp = 1.0 - HillAS66CND.value(-x)
//    println(x + "\t"+cndCody+"\t"+cndHart+"\t"+cndHartQL+"\t"+cndSchon+"\t"+cndHill+"\t"+cndOo)
    println(x+"\tCody\t"+cndCody)
    println(x+"\tSchonfelder\t"+cndSchon)
    println(x+"\tOoura\t"+cndOo)
    println(x+"\tSun\t"+cndSun)
        }
  }

  test("univariate-gaussian") {
    for (i <- 0 until 80) {
      val q =  i.toDouble /1.6
      val z = math.exp(-q*q*0.5)
      println(z+",")
    }
  }
}
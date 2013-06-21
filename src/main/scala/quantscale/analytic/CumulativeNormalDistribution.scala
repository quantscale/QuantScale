package quantscale.analytic
import quantscale.fdm.Epsilon
import quantscale.math.Function1D

/**
 * Reference: Cody, W.D. (1993). ALGORITHM 715: SPECFUN - A Portable
 * FORTRAN Package of Special Function Routines and Test Drivers" ACM
 * Transactions on Mathematical Software. 19, 22-32.
 *
 * This function evaluates the normal distribution function: The main
 * computation evaluates near-minimax approximations derived from those
 * in "Rational Chebyshev approximations for the error function" by W.
 * J. Cody, Math. Comp., 1969, 631-637. This transportable program uses
 * rational functions that theoretically approximate the normal
 * distribution function to at least 18 significant decimal digits. The
 * accuracy achieved depends on the arithmetic system, the compiler, the
 * intrinsic functions, and proper selection of the machine-dependent
 * constants.
 *
 *
 * Mathematical Constants:
 *
 * sqrpi = 1 / sqrt(2*pi), root32 = sqrt(32), thrsh = the argument for
 * which pnorm(thrsh,0,1) = 0.75.
 */

object CodyCND extends Function1D {

  private final val DBL_EPSILON = Epsilon.MACHINE_EPSILON; //2.2204460492503131e-16;
  private final val DBL_MIN = Double.MinValue; //2.2250738585072014e-308;

  private final val c = Array(0.39894151208813466764,
    8.8831497943883759412, 93.506656132177855979,
    597.27027639480026226, 2494.5375852903726711,
    6848.1904505362823326, 11602.651437647350124,
    9842.7148383839780218, 1.0765576773720192317e-8);

  private final val d = Array(22.266688044328115691,
    235.38790178262499861, 1519.377599407554805,
    6485.558298266760755, 18615.571640885098091,
    34900.952721145977266, 38912.003286093271411,
    19685.429676859990727);

  private final val p = Array(0.21589853405795699,
    0.1274011611602473639, 0.022235277870649807,
    0.001421619193227893466, 2.9112874951168792e-5,
    0.02307344176494017303);

  private final val q = Array(1.28426009614491121,
    0.468238212480865118, 0.0659881378689285515,
    0.00378239633202758244, 7.29751555083966205e-5);

  private final val a = Array(2.2352520354606839287,
    161.02823106855587881, 1067.6894854603709582,
    18154.981253343561249, 0.065682337918207449113);

  private final val b = Array(47.20258190468824187,
    976.09855173777669322, 10260.932208618978205,
    45507.789335026729956);

  final val one = 1.0;
  final val half = 0.5;
  final val zero = 0.0;
  final val sixten = 1.6;
  final val sqrpi = 0.39894228040143267794;
  final val thrsh = 0.66291;
  final val root32 = 5.656854248;

  def pnorms(x: Double): Double = {

    var xden: Double = 0;
    var temp: Double = 0;
    var xnum: Double = 0;
    var result: Double = 0;
    var ccum: Double = 0;
    var del: Double = 0;
    val min = DBL_MIN;
    val eps = DBL_EPSILON * .5;
    var xsq: Double = 0;
    var y: Double = 0;
    var i: Int = 0;

    y = Math.abs(x);
    if (y <= thrsh) {
      /* Evaluate pnorm for |z| <= 0.66291 */
      xsq = zero;
      if (y > eps) {
        xsq = x * x;
      }
      xnum = a(4) * xsq;
      xden = xsq;
      i = 1;
      while (i <= 3) {
        xnum = (xnum + a(i - 1)) * xsq;
        xden = (xden + b(i - 1)) * xsq;
        i += 1;
      }
      result = x * (xnum + a(3)) / (xden + b(3));
      temp = result;
      result = half + temp;
      ccum = half - temp;
    } else if (y <= root32) {
      /* Evaluate pnorm for 0.66291 <= |z| <= sqrt(32) */
      xnum = c(8) * y;
      xden = y;
      i = 1;
      while (i <= 7) {
        xnum = (xnum + c(i - 1)) * y;
        xden = (xden + d(i - 1)) * y;
        i += 1;
      }
      result = (xnum + c(7)) / (xden + d(7));
      xsq = fint(y * sixten) / sixten;
      del = (y - xsq) * (y + xsq);
      result = Math.exp(-xsq * xsq * half) * Math.exp(-del * half) * result;
      ccum = one - result;
      if (x > zero) {
        temp = result;
        result = ccum;
        ccum = temp;
      }
    } else {
      /* Evaluate pnorm for |z| > sqrt(32) */
      result = zero;
      xsq = one / (x * x);
      xnum = p(5) * xsq;
      xden = xsq;
      i = 1;
      while (i <= 4) {
        xnum = (xnum + p(i - 1)) * xsq;
        xden = (xden + q(i - 1)) * xsq;
        i += 1;
      }
      result = xsq * (xnum + p(4)) / (xden + q(4));
      result = (sqrpi - result) / y;
      xsq = fint(x * sixten) / sixten;
      del = (x - xsq) * (x + xsq);
      result = Math.exp(-xsq * xsq * half) * Math.exp(-del * half) * result;
      ccum = one - result;
      if (x > zero) {
        temp = result;
        result = ccum;
        ccum = temp;
      }
    }
    if (result < min) {
      result = 0.0;
    }
    if (ccum < min) {
      ccum = 0.0;
    }
    return result;
  }

  private def fint(x: Double): Double =

    if (x >= 0.0) StrictMath.floor(x) else -StrictMath.floor(-x);

  def value(x: Double): Double = pnorms(x);
}

object CodyInlinedCND extends Function1D {

  private final val DBL_EPSILON = Epsilon.MACHINE_EPSILON; //2.2204460492503131e-16;
  private final val DBL_MIN = Double.MinValue; //2.2250738585072014e-308;

  private final val c0 = 0.39894151208813466764
  val c1 = 8.8831497943883759412
  val c2 = 93.506656132177855979
  val c3 = 597.27027639480026226
  val c4 = 2494.5375852903726711
  val c5 = 6848.1904505362823326
  val c6 = 11602.651437647350124
  val c7 = 9842.7148383839780218
  val c8 = 1.0765576773720192317e-8

  private final val d0 = 22.266688044328115691
  val d1 = 235.38790178262499861
  val d2 = 1519.377599407554805
  val d3 = 6485.558298266760755
  val d4 = 18615.571640885098091
  val d5 = 34900.952721145977266
  val d6 = 38912.003286093271411
  val d7 = 19685.429676859990727

  private final val p0 = 0.21589853405795699
  val p1 = 0.1274011611602473639
  val p2 = 0.022235277870649807
  val p3 = 0.001421619193227893466
  val p4 = 2.9112874951168792e-5
  val p5 = 0.02307344176494017303

  private final val q0 = 1.28426009614491121
  val q1 = 0.468238212480865118
  val q2 = 0.0659881378689285515
  val q3 = 0.00378239633202758244
  val q4 = 7.29751555083966205e-5

  private final val a0 = 2.2352520354606839287
  val a1 = 161.02823106855587881
  val a2 = 1067.6894854603709582
  val a3 = 18154.981253343561249
  val a4 = 0.065682337918207449113

  private final val b0 = 47.20258190468824187
  val b1 = 976.09855173777669322
  val b2 = 10260.932208618978205
  val b3 = 45507.789335026729956

  final val one = 1.0;
  final val half = 0.5;
  final val zero = 0.0;
  final val sixten = 1.6;
  final val sqrpi = 0.39894228040143267794;
  final val thrsh = 0.66291;
  final val root32 = 5.656854248;

  def pnorms(x: Double): Double = {

    var xden: Double = 0;
    var temp: Double = 0;
    var xnum: Double = 0;
    var result: Double = 0;
    var ccum: Double = 0;
    val min = DBL_MIN;
    val eps = DBL_EPSILON * .5;
    var xsq: Double = 0;
    var y: Double = 0;

    y = Math.abs(x);
    if (y <= thrsh) {
      /* Evaluate pnorm for |z| <= 0.66291 */
      xsq = zero;
      if (y > eps) {
        xsq = x * x;
      }
      xnum = a4 * xsq;
      xden = xsq;
      xnum = (xnum + a0) * xsq
      xden = (xden + b0) * xsq;
      xnum = (xnum + a1) * xsq
      xden = (xden + b1) * xsq;
      xnum = (xnum + a2) * xsq
      xden = (xden + b2) * xsq;

      result = x * (xnum + a3) / (xden + b3);
      temp = result;
      result = half + temp;
      ccum = half - temp;
    } else if (y <= root32) {
      /* Evaluate pnorm for 0.66291 <= |z| <= sqrt(32) */
      xnum = c8 * y;
      xden = y;
      xnum = (xnum + c0) * y;
      xden = (xden + d0) * y;
      xnum = (xnum + c1) * y;
      xden = (xden + d1) * y;
      xnum = (xnum + c2) * y;
      xden = (xden + d2) * y;
      xnum = (xnum + c3) * y;
      xden = (xden + d3) * y;
      xnum = (xnum + c4) * y;
      xden = (xden + d4) * y;
      xnum = (xnum + c5) * y;
      xden = (xden + d5) * y;
      xnum = (xnum + c6) * y;
      xden = (xden + d6) * y;

      result = (xnum + c7) / (xden + d7);

      result = NormalDistribution.denormalizedValue(y) * result;
      ccum = one - result;
      if (x > zero) {
        temp = result;
        result = ccum;
        ccum = temp;
      }
    } else {
      /* Evaluate pnorm for |z| > sqrt(32) */
      result = zero;
      xsq = one / (x * x);
      xnum = p5 * xsq;
      xden = xsq;
      xnum = (xnum + p0) * xsq;
      xden = (xden + q0) * xsq;
      xnum = (xnum + p1) * xsq;
      xden = (xden + q1) * xsq;
      xnum = (xnum + p2) * xsq;
      xden = (xden + q2) * xsq;
      xnum = (xnum + p3) * xsq;
      xden = (xden + q3) * xsq;

      result = xsq * (xnum + p4) / (xden + q4);
      result = (sqrpi - result) / y;
      result *= NormalDistribution.denormalizedValue(y)
      ccum = one - result;
      if (x > zero) {
        temp = result;
        result = ccum;
        ccum = temp;
      }
    }
    if (result < min) {
      result = 0.0;
    }
    if (ccum < min) {
      ccum = 0.0;
    }
    return result;
  }

  private def fint(x: Double): Double =

    if (x >= 0.0) StrictMath.floor(x) else -StrictMath.floor(-x);

  def value(x: Double): Double = pnorms(x);
}

/**
 * The original Hill AS66 algorithm (without Alan Miller mistakes)
 */
object HillAS66CND extends Function1D {
  val CON = 1.28d
  val HALF = 0.5d
  //original val LTONE = 5.0d
  val ONE = 1.0d

  //original,val UTZERO = 12.5d 
  val ZERO = 0.0d

  //proposed by Hill in AS66
  val LTONE = 8.3
  val UTZERO = 37.52

  def value(X: Double): Double = {
    var Z = X
    var alnorm = 0.0
    var up = false
    if (Z < 0) {
      up = !up
      Z = -Z
    }
    if (Z <= LTONE || up && Z <= UTZERO) {
      var Y = HALF * Z * Z
      if (Z >= CON) {
        alnorm = 0.398942280385e0 * math.exp(-Y) /
          (Z - 3.8052e-8 + 1.00000615302e0 /
            (Z + 3.98064794e-4 + 1.98615381364e0 /
              (Z - 0.151679116635e0 + 5.29330324926e0 /
                (Z + 4.8385912808e0 - 15.1508972451e0 /
                  (Z + 0.742380924027e0 + 30.789933034e0 /
                    (Z + 3.99019417011e0))))))
      } else {
        alnorm = HALF - Z * (0.398942280444e0 - 0.399903438504e0 * Y /
          (Y + 5.75885480458e0 - 29.8213557808e0 /
            (Y + 2.62433121679e0 + 48.6959930692e0 /
              (Y + 5.92885724438e0))))
      }
    } else {
      alnorm = 0.0
    }

    if (!up) {
      alnorm = 1.0 - alnorm
    }
    return alnorm
  }
}

object SchonfelderCND extends Function1D {
  val RTWO = 1.414213562373095048801688724209e0
  val IM = 24
  val A = Array[Double](
    6.10143081923200417926465815756E-1, -4.34841272712577471828182820888E-1, 1.76351193643605501125840298123E-1, -6.0710795609249414860051215825E-2, 1.7712068995694114486147141191E-2,
    -4.321119385567293818599864968E-3, 8.54216676887098678819832055E-4, -1.27155090609162742628893940E-4, 1.1248167243671189468847072E-5, 3.13063885421820972630152E-7,
    -2.70988068537762022009086E-7, 3.0737622701407688440959E-8, 2.515620384817622937314E-9, -1.028929921320319127590E-9, 2.9944052119949939363E-11, 2.6051789687266936290E-11,
    -2.634839924171969386E-12, -6.43404509890636443E-13, 1.12457401801663447E-13, 1.7281533389986098E-14, -4.264101694942375E-15, -5.45371977880191E-16, 1.58697607761671E-16,
    2.0899837844334E-17, -5.900526869409E-18, -9.41893387554E-19, 2.14977356470E-19, 4.6660985008E-20, -7.243011862E-21, -2.387966824E-21, 1.91177535E-22, 1.20482568E-22, -6.72377E-25,
    -5.747997E-24, -4.28493E-25, 2.44856E-25, 4.3793E-26, -8.151E-27, -3.089E-27, 9.3E-29, 1.74E-28, 1.6E-29, -8.0E-30, -2.0E-30)

  def value(x: Double): Double = phi(x)
  /*
         * Normal distribution probabilities accurate to 1d-15. Reference: J.L.
         * Schonfelder, Math Comp 32(1978), pp 1232-1240.
         */
  def phi(Z: Double): Double = {
    var I = 0;

    var BM, B, BP = 0.0;
    var P, T, XA = 0.0;

    XA = math.abs(Z) / RTWO;
    if (XA > 100) {
      P = 0;
    } else {
      T = (8 * XA - 30) / (4 * XA + 15);
      BM = 0;
      B = 0;
      //unroll loop
      //      A24 //BP=0, B=0
      //      T*A24 + A23 //BP=0 B=A_24
      //      T*(T*A24 + A23)-A24+A22//BP=A_24, B=T*A24+A23
      //      T*(T*(T*A24 + A23)-A24+A22)-T*A24 + A23+A21//BP=T*A24+A23, B=T*(T*A24 + A23)-A24+A22

      I = IM
      while (I >= 0) {
        BP = B;
        B = BM;
        BM = T * B - BP + A(I);
        I -= 1
      }
      P = NormalDistribution.denormalizedValue(Z) * (BM - BP) / 4;
      //      P = math.exp(-XA * XA) * (BM - BP) / 4;
    }
    if (Z > 0)
      P = 1 - P;
    return P;
  }
}

object OouraCND extends Function1D {
  val RTWO = 1.414213562373095048801688724209e0

  /**
   * Takuya OOURA, Research Institute for Mathematical Sciences Kyoto
   * University derfc function
   * http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
   *
   * @param x
   * @return
   */
  def derfcOoura(x: Double): Double = {
    val y = derfOuraUnscaled(x)
    return if (x < 0) 2 - y else y
  }

  private def derfOuraUnscaled(x: Double): Double = {

    val t = 3.97886080735226 / (Math.abs(x * 0.7071067811865475244008) + 3.97886080735226);
    val u = t - 0.5;
    var y = (((((((((0.00127109764952614092 * u + 1.19314022838340944e-4) * u - 0.003963850973605135) * u - 8.70779635317295828e-4) * u + 0.00773672528313526668) * u + 0.00383335126264887303) * u - 0.0127223813782122755) * u - 0.0133823644533460069) * u + 0.0161315329733252248) * u + 0.0390976845588484035) * u + 0.00249367200053503304;
    y = ((((((((((((y * u - 0.0838864557023001992) * u - 0.119463959964325415) * u + 0.0166207924969367356) * u + 0.357524274449531043) * u + 0.805276408752910567) * u + 1.18902982909273333) * u + 1.37040217682338167) * u + 1.31314653831023098) * u + 1.07925515155856677) * u + 0.774368199119538609) * u + 0.490165080585318424) * u + 0.275374741597376782) * t *
      NormalDistribution.denormalizedValue(x)
    //    NormalDistribution.denormalizedValue(x*RTWO)
    //    Math.exp(-x * x);
    return y;
  }

  def value(x: Double): Double = {
    val y = derfOuraUnscaled(x) // * 0.7071067811865475244008);
    return if (x < 0) 0.5 * y else 1 - 0.5 * y;
  }

}

object HartCND extends Function1D {
  val P0 = 220.2068679123761e0
  val P1 = 221.2135961699311e0
  val P2 = 112.0792914978709e0
  val P3 = 33.91286607838300e0
  val P4 = 6.373962203531650e0
  val P5 = .7003830644436881e0
  val P6 = .3526249659989109e-01
  val Q0 = 440.4137358247522e0
  val Q1 = 793.8265125199484e0
  val Q2 = 637.3336333788311e0
  val Q3 = 296.5642487796737e0
  val Q4 = 86.78073220294608e0
  val Q5 = 16.06417757920695e0
  val Q6 = 1.755667163182642e0
  val Q7 = .8838834764831844e-1
  val CUTOFF = 7.071e0
  val ROOT2PI = 2.506628274631001e0

  def value(Z: Double): Double = {
    val ZABS = math.abs(Z)
    //       |Z| > 37.
    if (ZABS > 37.0) {
      //          PDF = 0.D0
      if (Z > 0.0) {
        //            P = 1.D0
        return 1.0
        //            Q = 0.D0
      } else {
        return 0.0
        //P = 0.D0
        //Q = 1.D0
      }
    }
    //       |Z| <= 37.
    val EXPNTL = NormalDistribution.denormalizedValue(ZABS) //math.exp(-0.5 * ZABS * ZABS)
    val PDF = EXPNTL / ROOT2PI
    //       |Z| < CUTOFF = 10/sqrt(2).
    var P = 0.0
    if (ZABS < CUTOFF) {
      P = EXPNTL * ((((((P6 * ZABS + P5) * ZABS + P4) * ZABS + P3) * ZABS +
        P2) * ZABS + P1) * ZABS + P0) / (((((((Q7 * ZABS + Q6) * ZABS +
          Q5) * ZABS + Q4) * ZABS + Q3) * ZABS + Q2) * ZABS + Q1) * ZABS +
          Q0)
    } //       |Z| >= CUTOFF.
    else {
      P = PDF / (ZABS + 1.0 / (ZABS + 2.0 / (ZABS + 3.0 / (ZABS + 4.0 /
        (ZABS + 0.65)))))
    }
    if (Z < 0.0) {
      //          Q = 1.D0 - P
      return P
    } else {
      //          Q = P
      //          P = 1.D0 - Q
      return 1.0 - P
    }
  }

  //Hart is not all that precise in the left tail, the following improves it
  def improve(x: Double, cn: Double): Double = {
    if (x < 0 && cn <= Epsilon.MACHINE_EPSILON) {

      // use asymptotic abram and stegun (26.2.12)
      // Taken from quantlib, which takes from Peter Jackel's book
      var sum = 1.0;
      val xsqr = x * x;
      var i = 1
      var g = 1.0
      var z, y = 0.0
      var a = Double.MaxValue
      var lasta = 0.0
      do {
        lasta = a;
        z = (4 * i - 3) / xsqr;
        y = z * ((4 * i - 1) / xsqr);
        a = g * (z - y);
        sum -= a;
        g *= y;
        i += 1
        a = StrictMath.abs(a);
      } while (lasta > a && a >= StrictMath.abs(sum * Epsilon.MACHINE_EPSILON));
      //      println(i)
      return -(NormalDistribution.denormalizedValue(x) / ROOT2PI) / x * sum;
    }
    return cn
  }
}

object SunCND extends Function1D {

  // Coefficients for approximation to  erfc in [1.25,1/.35]
  val eRa = Array(
    -9.86494403484714822705e-03,
    -6.93858572707181764372e-01,
    -1.05586262253232909814e01,
    -6.23753324503260060396e01,
    -1.62396669462573470355e02,
    -1.84605092906711035994e02,
    -8.12874355063065934246e01,
    -9.81432934416914548592e00)
  val eSa = Array(
    1.96512716674392571292e01,
    1.37657754143519042600e02,
    4.34565877475229228821e02,
    6.45387271733267880336e02,
    4.29008140027567833386e02,
    1.08635005541779435134e02,
    6.57024977031928170135e00,
    -6.04244152148580987438e-02)
  // Coefficients for approximation to  erfc in [1/.35,28]
  val eRb = Array(
    -9.86494292470009928597e-03,
    -7.99283237680523006574e-01,
    -1.77579549177547519889e01,
    -1.60636384855821916062e02,
    -6.37566443368389627722e02,
    -1.02509513161107724954e03,
    -4.83519191608651397019e02)
  val eSb = Array(
    3.03380607434824582924e01,
    3.25792512996573918826e02,
    1.53672958608443695994e03,
    3.19985821950859553908e03,
    2.55305040643316442583e03,
    4.74528541206955367215e02,
    -2.24409524465858183362e01)
  val RTWO = 1.414213562373095048801688724209e0

  def value(z: Double): Double = {
//    return 0.5 * erfc(-x / RTWO)
    val x = -z/RTWO
    var s, retval, R, S: Double = 0.0
    val abs_x = if (x >= 0.0) x else -x;
    if (abs_x < 1.25)
      retval = 0.5*(1.0 - erf(abs_x));
    else if (abs_x > 28.0)
      retval = 0.0;
    else { // 1.25 < |x| < 28 
      s = 1.0 / (abs_x * abs_x);
      if (abs_x < 2.8571428) { // ( |x| < 1/0.35 ) 
        R = eRa(0) + s * (eRa(1) + s * (eRa(2) + s * (eRa(3) + s * (eRa(4) + s * (eRa(5) + s * (eRa(6) + s * eRa(7)))))));
        S = 1.0 + s * (eSa(0) + s * (eSa(1) + s * (eSa(2) + s * (eSa(3) + s * (eSa(4) + s * (eSa(5) + s * (eSa(6) + s * eSa(7))))))));
      } else { // ( |x| > 1/0.35 )
        R = eRb(0) + s * (eRb(1) + s * (eRb(2) + s * (eRb(3) + s * (eRb(4) + s * (eRb(5) + s * eRb(6))))));
        S = 1.0 + s * (eSb(0) + s * (eSb(1) + s * (eSb(2) + s * (eSb(3) + s * (eSb(4) + s * (eSb(5) + s * eSb(6)))))));
      }
      retval = NormalDistribution.denormalizedValue(z)*math.exp(- 0.5625 + R / S) *0.5/ abs_x;
    }
    return if (x >= 0.0) retval else 1.0 - retval
  }

  def erfc(x: Double): Double = {
    var s, retval, R, S: Double = 0.0
    val abs_x = if (x >= 0.0) x else -x;
    if (abs_x < 1.25)
      retval = 1.0 - erf(abs_x);
    else if (abs_x > 28.0)
      retval = 0.0;
    else { // 1.25 < |x| < 28 
      s = 1.0 / (abs_x * abs_x);
      if (abs_x < 2.8571428) { // ( |x| < 1/0.35 ) 
        R = eRa(0) + s * (eRa(1) + s * (eRa(2) + s * (eRa(3) + s * (eRa(4) + s * (eRa(5) + s * (eRa(6) + s * eRa(7)))))));
        S = 1.0 + s * (eSa(0) + s * (eSa(1) + s * (eSa(2) + s * (eSa(3) + s * (eSa(4) + s * (eSa(5) + s * (eSa(6) + s * eSa(7))))))));
      } else { // ( |x| > 1/0.35 )
        R = eRb(0) + s * (eRb(1) + s * (eRb(2) + s * (eRb(3) + s * (eRb(4) + s * (eRb(5) + s * eRb(6))))));
        S = 1.0 + s * (eSb(0) + s * (eSb(1) + s * (eSb(2) + s * (eSb(3) + s * (eSb(4) + s * (eSb(5) + s * eSb(6)))))));
      }
//      retval = NormalDistribution.denormalizedValue(x*RTWO)*math.exp(- 0.5625 + R / S) / abs_x;
      retval = math.exp(-x * x - 0.5625 + R / S) / abs_x;
    }
    return if (x >= 0.0) retval else 2.0 - retval
  }

  // ====================================================
  // Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
  //
  // Developed at SunSoft, a Sun Microsystems, Inc. business.
  // Permission to use, copy, modify, and distribute this
  // software is freely granted, provided that this notice 
  // is preserved.
  // ====================================================
  //
  //                           x
  //                    2      |\
  //     erf(x)  =  ---------  | exp(-t*t)dt
  //                 sqrt(pi) \| 
  //                           0
  //
  //     erfc(x) =  1-erf(x)
  //  Note that 
  //              erf(-x) = -erf(x)
  //              erfc(-x) = 2 - erfc(x)
  //
  // Method:
  //      1. For |x| in [0, 0.84375]
  //          erf(x)  = x + x*R(x^2)
  //          erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
  //                  = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
  //         where R = P/Q where P is an odd poly of degree 8 and
  //         Q is an odd poly of degree 10.
  //                                               -57.90
  //                      | R - (erf(x)-x)/x | <= 2
  //      
  //
  //         Remark. The formula is derived by noting
  //          erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
  //         and that
  //          2/sqrt(pi) = 1.128379167095512573896158903121545171688
  //         is close to one. The interval is chosen because the fix
  //         point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
  //         near 0.6174), and by some experiment, 0.84375 is chosen to
  //         guarantee the error is less than one ulp for erf.
  //
  //      2. For |x| in [0.84375,1.25], let s = |x| - 1, and
  //         c = 0.84506291151 rounded to single (24 bits)
  //              erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
  //              erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
  //                        1+(c+P1(s)/Q1(s))    if x < 0
  //              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
  //         Remark: here we use the taylor series expansion at x=1.
  //              erf(1+s) = erf(1) + s*Poly(s)
  //                       = 0.845.. + P1(s)/Q1(s)
  //         That is, we use rational approximation to approximate
  //                      erf(1+s) - (c = (single)0.84506291151)
  //         Note that |P1/Q1|< 0.078 for x in [0.84375,1.25]
  //         where 
  //              P1(s) = degree 6 poly in s
  //              Q1(s) = degree 6 poly in s
  //
  //      3. For x in [1.25,1/0.35(~2.857143)], 
  //              erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
  //              erf(x)  = 1 - erfc(x)
  //         where 
  //              R1(z) = degree 7 poly in z, (z=1/x^2)
  //              S1(z) = degree 8 poly in z
  //
  //      4. For x in [1/0.35,28]
  //              erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
  //                      = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
  //                      = 2.0 - tiny            (if x <= -6)
  //              erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
  //              erf(x)  = sign(x)*(1.0 - tiny)
  //         where
  //              R2(z) = degree 6 poly in z, (z=1/x^2)
  //              S2(z) = degree 7 poly in z
  //
  //      Note1:
  //         To compute exp(-x*x-0.5625+R/S), let s be a single
  //         precision number and s := x; then
  //              -x*x = -s*s + (s-x)*(s+x)
  //              exp(-x*x-0.5626+R/S) = 
  //                      exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
  //      Note2:
  //         Here 4 and 5 make use of the asymptotic series
  //                        exp(-x*x)
  //              erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
  //                        x*sqrt(pi)
  //         We use rational approximation to approximate
  //              g(s)=f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
  //         Here is the error bound for R1/S1 and R2/S2
  //              |R1/S1 - f(x)|  < 2**(-62.57)
  //              |R2/S2 - f(x)|  < 2**(-61.52)
  //
  //      5. For inf > x >= 28
  //              erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
  //              erfc(x) = tiny*tiny (raise underflow) if x > 0
  //                      = 2 - tiny if x<0
  //
  //      7. Special case:
  //              erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
  //              erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2, 
  //              erfc/erf(NaN) is NaN
  //

  // Coefficients for approximation to  erf on [0,0.84375]
  val e_efx = 1.28379167095512586316e-01;
  //double efx8=1.02703333676410069053e00;
  val ePp = Array(
    1.28379167095512558561e-01,
    -3.25042107247001499370e-01,
    -2.84817495755985104766e-02,
    -5.77027029648944159157e-03,
    -2.37630166566501626084e-05)
  val eQq = Array(
    3.97917223959155352819e-01,
    6.50222499887672944485e-02,
    5.08130628187576562776e-03,
    1.32494738004321644526e-04,
    -3.96022827877536812320e-06)
  // Coefficients for approximation to  erf  in [0.84375,1.25] 
  val ePa = Array(
    -2.36211856075265944077e-03,
    4.14856118683748331666e-01,
    -3.72207876035701323847e-01,
    3.18346619901161753674e-01,
    -1.10894694282396677476e-01,
    3.54783043256182359371e-02,
    -2.16637559486879084300e-03)
  val eQa = Array(
    1.06420880400844228286e-01,
    5.40397917702171048937e-01,
    7.18286544141962662868e-02,
    1.26171219808761642112e-01,
    1.36370839120290507362e-02,
    1.19844998467991074170e-02)
  val e_erx = 8.45062911510467529297e-01

  /**
   * Error function.
   * Based on C-code for the error function developed at Sun Microsystems.
   * @author Jaco van Kooten
   */
  def erf(x: Double): Double = {

    var P, Q, s, retval: Double = 0.0
    val abs_x = if (x >= 0.0) x else -x
    if (abs_x < 0.84375) { // 0 < |x| < 0.84375
      if (abs_x < 3.7252902984619141e-9) // |x| < 2**-28
        retval = abs_x + abs_x * e_efx;
      else {
        s = x * x;
        P = ePp(0) + s * (ePp(1) + s * (ePp(2) + s * (ePp(3) + s * ePp(4))));
        Q = 1.0 + s * (eQq(0) + s * (eQq(1) + s * (eQq(2) + s * (eQq(3) + s * eQq(4)))));
        retval = abs_x + abs_x * (P / Q);
      }
    } else if (abs_x < 1.25) { // 0.84375 < |x| < 1.25
      s = abs_x - 1.0;
      P = ePa(0) + s * (ePa(1) + s * (ePa(2) + s * (ePa(3) + s * (ePa(4) + s * (ePa(5) + s * ePa(6))))));
      Q = 1.0 + s * (eQa(0) + s * (eQa(1) + s * (eQa(2) + s * (eQa(3) + s * (eQa(4) + s * eQa(5))))));
      retval = e_erx + P / Q;
    } else if (abs_x >= 6.0)
      retval = 1.0;
    else // 1.25 < |x| < 6.0 
      retval = 1.0 - erfc(abs_x);
    return if (x >= 0.0) retval else -retval
  }
}

object CumulativeNormalDistribution extends Function1D {
  /**
   * tip: transform to -x for x big as there is more precision around 0 than around 1.
   */
  def value(x: Double): Double = {
    return CodyCND.pnorms(x);
  }

}
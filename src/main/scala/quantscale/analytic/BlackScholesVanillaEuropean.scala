package quantscale.analytic

import quantscale.fdm.Epsilon
import quantscale.math.AS241InvCND

case class Price() extends BSMMeasure(Double.NaN)
case class DeltaGreek() extends BSMMeasure(Double.NaN)
case class GammaGreek() extends BSMMeasure(Double.NaN)
case class VegaGreek() extends BSMMeasure(Double.NaN)
case class ThetaGreek() extends BSMMeasure(Double.NaN)
case class RhoGreek() extends BSMMeasure(Double.NaN)
case class Rho2Greek() extends BSMMeasure(Double.NaN)

class BSMMeasure(var value: Double = Double.NaN) {
}

class DiscountFactor(val rate: Double, val time: Double) {
  def value = math.exp(-rate*time)
}

class Variance(val vol : Double, val time: Double) {
  def value = vol*vol*time
}

object BlackScholesVanillaEuropean {

  def df(rate: Double, time: Double): Double = {
    return Math.exp(-rate * time);
  }

  def variance(vol: Double, time: Double): Double = {
    return vol * vol * time;
  }

  def priceEuropeanVanilla(
    isCall: Boolean,
    strike: Double,
    spot: Double,
    variance: Double,
    driftDf: Double,
    discountDf: Double): Double = {

    val sign = if (isCall) 1 else -1;
    val sqrtVar = math.sqrt(variance);
    val forward = spot / driftDf;
    val d1 = 1.0 / sqrtVar * math.log(forward / strike) + 0.5 * sqrtVar;
    val d2 = d1 - sqrtVar;
    val nd1 = CumulativeNormalDistribution.value(sign * d1);
    val nd2 = CumulativeNormalDistribution.value(sign * d2);
    val price = sign * discountDf * (forward * nd1 - strike * nd2);
    return price;
  }

  def priceEuropeanVanillaAdjoint(isCall: Boolean,
                                  strike: Double,
                                  spot: Double,
                                  variance: Variance,
                                  driftDf: DiscountFactor,
                                  discountDf: DiscountFactor,
                                  measures: Array[BSMMeasure]
                                   ) {

    val vm0 = spot
    val vm1 = variance.vol
    val vm2 = variance.time
    val vm3 = driftDf.rate
    val vm4 = driftDf.time
    val vm5 = discountDf.rate
    val vm6 = discountDf.time

    val sign = if (isCall) 1 else -1
    val sqrtvm2 = math.sqrt(vm2)
    val sqrtVar = vm1*sqrtvm2
    val cf = math.exp(vm3*vm4)
    val forward = vm0 * cf
    val log = math.log(forward / strike)
    val d1 = 1.0 / sqrtVar * log + 0.5 * sqrtVar
    val d2 = d1 - sqrtVar
    val cnd1 = CumulativeNormalDistribution.value(sign * d1)
    val cnd2 = CumulativeNormalDistribution.value(sign * d2)
    val discountFactor = math.exp(-vm5*vm6)
    val price = sign * discountFactor * (forward * cnd1 - strike * cnd2)

    val priceBar = 1.0
    val discountBar = priceBar*sign*(forward * cnd1 - strike * cnd2)
    val vm5Bar = -discountBar*vm6*discountFactor
    val vm6Bar = -discountBar*vm5*discountFactor
    val cnd1Bar = priceBar*sign*discountFactor*forward
    val cnd2Bar = -priceBar*sign*discountFactor*strike
    val d2Bar = cnd2Bar * sign * NormalDistribution.value(sign*d2)
    val d1Bar = cnd1Bar * sign * NormalDistribution.value(sign*d1) +d2Bar
    val sqrtVarBar = -d2Bar + d1Bar*(0.5-log/(sqrtVar*sqrtVar))
    val logBar = d1Bar/sqrtVar
    val forwardBar =logBar/forward + priceBar*sign * discountFactor * cnd1
    val vm1Bar = sqrtVarBar*sqrtvm2
    val sqrtvm2Bar = sqrtVarBar*vm1
    val vm2Bar = sqrtvm2Bar*0.5/sqrtvm2
    val vm0Bar = forwardBar*cf
    val cfBar = forwardBar*vm0
    val vm3Bar = cfBar*vm4*cf
    val vm4Bar = cfBar*vm3*cf


    for (measure <- measures) {
      measure match {
        case Price() =>
          measure.value = price
        case DeltaGreek() =>
          measure.value = vm0Bar
        case VegaGreek() =>
          measure.value = vm1Bar
        case RhoGreek() => measure.value = vm5Bar
        case Rho2Greek() => measure.value = vm3Bar
        case ThetaGreek() => measure.value = vm2Bar + vm4Bar + vm6Bar
      }
    }
  }
  def priceEuropeanVanilla(
    isCall: Boolean,
    strike: Double,
    spot: Double,
    variance: Double,
    driftDf: Double,
    discountDf: Double,
    measures: Array[BSMMeasure]) {

    val sign = if (isCall) 1 else -1;
    val sqrtVar = math.sqrt(variance);
    val forward = spot / driftDf;
    val d1 = 1.0 / sqrtVar * math.log(forward / strike) + 0.5 * sqrtVar;
    val d2 = d1 - sqrtVar;
    val cnd1 = CumulativeNormalDistribution.value(sign * d1);
    val cnd2 = CumulativeNormalDistribution.value(sign * d2);

    for (measure <- measures) {
      measure match {
        case Price() =>
          measure.value = sign * discountDf * (forward * cnd1 - strike * cnd2)
        case GammaGreek() =>
          val nd1 = NormalDistribution.value(d1)
          measure.value = nd1 * discountDf / (driftDf * spot * sqrtVar)
        case DeltaGreek() =>
          measure.value = discountDf / driftDf * (if (isCall) cnd1 else (cnd1 - 1))
      }
    }
  }
  
}
//
//class JaeckelBlackImpliedVolatility {
//        
//        val tolerance = Math.pow(Epsilon.MACHINE_EPSILON, 0.75);
//        val maxNumberOfIterations = 1000;
//        val useHalleysMethod = true;
//    
//    /**
//     * Computes the implied volatility of Black model given the normalised premium (notional of 1).
//     * 
//     * @param prem - premium of option on a notional of 1
//     * @param F - forward price of underlying
//     * @param K - strike of underlying
//     * @param T - time to expiry of option
//     * @param isCall - is option a call option
//     * @param df - discounting factor
//     *  
//     * @return implied Black vol
//     */
//    def impliedVolatility(prem : Double,
//                                                                F : Double,
//                                                                K:Double,
//                                                                T:Double,
//                                                                 isCall : Boolean,
//                                                                 df: Double) : Double = 
//    {
//        var theta = if (isCall)  1 else -1;
//        
//        // Log moneyness & normalized target option premium 
//        val x = Math.log(F/K); 
//        var beta=prem/(df*Math.sqrt(F*K));
//        val sqrT = Math.sqrt(T);
//        
//        // Closed-form for ATM options
//        if(Math.abs(x)<Epsilon.MACHINE_EPSILON_SQRT) {
//               val sigma = -2*AS241InvCND.value(0.5*(1.-beta));
//                return sigma/sqrT;
//        }
//        
//        // Operate on Out-Of-The-Money
//        if(isCall && x>0 || !isCall && x<0) {
//                beta -= normalizedIntrinsicValue(x, theta);
//                theta = 1-2*H(x);
//        }
//        val intrinsicValue=normalizedIntrinsicValue(x, theta);
//        
//        // Vol and normalized option value at the inflection point
//        val sigma_c = Math.sqrt(2*Math.abs(x));
//        val b_c = normalizedBlackPremium(sigma_c, x, theta);
//        
//        // Initial guess
//        var sigma = initialGuess(x, beta, b_c, theta);
//        var b = normalizedBlackPremium(sigma, x, theta);
//        var v = objectiveFunction(x, intrinsicValue, beta, b, b_c, theta);
//        
//        var i=0;
//        while(Math.abs(v)>tolerance && i<maxNumberOfIterations) {
//                 val Dsigma = useHalleysMethod ? 
//                                 iterationStepUsingHalleysMethod(x, intrinsicValue, beta, b, b_c, sigma, theta) : 
//                                 iterationStep(x, intrinsicValue, beta, b, b_c, sigma, theta);
//                 sigma += Dsigma;
//                 b = normalizedBlackPremium(sigma, x, theta);
//                 v = objectiveFunction(x, intrinsicValue, beta, b, b_c, theta);
//                 i+=1
//        }
//        return if (i<maxNumberOfIterations) sigma/sqrT  else Double.NaN 
//    }
//    
//    /**
//     * Computes the normlaizedblack option value.
//     * @param sigma - the term Black vol: vol*timeToExpiry
//     * @return premium of european option evaluated using Black's model
//     */
//    private def normalizedBlackPremium( sigma : Double, x: Double, theta: Double) : Double =
//    {   
//        if(Math.abs(x-0.)<Epsilon.MACHINE_EPSILON_SQRT) {
//                return 1.-2 * CumulativeNormalDistribution.value(-0.5*sigma);  // At-The-Money
//        }
//        return theta*Math.exp(0.5*x)*CumulativeNormalDistribution.value(theta*(x/sigma+0.5*sigma))-
//                                theta*Math.exp(-0.5*x)*CumulativeNormalDistribution.value(theta*(x/sigma-0.5*sigma));   
//    }
//    
//    /**
//     * Computes the normalized intrinsic value.
//     * 
//     * @param x - the log moneyness
//     * @param theta - is 1 for a call option, -1 for a put option  
//     * @return normalized intrinsic value
//     */
//    private def normalizedIntrinsicValue(x: Double, theta: Int) : Double =
//    {
//        return if (theta * x >= 0)  theta * (Math.exp(0.5*x)-Math.exp(-0.5*x))  else 0.;
//    }
//    
//    private def sigmaLow(x:Double, beta:Double,  b_c:Double, theta:Int) : Double =  {
//        val i = normalizedIntrinsicValue(x, theta);
//        val t1=Math.log((beta-i)/(b_c-i));
//        val t2=Math.abs(x)-4*t1;
//        return Math.sqrt((2*x*x)/t2);
//    }
//    
//    private def sigmaHigh(x:Double, beta:Double,  b_c:Double, theta:Int) : Double = {
//        val t1 = CumulativeNormalDistribution.value(-Math.sqrt(0.5*Math.abs(x)));
//        val t2 = Math.exp(0.5*theta*x);
//        val t3 = (t2-beta)/(t2-b_c);
//        return -2 * AS241InvCND.value(t3*t1);
//    }
//    
//    private def gamma(x:Double, b_c:Double, theta:Int) : Double = {
//        val sigmaStar = sigmaHigh(x, 0., b_c, theta);
//        val bStar = normalizedBlackPremium(sigmaStar, x, theta);
//        val sigmaLowStar = sigmaLow(x, bStar, b_c, theta);
//        val sigmaHighStar = sigmaHigh(x, bStar, b_c, theta);
//        val t2 = (sigmaStar-sigmaLowStar)/(sigmaHighStar-sigmaLowStar);
//        if(t2<0.) {
//                return 0.;
//        }
//        return Math.log(t2)/Math.log(bStar/b_c);
//    }
//    
//    private double interpolatedSigmaLow(double x, double beta, double b_c, int theta) {
//        final double normalizedBeta = beta/b_c;
//        final double w = Math.min(Math.pow(normalizedBeta, gamma(x, b_c, theta)), 1.);
//        return (1.-w)*sigmaLow(x, beta, b_c, theta)+w*sigmaHigh(x, beta, b_c, theta);
//    }
//    
//    private double initialGuess(double x, double beta, double b_c, int theta) {
//        return beta<=b_c ? interpolatedSigmaLow(x, beta, b_c, theta) :  sigmaHigh(x, beta, b_c, theta);
//    }
//    
//    private double objectiveFunction(double x, double intrinsicValue, double targetOptionPrice, double currentOptionPrice, double b_c, int theta) {
//        return targetOptionPrice<=b_c ? 
//                        1./Math.log(currentOptionPrice-intrinsicValue) - 1./Math.log(targetOptionPrice-intrinsicValue) :
//                        currentOptionPrice - targetOptionPrice;
//    }
//    
//    private double iterationStep(double x, double intrinsicValue, double targetOptionValue, double b, double b_c, double sigma, int theta)
//    {
//        final double bDerivative = bPrime(x, sigma);
//        if(targetOptionValue>=b_c) {
//                return (targetOptionValue-b)/bDerivative;
//        }
//        else {
//                return Math.log((targetOptionValue-intrinsicValue)/(b-intrinsicValue))*
//                           (Math.log(b-intrinsicValue)/Math.log(targetOptionValue-intrinsicValue))*(b-intrinsicValue)/bDerivative;
//        }
//    }
//    
//    private double iterationStepUsingHalleysMethod(double x,
//                                                                                           double intrinsicValue,
//                                                                                           double targetOptionValue,
//                                                                                           double b,
//                                                                                           double b_c,
//                                                                                           double sigma,
//                                                                                           int theta)
//    {
//        double fDoublePrimeOverfPrime = bDoublePrimeOverbPrime(x, sigma);
//        if(targetOptionValue<b_c) {
//                final double logBetaMinusi = Math.log(b-intrinsicValue);
//                fDoublePrimeOverfPrime -= ((2+logBetaMinusi)/logBetaMinusi)*bPrime(x, sigma)/(b-intrinsicValue);
//        }
//        final double nu_n = iterationStep(x, intrinsicValue, targetOptionValue, b, b_c, sigma, theta);
//        final double nuHat_n = Math.max(nu_n, -0.5*sigma);
//        final double etaHat_n = Math.max(0.5*nuHat_n*fDoublePrimeOverfPrime, -0.75);
//        return Math.max(nuHat_n/(1.+etaHat_n), -0.5*sigma);
//    }
//    
//    private double bPrime(double x, double sigma)
//    {
//        final double sigmaSquared = sigma*sigma;
//        return Math.exp(-0.5*x*x/sigmaSquared-0.125*sigmaSquared)*MathConstants.M_ONE_DIV_SQRT_TWO_PI;
//    }
//    
//    private double bDoublePrimeOverbPrime(double x, double sigma)
//    {
//        return x*x / (sigma*sigma*sigma) - 0.25*sigma;
//    }
// 
//    // Heaviside function; should be somewhere more centrally
//    private int H(double x) {
//        return x>=0 ? 1 : 0; 
//    }
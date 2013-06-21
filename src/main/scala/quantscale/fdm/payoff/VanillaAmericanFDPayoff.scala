package quantscale.fdm.payoff
import quantscale.fdm.Epsilon

class VanillaAmericanFDPayoff(
  isCall: Boolean,
  strike: Double,
  firstExerciseTime: Double,
  expiryTime: Double) extends AmericanFDPayoff {

  override def copy(): FDPayoff = {
    return new VanillaAmericanFDPayoff(isCall, strike, firstExerciseTime, expiryTime)
  }

  private val putCallSign: Int = if (isCall) 1 else -1

  private var lowerBoundV: Array[Double] = null

  var exerciseBoundary: BoundaryListener = null

  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous

  override def initState(underlyingLevel: Array[Double]) {
    super.initState(underlyingLevel);
    lowerBoundV = null;
  }

  private def computeIntrinsicValue(out: Array[Double]) {
    var i: Int = 0;
    while (i < state.size) {
      out(i) = math.max(putCallSign * (space(i) - strike), 0);
      i += 1;
    }
  }

  def lowerBound: Array[Double] = {
    if (lowerBoundV == null) {
      lowerBoundV = new Array[Double](space.length);
      computeIntrinsicValue(lowerBoundV);
    }
    return lowerBoundV;
  }

  override def eval(): Unit = {
    if (math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      computeIntrinsicValue(state.price);
      _isDiscontinuous = true
    } else if (time > firstExerciseTime - Epsilon.MACHINE_EPSILON_SQRT) {
      var i: Int = 0;
      lowerBoundV = lowerBound;
      var iBoundary = 0
      _isDiscontinuous = false
      var isAbove = 0
      while (i < state.size) {
        val currentLowerBound = lowerBoundV(i)
        if (currentLowerBound > state.price(i)) {
          state.price(i) = currentLowerBound
          if (currentLowerBound > state.price(i) + Epsilon.MACHINE_EPSILON_SQRT) {
            //          if (isAbove == 0) isAbove = 1 else if (isAbove> 0) {iBoundary = i}  
            _isDiscontinuous = true
          }
        } else {
          if (currentLowerBound < state.price(i) - Epsilon.MACHINE_EPSILON_SQRT) {
            if (isAbove == 0) {
              isAbove = 1
              iBoundary = i
            }
          }
        }
        //        state(i) = math.max(lowerBoundV(i), state(i));
        i += 1;
      }
      //we know that the exercise boundary is between iBoundary and iBoundary+1
      var guess = space(iBoundary)
      if (exerciseBoundary != null && iBoundary > 1 && iBoundary < space.length - 1) {
        if (isAbove > 0) {
          //find a and b in a*x+b
          val aBound = (lowerBound(iBoundary) - lowerBound(iBoundary - 1)) / (space(iBoundary) - space(iBoundary - 1))
          val bBound = lowerBound(iBoundary) - aBound * space(iBoundary)
          val aState = (state.price(iBoundary + 1) - state.price(iBoundary)) / (space(iBoundary + 1) - space(iBoundary))
          val bState = state.price(iBoundary) - aState * space(iBoundary)
          guess = -(bState - bBound) / (aState - aBound)
          var midStart = (space(iBoundary) + space(iBoundary + 1)) / 2
          var midEnd = 0.
          if (guess > midStart) {
            iBoundary += 1
          }
          midStart = (space(iBoundary) + space(iBoundary - 1)) / 2
          midEnd = (space(iBoundary + 1) + space(iBoundary)) / 2
          var average = (0.5 * aBound * (guess * guess - midStart * midStart) + bBound * (guess - midStart) +
            0.5 * aState * (midEnd * midEnd - guess * guess) + bState * (midEnd - guess)) / (midEnd - midStart)
          //slight improvement, but not great as it creates a spike in the gamma - maybe go to a quadratic fit
          //          state(iBoundary) = average
        } else if (isAbove < 0) {
          val aBound = (lowerBound(iBoundary + 1) - lowerBound(iBoundary)) / (space(iBoundary + 1) - space(iBoundary))
          val bBound = lowerBound(iBoundary) - aBound * space(iBoundary)
          val aState = (state.price(iBoundary) - state.price(iBoundary - 1)) / (space(iBoundary) - space(iBoundary - 1))
          val bState = state.price(iBoundary) - aState * space(iBoundary)
          guess = -(bState - bBound) / (aState - aBound)
          var midStart = (space(iBoundary) + space(iBoundary + 1)) / 2
          var midEnd = 0.
          if (guess > midStart) {
            iBoundary += 1
          }
          midStart = (space(iBoundary) + space(iBoundary - 1)) / 2
          midEnd = (space(iBoundary + 1) + space(iBoundary)) / 2
          var average = (0.5 * aState * (guess * guess - midStart * midStart) + bState * (guess - midStart) +
            0.5 * aBound * (midEnd * midEnd - guess * guess) + bBound * (midEnd - guess)) / (midEnd - midStart)
          //          state(iBoundary) = average
        }
      }
      if (exerciseBoundary != null) exerciseBoundary.setValue(time, guess, state.price(iBoundary))
    } else {
      _isDiscontinuous = false
    }
  }

}
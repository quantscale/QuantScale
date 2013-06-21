package quantscale.fdm.payoff
import quantscale.fdm.Epsilon
import java.util.Arrays

class VanillaBermudanFDPayoff(
  isCall: Boolean,
  strike: Double,
  exerciseTime: Array[Double]) extends AmericanFDPayoff {

      override def copy(): FDPayoff = {
        return new VanillaBermudanFDPayoff(isCall, strike, exerciseTime)
      }
  private val putCallSign: Int = if (isCall) 1 else -1;

  private var lowerBoundV: Array[Double] = null;

  private val expiryTime = exerciseTime.last

  private var zeroBoundV: Array[Double] = null;

  var exerciseBoundary: BoundaryListener = null

  private var _isDiscontinuous = false

  override def isDiscontinuous = _isDiscontinuous

  override def initState(underlyingLevel: Array[Double]) {
    super.initState(underlyingLevel);
    lowerBoundV = null;
    zeroBoundV = null;
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
      //      println("C "+Arrays.toString(lowerBoundV));
    }
    if (zeroBoundV == null) {
      zeroBoundV = new Array[Double](space.length);
    }
    //      println("L "+Arrays.toString(lowerBoundV));
    if (isExerciseTime) {
      //      println("E "+Arrays.toString(lowerBoundV));
      return lowerBoundV;
    } else {
      return zeroBoundV
    }

  }

  def isExerciseTime(): Boolean = {
    var i = exerciseTime.length - 1
    while (i >= 0) {
      if (math.abs(time - exerciseTime(i)) < Epsilon.MACHINE_EPSILON_SQRT) {
        //        println("found exercise time "+time)
        return true
      }
      i -= 1
    }
    return false
  }

  override def eval(): Unit = {
    if (math.abs(expiryTime - time) < Epsilon.MACHINE_EPSILON_SQRT) {
      computeIntrinsicValue(state.price);
      _isDiscontinuous = true
    } else {
      val _lowerBound = lowerBound;
      var iBoundary = 0
      var i = 0
      _isDiscontinuous = false
      while (i < state.size) {
        if (_lowerBound(i) > state.price(i)) {
          state.price(i) = _lowerBound(i)
          if (iBoundary == 0) iBoundary = i
          if (_lowerBound(i) > state.price(i) + Epsilon.MACHINE_EPSILON_SQRT) _isDiscontinuous = true
        }
        i += 1;
      }
      if (exerciseBoundary != null) exerciseBoundary.setValue(time, space(iBoundary), state.price(iBoundary))
    }
  }

}
package quantscale.fdm.payoff;
import quantscale.fdm.transform.CoordinateTransformation
import quantscale.fdm.transform.IdentityTransformation
import quantscale.fdm.State

abstract class FDPayoff() {
    var state : State = null; //state can have n dimensions
    var space : Array[Double] = null;
    var originalSpace : Array[Double] = null;
    var time : Double = 0;
    
    var spaceTransform : CoordinateTransformation = IdentityTransformation;
    
    def copy() : FDPayoff
    
    def stateDimensions(): Int = {return 1}
    
    def initState(underlyingLevel : Array[Double]) {
        state = new State(stateDimensions, underlyingLevel.size)
        originalSpace = underlyingLevel;
        space = spaceTransform.transform(underlyingLevel);
    }
    
    def setTime(t:Double) {
      time = t;
    }
    
    /**
     * Done backward in time
     */
    def eval(); 
    
    /**
     * discontinuity at current time?
     */
    def isDiscontinuous = false
    
}
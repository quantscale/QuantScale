package quantscale.fdm
import org.slf4j.LoggerFactory
import java.util.Arrays
import quantscale.fdm.payoff.FDPayoff
import quantscale.fdm.payoff.FDPayoffSmoother
import quantscale.fdm.listener.FDMListener1D
import quantscale.fdm.mesh.MeshUtil
import quantscale.fdm.method.Parabolic1DMethod

class FDMSolver1D(val spec: Parabolic1DFDSpec, methodV: Parabolic1DMethod, val solver: TridiagonalSolver) {
  final val logger = LoggerFactory.getLogger(getClass());

  var smoother: FDPayoffSmoother = null;
  var listener: FDMListener1D = null;

  var grid = spec.grid;
  val method = methodV;
  var price: Array[Double] = null;
  
  def copy() : FDMSolver1D = {
    var c = new FDMSolver1D(spec.copy(), methodV.copy(), solver.copy())
    c.smoother = smoother
    return c
  }
  
  //init method matrix with grid, what about mu, sigma, r?
  def solve(payoff: FDPayoff) {

    solver.init(grid.spaceSize);
    method.solver = solver;
    var timeIterator = grid.timeIterator; //could be adaptive with error guess
    var currentTime = timeIterator.next();
    var previousTime = currentTime;
    //init price vector
    payoff.spaceTransform = grid.spaceTransform;
    payoff.initState(grid.spaceVector); //returns the payoff state
    payoff.setTime(currentTime);
    payoff.eval(); //smoothing in the payoff, eventually

    if (smoother != null) {
      smoother.makeSmooth(payoff);
    }
    if (logger.isDebugEnabled()) {
      logger.debug("currentTime=" + currentTime);
      logger.debug("payoff state=" + payoff.state)
      logger.debug("x=" + Arrays.toString(grid.spaceVector));
    }
    method.initSystem(spec);
    while (timeIterator.hasNext) {
      currentTime = timeIterator.next();
      if (logger.isDebugEnabled()) {
        logger.debug("currentTime=" + currentTime);
      }
      
      payoff.setTime(currentTime);
      method.solve(currentTime, previousTime - currentTime, payoff.state); //new price in price
     
      payoff.eval(); //eventually update price (barrier, etc.)
      if (logger.isDebugEnabled()) {
        logger.debug("payoff state=" + payoff.state)
      }
      if (listener != null) listener.update(currentTime, payoff.state, spec.grid.spaceVector);
      previousTime = currentTime;
    }
    method.shutdown()
    price = payoff.state.price;
  }

  def price(level: Double): Double = {
    if (price == null) {
      throw new RuntimeException("price is undefined, use solve before");
    }
    //interpolate to get price@level.
    val i = MeshUtil.findIndex(grid.spaceVector, level);
    return MeshUtil.interpolateLinearly(grid.spaceVector, price, level, i);
  }
}
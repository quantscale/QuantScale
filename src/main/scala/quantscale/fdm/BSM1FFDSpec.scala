package quantscale.fdm

abstract class BSM1FFDSpec(gridV : Grid) {
   val grid : Grid = gridV;
   
   //discretized rates, eventually.
   
   def vol(timeEnd: Double, dt: Double, space: Double) : Double;
   def drift(timeEnd: Double, dt: Double): Double;
   def discountRate(timeEnd: Double, dt: Double): Double;
}
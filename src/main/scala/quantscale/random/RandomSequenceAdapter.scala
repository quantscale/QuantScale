package quantscale.random

class RandomSequenceAdapter(rng: RandomGeneratorDouble, val dimension : Int) extends RandomSequenceGenerator {
  
  
   def nextSequence(sequence : Array[Double]) {
     var i = 0
     while (i < dimension) {
       sequence(i) = rng.nextDoubleOpen()
       i+=1
     }
   }

}
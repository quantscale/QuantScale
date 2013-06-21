package quantscale.mc

class ConstantDiscounter(val rate : Double) extends Discounter {
  //time is the monte carlo time
    def df(time : Double) : Double = {
        return math.exp(-time*rate)
    }
}
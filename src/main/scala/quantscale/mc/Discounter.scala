package quantscale.mc

trait Discounter {
    def df(time : Double) : Double
}
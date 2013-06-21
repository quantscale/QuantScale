package quantscale.random

trait RandomSequenceFactory {
  def makeRandomSequenceGenerator(dimension: Int): RandomSequenceGenerator
}

class PseudoRandomSequenceFactory(rng: RandomGeneratorDouble) extends RandomSequenceFactory {
  def makeRandomSequenceGenerator(dimension: Int): RandomSequenceGenerator = {
    return new RandomSequenceAdapter(rng, dimension)
  }
}
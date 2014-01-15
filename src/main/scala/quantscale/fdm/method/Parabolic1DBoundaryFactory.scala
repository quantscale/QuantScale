package quantscale.fdm.method

import quantscale.fdm.OperatorLine

trait ParabolicODEBoundaryFactory {
  def makeODELine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, line: OperatorLine)
}

trait Parabolic1DBoundaryFactory extends ParabolicODEBoundaryFactory {
  def makeLine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, multiplier: Double, line: OperatorLine)
}

object ForwardLinearOrder1Parabolic1DBoundaryFactory extends Parabolic1DBoundaryFactory with ParabolicODEBoundaryFactory {
  def makeLine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, multiplier: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD1Forward(i, x, b * multiplier).plus(i, 1.0 + multiplier * c)
  }

  def makeODELine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD1Forward(i, x, b).plus(i, c)
  }
}

object ForwardPartialOrder2Parabolic1DBoundaryFactory extends Parabolic1DBoundaryFactory {
  def makeLine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, multiplier: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD2Forward(i, x, a * multiplier)
      .plusD1Forward2(i, x, b * multiplier)
      .plus(i, 1.0 + multiplier * c)
  }

  def makeODELine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD2Forward(i, x, a)
      .plusD1Forward2(i, x, b)
      .plus(i, c)
  }

}

object BackwardLinearOrder1Parabolic1DBoundaryFactory extends Parabolic1DBoundaryFactory {
  def makeLine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, multiplier: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD1Backward(i, x, b * multiplier).plus(i, 1.0 + multiplier * c)
  }

  def makeODELine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD1Backward(i, x, b).plus(i, c)

  }
}

object BackwardPartialOrder2Parabolic1DBoundaryFactory extends Parabolic1DBoundaryFactory {
  def makeLine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, multiplier: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD2Backward(i, x, a * multiplier)
      .plusD1Backward2(i, x, b * multiplier)
      .plus(i, 1.0 + multiplier * c)
  }

  def makeODELine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD2Backward(i, x, a)
      .plusD1Backward2(i, x, b)
      .plus(i, c)
  }
}

object BackwardOrder1Parabolic1DBoundaryFactory extends Parabolic1DBoundaryFactory {
  def makeLine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, multiplier: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD2Backward(i, x, a * multiplier)
      .plusD1Backward(i, x, b * multiplier)
      .plus(i, 1.0 + multiplier * c)
  }

  def makeODELine(i: Int, x: Array[Double], a: Double, b: Double, c: Double, line: OperatorLine) {
    line.fill(0.0)
    line.plusD2Backward(i, x, a)
      .plusD1Backward(i, x, b)
      .plus(i, c)
  }
}  
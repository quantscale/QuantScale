package quantscale.mc

import quantscale.random.RandomSequenceFactory
import quantscale.random.RandomSequenceGenerator
import quantscale.math.AS241InvCND
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.SingularValueDecomposition
import quantscale.math.SVDSqrt

class BSMMCSpec(val assets: Array[BSM1FMCSpec], val correlation: Array[Array[Double]]) {
  def this(asset: BSM1FMCSpec) {
    this(Array(asset), Array(Array(1.0)))
  }
}

trait BSMPathListener {
  def init(path: MCPath, sqrtCovariance: Array[RealMatrix])
  def startPath(x0: Array[Double])
  def nextStep(i: Int, j: Int, diffusion: Double)
}

class BSM1FMartingaleQuotePower(maxPower: Int, quoteIndex: Int) extends BSMPathListener {
  var martingalePowerPath: Array[MCPath] = null //power-1
  private var martingalePowerStart: Array[Double] = null //power-1,
  private var _sqrtCovariance: Array[RealMatrix] = null

  def init(path: MCPath, sqrtCovariance: Array[RealMatrix]) {
    _sqrtCovariance = sqrtCovariance
    val dimension = _sqrtCovariance(0).getRowDimension()
    martingalePowerPath = Array.ofDim[MCPath](maxPower)
    var p = 0
    while (p < maxPower) {
      martingalePowerPath(p) = new MCPath(path.time, path.dimension)
      p += 1
    }
    martingalePowerStart = Array.ofDim[Double](maxPower)
  }

  def startPath(x0: Array[Double]) {
    var p = 0
    while (p < martingalePowerStart.length) {
      martingalePowerStart(p) = (p + 1) * x0(0)
      p += 1
    }
  }

  def nextStep(i: Int, j: Int, diffusion: Double) {
    var p = 0
    while (p < maxPower) {
      var previousStep = if (i == 0) martingalePowerStart(p) else martingalePowerPath(p).values(i - 1)(j)
      var sq = _sqrtCovariance(i).getEntry(quoteIndex, quoteIndex)
      sq *= sq
      martingalePowerPath(p).values(i)(j) = previousStep - 0.5 * (p + 1) * (p + 1) * sq + diffusion * (p + 1)
      //TODO try diffusion*diffusion
      p += 1
    }
  }
}

class BSM2FMartingaleQuotePower(maxPower: Int) extends BSMPathListener {
  var martingalePowerPath: Array[MCPath] = null //2
  var martingalePowerPathCross: Array[MCPath] = null //1

  private var martingalePowerStart: Array[Array[Double]] = null //power-1, dimension
  private var martingalePowerStartCross: Array[Double] = null //power-1,

  private var _sqrtCovariance: Array[RealMatrix] = null

  def init(path: MCPath, sqrtCovariance: Array[RealMatrix]) {
    _sqrtCovariance = sqrtCovariance
    val dimension = _sqrtCovariance(0).getRowDimension()

    martingalePowerPath = Array.ofDim[MCPath](maxPower)
    var p = 0
    while (p < maxPower) {
      martingalePowerPath(p) = new MCPath(path.time, path.dimension)
      p += 1
    }
    martingalePowerPathCross = Array.ofDim[MCPath](1)
    martingalePowerPathCross(0) = new MCPath(path.time, 1)
    martingalePowerStart = Array.ofDim[Double](maxPower, path.dimension)
    martingalePowerStartCross = Array.ofDim[Double](1)
  }

  def startPath(x0: Array[Double]) {
    var p = 0
    while (p < martingalePowerStart.length) {
      var d = 0
      while (d < x0.length) {
        martingalePowerStart(p)(d) = (p + 1) * x0(d)
        d += 1
      }
      p += 1
    }
    martingalePowerStartCross(0) = x0(0) + x0(1)
  }

  def nextStep(i: Int, j: Int, diffusion: Double) {
    var p = 0
    while (p < maxPower) {
      val previousStep = if (i == 0) martingalePowerStart(p)(j) else martingalePowerPath(p).values(i - 1)(j)
      var sq = _sqrtCovariance(i).getEntry(j, j) //an interesting biased other choice is diffusion
      sq *= sq
      martingalePowerPath(p).values(i)(j) = previousStep - 0.5 * (p + 1) * (p + 1) * sq + diffusion * (p + 1)
      p += 1
    }
    if (j == 0) {
      val previousCross = if (i == 0) martingalePowerStartCross(0) else martingalePowerPathCross(0).values(i - 1)(0)
      val sig11 = _sqrtCovariance(i).getEntry(0, 0)
      val sig12 = _sqrtCovariance(i).getEntry(0, 1)
      val sig21 = _sqrtCovariance(i).getEntry(1, 0)
      val sig22 = _sqrtCovariance(i).getEntry(1, 1)
      martingalePowerPathCross(0).values(i)(0) = previousCross - 0.5 * (sig11 * sig11 + sig12 * sig12 + 2 * sig11 * sig21 + 2 * sig12 * sig22 + sig21 * sig21 + sig22 * sig22) + diffusion
    } else if (j == 1) {
      martingalePowerPathCross(0).values(i)(0) += diffusion
    }
  }
}

class BSMPathGenerator(spec: BSMMCSpec, rngFactory: RandomSequenceFactory, brownianTransformFactory: BrownianTransformFactory) extends PathGenerator {
  private var _rng: RandomSequenceGenerator = null
  private var _path: MCPath = null
  private var _z: Array[Double] = null
  private var _zTransformed: Array[Array[Double]] = null
  private var _dimension = spec.assets.length
  private var _timeDimension = 0
  private var _x0 = Array.ofDim[Double](_dimension)
  private var _a, _b: Array[Array[Double]] = null
  private var _brownianTransform: BrownianTransform = null
  private var _sqrtCovariance: Array[RealMatrix] = null
  var sqrtMatrix = new SVDSqrt

  var pathListener: BSMPathListener = null

  def init(evaluationTimes: Array[Double]) {
    val assets = spec.assets
    _timeDimension = evaluationTimes.length
    _rng = rngFactory.makeRandomSequenceGenerator(_timeDimension * _dimension)

    _path = new MCPath(evaluationTimes, _dimension)
    _z = Array.ofDim[Double](_rng.dimension)
    _zTransformed = Array.ofDim[Double](_timeDimension, _dimension)

    var i = 0
    while (i < _dimension) {
      _x0(i) = math.log(assets(i).initialPrice)
      i += 1
    }

    _a = Array.ofDim(_timeDimension, _dimension)
    _b = Array.ofDim(_timeDimension, _dimension)

    i = 0
    var previousTime = 0.0

    while (i < _timeDimension) {
      var j = 0
      val time = evaluationTimes(i)
      while (j < _dimension) {
        val logDriftCf = math.log(assets(j).driftCf(time) / assets(j).driftCf(previousTime))
        val variance = math.max(1e-8, assets(j).variance(time) - assets(j).variance(previousTime))
        _a(i)(j) = logDriftCf - 0.5 * variance
        _b(i)(j) = Math.sqrt(variance)
        j += 1
      }
      previousTime = time
      i += 1
    }

    i = 0
    previousTime = 0.0
    _sqrtCovariance = Array.ofDim[RealMatrix](_timeDimension)
    while (i < _timeDimension) {
      var j = 0
      val time = evaluationTimes(i)
      val m = new Array2DRowRealMatrix(_dimension, _dimension)
      while (j < _dimension) {
        var k = 0
        while (k <= j) {
          val covar = spec.correlation(j)(k) * _b(i)(j) * _b(i)(k)
          m.setEntry(j, k, covar)
          m.setEntry(k, j, covar)
          k += 1
        }
        j += 1
      }
      _sqrtCovariance(i) = sqrtMatrix.sqrt(m)
      previousTime = time
      i += 1
    }
    if (pathListener != null) pathListener.init(_path, _sqrtCovariance)
    if (brownianTransformFactory != null) {
      _brownianTransform = brownianTransformFactory.makeBrownianTransform(evaluationTimes)
    }
  }

  def nextPath(): MCPath = {
    _rng.nextSequence(_z)
    var i = 0
    while (i < _z.length) {
      _z(i) = AS241InvCND.value(_z(i))
      i += 1
    }
    //BB / PCA
    if (_brownianTransform == null) {
      i = 0
      while (i < _timeDimension) {
        System.arraycopy(_z, i * _dimension, _zTransformed(i), 0, _dimension)
        i += 1
      }
    } else {
      _brownianTransform.transform(_z, _dimension, _zTransformed)
    }

    // X = ln S, dX = (drift - .5*sig^2) dt + sig dW
    var previousPath = _x0
    i = 0
    if (pathListener != null) pathListener.startPath(_x0)
    while (i < _timeDimension) {
      val a = _a(i)
      val sqrt = _sqrtCovariance(i)
      val zi = _zTransformed(i)
      var j = 0
      while (j < _dimension) {
        var k = 0
        var b = 0.0
        while (k < _dimension) {
          b += sqrt.getEntry(j, k) * zi(k);
          k += 1
        }
        _path.values(i)(j) = previousPath(j) + a(j) + b
        if (pathListener != null) {
          pathListener.nextStep(i, j, b)
        }
        j += 1
      }
      previousPath = _path.values(i)
      i += 1
    }
    return _path
  }
}
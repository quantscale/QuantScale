package quantscale.mc

import quantscale.math.QRRegression

class GlassermanYuSimulation(pathGenerator: PathGenerator, pathPricer: LongstaffSchwartzPathPricer) {

  var isITMRegression = true

  var trainingMean = Double.NaN
  var lowMean = Double.NaN
  var upMean = Double.NaN
  
  private def isInTheMoney() = false

  def evaluateForwardBackward(numTraining: Int, numSimulations: Int): Double = {
    var i = 0

    var path: MCPath = null
    pathGenerator.init(pathPricer.evaluationTimes)
    //compute exercise flows
    var exerciseTimes = pathPricer.exerciseTimes
    var exerciseFlows = Array.ofDim[Double](numTraining, exerciseTimes.length)
    var independentVariablesDimension = pathPricer.independentVariablesDimension
    var independentVariables = Array.ofDim[Double](numTraining, exerciseTimes.length, independentVariablesDimension)
    while (i < numTraining) {
      path = pathGenerator.nextPath()
      pathPricer.eval(path, exerciseFlows(i), independentVariables(i))
      i += 1
    }

    //regress
    var regressionList = Array.ofDim[QRRegression](exerciseTimes.length)
    var indexExercise = exerciseTimes.length - 2
    var sumTraining = 0.0
    while (indexExercise >= 0) {
      var regression = new QRRegression()
      regressionList(indexExercise) = regression
      var numInTheMoney = 0
      if (isITMRegression) {
        i = 0
        while (i < numTraining) {
          if (exerciseFlows(i)(indexExercise) > 0) {
            numInTheMoney += 1
          }
          i += 1
        }
      } else {
        numInTheMoney = numTraining
      }
      var x = Array.ofDim[Double](numInTheMoney, independentVariablesDimension)
      var y = Array.ofDim[Double](numInTheMoney)

      var itmIndex = 0
      i = 0
      //eventually move following to pathPricer assuming exerciseFlow is the discounted value to previous ex date

      while (itmIndex < numInTheMoney) {
        if ((!isITMRegression) || (exerciseFlows(i)(indexExercise) > 0)) {
          x(itmIndex) = independentVariables(i)(indexExercise + 1) //indexExercise = now
          y(itmIndex) = exerciseFlows(i)(indexExercise + 1)
          itmIndex += 1
        }
        i += 1
      }
      regression.regress(x, y)

      //update exerciseFlows
      i = 0
      while (i < numTraining) {
        val immediateValue = exerciseFlows(i)(indexExercise)
        if (!isITMRegression || immediateValue > 0) {
          val continuationValue = regression.value(independentVariables(i)(indexExercise))
          if (continuationValue > immediateValue) {
            exerciseFlows(i)(indexExercise) = exerciseFlows(i)(indexExercise + 1)
          }
        } else {
          if (isITMRegression) { //probably not necessary
            exerciseFlows(i)(indexExercise) = exerciseFlows(i)(indexExercise + 1)
          }
        }
        if (indexExercise == 0) {
          sumTraining += exerciseFlows(i)(0)
        }
        i += 1

      }
      indexExercise -= 1
    }
    trainingMean = sumTraining / numTraining
//    println("GY calibration value = " + trainingMean)
    //eval forward remaining path and compute mean exercise value using those only
    i = 0
    var sumExercise = 0.0
    var sumPayments = 0.0
    var sumUpperBound = 0.0
    var tmpExerciseFlow = Array.ofDim[Double](exerciseTimes.length)
    var tmpIndependentVariable = Array.ofDim[Double](exerciseTimes.length, independentVariablesDimension)
    while (i < numSimulations) {
      path = pathGenerator.nextPath()
      sumPayments += pathPricer.eval(path, tmpExerciseFlow, tmpIndependentVariable)
      indexExercise = 0
      var isMinTimeFound = false
      var lowerBound = tmpExerciseFlow(exerciseTimes.length - 1)
      var M = 0.0
      var U = 0.0
      while ((indexExercise < exerciseTimes.length - 1)) {
        var currentExerciseValue = tmpExerciseFlow(indexExercise)
        if (!isITMRegression || currentExerciseValue > 0) {
          val continuationValue = regressionList(indexExercise).value(tmpIndependentVariable(indexExercise)) //indexExercise = now
          val estimatedFutureValue = regressionList(indexExercise).value(tmpIndependentVariable(indexExercise + 1))
          if (!isMinTimeFound && currentExerciseValue > continuationValue) {
            isMinTimeFound = true
            lowerBound = currentExerciseValue
          }

          U = math.max(U, currentExerciseValue - M)
          M += (estimatedFutureValue - continuationValue)
        }
        //        U = math.max(currentExerciseValue,M)
        indexExercise += 1
      }
      U = math.max(U, tmpExerciseFlow(exerciseTimes.length - 1) - M)

      sumUpperBound += U
      //      indexExercise = exerciseTimes.length - 2
      //      var exerciseValue = 0.0; //tmpExerciseFlow(exerciseTimes.length-1)
      //      while (indexExercise >= 0) {
      //        var currentExerciseValue = tmpExerciseFlow(indexExercise)
      //        if (currentExerciseValue > 0) {
      //          val continuationValue = regressionList(indexExercise).value(tmpIndependentVariable(indexExercise)) //indexExercise = now
      //          if (continuationValue > currentExerciseValue) {
      //            currentExerciseValue = tmpExerciseFlow(indexExercise + 1)
      //          } 
      //        } else {
      //          currentExerciseValue = tmpExerciseFlow(indexExercise + 1)
      //        }
      //        exerciseValue = currentExerciseValue
      //        indexExercise -= 1
      //      }
      sumExercise += lowerBound //lowerBound
      i += 1
    }
    upMean = (sumPayments+sumUpperBound) / numSimulations
//    println("GY upper bound = " + upMean)
    lowMean = (sumPayments + sumExercise) / numSimulations
    return lowMean
  }
}
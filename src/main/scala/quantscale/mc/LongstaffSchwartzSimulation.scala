package quantscale.mc

import quantscale.math.QRRegression

class LongstaffSchwartzSimulation(pathGenerator: PathGenerator, pathPricer: LongstaffSchwartzPathPricer) {
    
  var trainingMean : Double = Double.NaN
  
  def evaluateForwardBackward(numTraining: Int, numSimulations: Int): Double = {
    var i = 0
    
    var path: MCPath = null
    pathGenerator.init(pathPricer.evaluationTimes)
    //compute exercise flows
    var exerciseTimes = pathPricer.exerciseTimes
    var exerciseFlows = Array.ofDim[Double](numTraining, exerciseTimes.length)
    var independentVariablesDimension = pathPricer.independentVariablesDimension
    var independentVariables = Array.ofDim[Double](numTraining, exerciseTimes.length, independentVariablesDimension)
    var sumTraining = 0.0
    
    while (i < numTraining) {
      path = pathGenerator.nextPath()
      pathPricer.eval(path, exerciseFlows(i), independentVariables(i))
      sumTraining += exerciseFlows(i)(exerciseTimes.length-1)
      i += 1
    }
//     println("LS european calibration value = " + sumTraining / numTraining)

    //regress
    var regressionList = Array.ofDim[QRRegression](exerciseTimes.length)
    var indexExercise = exerciseTimes.length - 2
    
    sumTraining = 0.0
    while (indexExercise >= 0) {
      var regression = new QRRegression()
      regressionList(indexExercise) = regression
      var numInTheMoney = 0
      i = 0
      while (i < numTraining) {
        if (exerciseFlows(i)(indexExercise) > 0) {
          numInTheMoney += 1
        }
        i += 1
      }
      var x = Array.ofDim[Double](numInTheMoney, independentVariablesDimension)
      var y = Array.ofDim[Double](numInTheMoney)

      var itmIndex = 0
      i = 0
      //eventually move following to pathPricer assuming exerciseFlow is the discounted value to previous ex date

      while (itmIndex < numInTheMoney) {
        if (exerciseFlows(i)(indexExercise) > 0) {
          x(itmIndex) = independentVariables(i)(indexExercise) //indexExercise = now
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
        if (immediateValue > 0) {
          val continuationValue = regression.value(independentVariables(i)(indexExercise)) //indexExercise = now
          if (continuationValue > immediateValue) {
              exerciseFlows(i)(indexExercise) = exerciseFlows(i)(indexExercise + 1)
          }
        } else {
          exerciseFlows(i)(indexExercise) = exerciseFlows(i)(indexExercise + 1)
        }
         if (indexExercise == 0) {
          sumTraining += exerciseFlows(i)(0)
        }
        i += 1
      }
      indexExercise -= 1
    }
     trainingMean = sumTraining / numTraining
//        println("LS calibration value = " + trainingMean)
        
    //eval forward remaining path and compute mean exercise value using those only
    i = 0
    var sumExercise = 0.0
    var sumEuropean = 0.0
    var sumPayments = 0.0
    var tmpExerciseFlow = Array.ofDim[Double](exerciseTimes.length)
    var tmpIndependentVariable = Array.ofDim[Double](exerciseTimes.length, independentVariablesDimension)
    while (i < numSimulations) {
      path = pathGenerator.nextPath()
      sumPayments += pathPricer.eval(path, tmpExerciseFlow, tmpIndependentVariable)
      //the following could move to pathPricer as it is independent of i
      indexExercise = exerciseTimes.length - 2
      var exerciseValue = 0.0; //tmpExerciseFlow(exerciseTimes.length-1)
      while (indexExercise >= 0) {
        var currentExerciseValue = tmpExerciseFlow(indexExercise)
        if (currentExerciseValue > 0) {
          val continuationValue = regressionList(indexExercise).value(tmpIndependentVariable(indexExercise)) //indexExercise = now
          if (continuationValue > currentExerciseValue) {
            tmpExerciseFlow(indexExercise) = tmpExerciseFlow(indexExercise + 1)
          } 
        } else {
          tmpExerciseFlow(indexExercise) = tmpExerciseFlow(indexExercise + 1)
        }
        exerciseValue = tmpExerciseFlow(indexExercise)
        indexExercise -= 1
      }
      sumExercise += exerciseValue
      sumEuropean+= tmpExerciseFlow(exerciseTimes.length-1)
      i += 1
    }
//    println("LS european value="+sumEuropean/numSimulations)
    return (sumPayments + sumExercise) / numSimulations
  }
}
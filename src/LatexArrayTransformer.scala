import scala.collection.mutable.ListBuffer
object LatexArrayTransformer extends App {

  val source = scala.io.Source.fromFile("scala_trbdf2vscrank.txt")
  val lines = source.getLines()
  var i = 0 //first col line index
  var j = 0 //second col line index
  var tSteps = 0
  val newLines = ListBuffer[String]()
  lines.foreach(line => {
    if (i == 0) {
      newLines += line
    } else {
      val columns = line.split(" ")
      //& TRBDF2 & 6.06077167 & 2.67E-02 & 0.53 \\ \hline
      val tStepsNew = Integer.parseInt(columns(0));
      if (tStepsNew >= 600) {
        j +=1
        var existingLine = newLines(j)
        existingLine = existingLine.stripSuffix("\\\\");
        if (tStepsNew != tSteps) {
          tSteps = tStepsNew
          var s = "& \\multicolumn{1}{|r|}{%d}".format(tSteps)
          s += " & %.7f &  %.2e & %.5f \\\\".format(
            columns(2).toDouble,
            columns(3).toDouble,
            columns(4).toDouble)
          existingLine += s
          
        } else {
          existingLine += "&& %.7f &  %.2e & %.5f \\\\".format(
            columns(2).toDouble,
            columns(3).toDouble,
            columns(4).toDouble)
        }
        newLines(j) = existingLine
      } else {
        if (tStepsNew != tSteps) {
          tSteps = tStepsNew
          var s = "\\hline \\multicolumn{1}{|r|}{%d}".format(tSteps)
          s += "& %s & %.7f &  %.2e & %.5f \\\\".format(
            columns(1).replace("_", "\\_"),
            columns(2).toDouble,
            columns(3).toDouble,
            columns(4).toDouble)
          newLines += s
        } else {
          newLines += "& %s & %.7f &  %.2e & %.5f \\\\".format(
            columns(1).replace("_", "\\_"),
            columns(2).toDouble,
            columns(3).toDouble,
            columns(4).toDouble)
        }
      }
    }
    i += 1
  })

  source.close()
  while (j<newLines.length) {
    var l = newLines(j).stripSuffix("\\\\")
    l += "& & & & \\\\"
    newLines(j) = l
    j+=1
  }
  newLines.foreach(line => println(line))

}
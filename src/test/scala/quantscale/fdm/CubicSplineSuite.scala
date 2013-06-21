package test.quantscale.fdm
import org.junit.runner.RunWith
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner
import quantscale.math.CubicSpline
import quantscale.math.NotAKnot

@RunWith(classOf[JUnitRunner])
class CubicSplineSuite extends FunSuite {

   test("NotAKnotAgainstScilab")  {
        val x = Array(-9.0, -5., -3., -2., 0., 1., 2., 4., 7., 9.)

        val y = Array(0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,  0.02, 0.0121951)

        val spline = CubicSpline.makeCubicSpline(x,y, new NotAKnot, new  NotAKnot)
        val z = new Array[Double](21);
        val expected = Array(0.0121951, -0.0436779, -0.0510208, -0.0253661,
            0.0177536,

            0.0628058, 0.0942579, 0.1001334, 0.2794847, 0.7838355, 1.,

            0.5539709, 0.2352764, 0.1163079, 0.0688725, 0.0487366,

            0.0347116, 0.0250999, 0.0188466, 0.0148966, 0.0121951 ) // values
        // obtained
        // from
        // Scilab.
        z(0) = x(0)
        for (i <- 1 until  21) {
            z(i) = z(i - 1) + (18 / 20.0);
            val value = spline.value(z(i));
            println( z(i) + " " + value + " " + expected(i));
            //        Log.debug(getLogCategory(), str);
            //        logger.info(str);
            //        console.writer().println(str);

            assert(Math.abs(expected(i)- value)< 1e-6,"matching scilab");

        }

    }

    test("NaturalExtrapolationAgainstScilab") {
        val x = Array( -9., -5., -3, -2, 0, 1, 2, 4, 7, 9 )

        val y = Array( 0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,

            0.02, 0.0121951 )

        val spline = CubicSpline.makeCubicSpline(x, y);

        val z = Array( -20., -10, -9, -8.1, 9, 10, 20 )

        val expected = Array (-1.1106312, 0.0192726, 0.0121951, 0.0056699,
            0.0121951, 0.0090887, -0.3722188 )

        // values obtained from Scilab.
        /*
         * deff('[y]=f(t)','y=ones(t)./(1+t.*t)'); t = [-9 -5 -3 -2 0 1 2 4 7
         * 9]; y=f(t); d = splin(t,y,"natural"); interp([-20 -10 -9 -8.1 9 10
         * 20],t,y,d, "natural")
         */

        for (i <- 0 until z.length) {
            val value = spline.value(z(i));
            println (z(i) + " " + value + " " + expected(i)+" "+math.abs(expected(i)-value))
            assert(math.abs(expected(i)- value)< 1e-5,"matching scilab");

        }

    }

//    public void testNotAKnotExtrapolationAgainstScilab()
//        throws InterpolatorException {
//        double[] x = { -9, -5, -3, -2, 0, 1, 2, 4, 7, 9 };
//
//        double[] y = { 0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,
//
//            0.02, 0.0121951 };
//
//
//        double z[] = { -20, -10, -9, -8.1, 9, 10, 20 };
//
//        double expected[] = { 10.566406, 0.1500863, 0.0121951, -0.0436779,
//            0.0121951,
//
//            0.0093745, -0.3240691 };
//
//        // values obtained from Scilab.
//        /*
//         * deff('[y]=f(t)','y=ones(t)./(1+t.*t)'); t = [-9 -5 -3 -2 0 1 2 4 7
//         * 9]; y=f(t); d = splin(t,y); interp([-20 -10 -9 -8.1 9 10 20],t,y,d,
//         * "natural")
//         */
//
//        for (int i = 0; i < z.length; i++) {
//            double value = spline.value(z[i]);
//            Log.debug(getLogCategory(), z[i] + " " + value + " " + expected[i]);
//            //        Log.debug(getLogCategory(), str);
//            //        logger.info(str);
//            //        console.writer().println(str);
//
//            assertEquals("matching scilab", expected[i], value, 1e-5);
//
//        }
//
//    }

    test("NaturalAgainstScilab") {
        val x = Array( -9., -5, -3, -2, 0, 1, 2, 4, 7, 9 )

        val y = Array( 0.0121951, 0.0384615, 0.1, 0.2, 1., 0.5, 0.2, 0.0588235,

            0.02, 0.0121951 )

        val spline = CubicSpline.makeCubicSpline(x, y);
        val z = new Array[Double](21)
        val expected = Array( 0.0121951, 0.0056699, 0.0031232, 0.0085338,
            0.0258802,

            0.0581919, 0.0914586, 0.1006235, 0.2793195, 0.7835630, 1.,

            0.5539767, 0.2352723, 0.1163160, 0.0688776, 0.0487282,

            0.0346897, 0.0250792, 0.0188565, 0.0149455, 0.0121951 ) // values
        // obtained
        // from
        // Scilab.
        z(0) = x(0);
        var maxError = 0.
        for (i <- 1 until 21) {
            z(i) = z(i - 1)+ (18 / 20.0);
            val value = spline.value(z(i));
            println( z(i) + " " + value + " " + expected(i))
            val newError = math.abs(value - expected(i));
            maxError = Math.max(maxError, newError);
            assert(math.abs(expected(i)-value)< 1e-6,"matching scilab");
        }
        println("maxError=" + maxError);
    }

}
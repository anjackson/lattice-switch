/**-------------------------------------------------------------
 * jBinLats - Utils.java
 * net.anjackson.maths.Utils
 * 
 * Created on 14-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.maths;

import java.lang.reflect.Method;

/**
 * This class contains some simple helper methods for numberical algorithms.
 * It take the place of the functions defined in the C header files of the 
 * Numerical Recipies code.
 *
 * @author ajackso1
 * @version $Id: Utils.java 628 2006-03-12 16:28:14Z anj $
 *
 */
public class Utils {

	/**
	 * Computes (a^2 + b^2)^1/2 without destructive underflow or overflow.
	 * @param a first parameter of Pythagoras eqn.
	 * @param b second parameter of Pythagoras eqn.
	 * @return (a^2 + b^2)^1/2
	 */
	public static double dpythag(double a, double b) {
		try {
			// Java 1.5+ has a builtin method for this: Math.hypot(a,b);
			Class mathc = (Class)Math.class;
			Method hypot = mathc.getDeclaredMethod("hypot", new Class[] {double.class, double.class});
			Double ans = (Double)hypot.invoke(mathc, new Object[] {new Double(a),new Double(b)});
			return ans.doubleValue();
		} catch (Exception e) {
			// If there was an exception, the method call failed.
			// TASK This is ridiculously complex for this problem.  Best to remove it for <1.5 compatibility.
		}
		// If we made it this far, we fall back on this safe calc:
		double absa,absb;
		absa=Math.abs(a);
		absb=Math.abs(b);
		if (absa > absb) {
			return absa*Math.sqrt(1.0+DSQR(absb/absa));
		} else {
			return (absb == 0.0 ? 0.0 : absb*Math.sqrt(1.0+DSQR(absa/absb)));
		}
	}

	public static double SIGN(double a, double b) {
		return ((b) >= 0.0 ? Math.abs(a) : -Math.abs(a));
	}

	public static float SQR(float a) {
		return (a == 0.0f ? 0.0f : a*a);
	}

	public static double DSQR(double a) {
		return (a == 0.0 ? 0.0 : a*a);
	}

	public static double DMAX(double a, double b) {
		return ((a) > (b) ? (a) : (b));
	}

	public static double DMIN(double a, double b) {
		return ((a) < (b) ? (a) : (b));
	}

	public static float FMAX(float a, float b) {
		return ((a) > (b) ? (a) : (b));
	}

	public static float FMIN(float a, float b) {
		return ((a) < (b) ? (a) : (b));
	}

	public static long LMAX(long a, long b) {
		return ((a) > (b) ? (a) : (b));
	}

	public static long LMIN(long a, long b) {
		return ((a) < (b) ? (a) : (b));
	}

	public static int IMAX(int a, int b) {
		return ((a) > (b) ? (a) : (b));
	}

	public static int IMIN(int a, int b) {
		return ((a) < (b) ? (a) : (b));
	}

}

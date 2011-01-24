/**-------------------------------------------------------------
 * jBinLats - StdinReader.java
 * net.anjackson.utils.StdinReader
 * 
 * Created on 13-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.utils;

import java.io.*;

/**
 * Simple wrapper and utility class for System.in.
 * a.k.a. stdin
 * 
 * The utility classes parse the input as various types, reporting errors 
 * when appropriate.
 * 
 * @author ajackso1
 * @version $Id: StdinReader.java 516 2006-01-25 18:17:55Z anj $
 *
 */
public class StdinReader {

	// The stream itself:
	static InputStreamReader raw_stdin = new InputStreamReader (System.in);
	static BufferedReader stdin = new BufferedReader (raw_stdin);

	// Read a String from standard system input
	public static String getString() {
		try {
			String line = stdin.readLine();
			//System.out.println("StdioReader: "+line);
			return line;
		} catch( Exception e ) {
			//System.out.println("getString() exception, returning empty string");
			return "";
		}
	}

	// Read a char from standard system input
	public static char getChar() {
		return ParseUtils.getChar(getString());
	}

	// Read a Number as a String from standard system input
	// Return the Number
	public static Number getNumber() {
		return ParseUtils.getNumber(getString());
	}

	// Read an integer from standard system input
	public static int getInt() {
		return ParseUtils.getInt(getString());
	}

	// Read a long integer from standard system input
	public static long getLong() {
		return ParseUtils.getLong(getString());
	}

	// Read a float from standard system input
	public static float getFloat() {
		return ParseUtils.getFloat(getString());
	}

	// Read a double from standard system input
	public static double getDouble() {
		return ParseUtils.getDouble(getString());
	}
	
	// Read an array of integers, space or comma seperated.
	public static int[] getIntegerArray() {
		return ParseUtils.getIntegerArray(getString());
	}
	
	/**
	 * Reads in a y/n or true/false or 1/0 and returns boolean true/false
	 * @return a boolean version of the input
	 */
	public static boolean getBoolean() {
		return ParseUtils.getBoolean(getString());
	}
}

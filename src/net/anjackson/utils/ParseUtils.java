/**-------------------------------------------------------------
 * jBinLats - ParseUtils.java
 * net.anjackson.utils.ParseUtils
 * 
 * Created on 09-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.utils;

import java.text.NumberFormat;

/**
 * Utility class for parsing strings into other types.
 * Any exceptions are caught and reported.
 * 
 * @author ajackso1
 * @version $Id: ParseUtils.java 525 2006-01-27 18:31:25Z anj $
 *
 */
public class ParseUtils {
	// Read a char from standard system input
	public static char getChar(String s) {
		if( s.length() >= 1 )
			return s.charAt(0);
		else
			return '\n';
	}

	// Read a Number as a String from standard system input
	// Return the Number
	public static Number getNumber(String numberString) {
		try {
			numberString = numberString.trim().toUpperCase();
			return NumberFormat.getInstance().parse( numberString );
		} catch( Exception e ) {
			// if any exception occurs, just give a 0 back
			System.out.println("getNumber() exception parsing '"+numberString+"', returning 0");
			return new Integer( 0 );
		}
	}

	// Read an integer from standard system input
	public static int getInt(String s) {
			return getNumber(s).intValue();
	}

	// Read a long integer from standard system input
	public static long getLong(String s) {
			return getNumber(s).longValue();
	}

	// Read a float from standard system input
	public static float getFloat(String s) {
			return getNumber(s).floatValue();
	}

	// Read a double from standard system input
	public static double getDouble(String s) {
			return getNumber(s).doubleValue();
	}
	
	// Read an array of integers, space or comma seperated.
	public static int[] getIntegerArray( String numberString ) {
		try {
			numberString = numberString.trim().toUpperCase();
			// Split the string on ' ' or ',':
			String[] intss = numberString.split("[\\s,]+");
			// Count the number of numbers:
			int[] ints = new int[intss.length];		
			// Put them into an integer array and return it:
			for( int i = 0; i < intss.length; i++ ) {
				ints[i] = Integer.parseInt(intss[i]);
			}
			return ints;
			
		} catch( Exception e ) {
			// if any exception occurs, just give a null back
			System.out.println("getIntegerArray() exception, returning null");
			return null;
		}
	}
	
	/**
	 * Reads in a y/n or true/false or 1/0 and returns boolean true/false
	 * @return a boolean version of the input
	 */
	public static boolean getBoolean( String booleanString ) {
		try {
			booleanString = booleanString.trim().toUpperCase();
			//System.out.println("StdioReader.getBoolean: Parsing "+booleanString);
			if( booleanString.length() == 1 ) {
				if( booleanString.charAt(0) == 'Y' || 
					booleanString.charAt(0) == '1' ) {
					return true;
				} else if( booleanString.charAt(0) == 'N' || 
					booleanString.charAt(0) == '0' ) {
					return false;
				}
				
			} else {
				if( "TRUE".equals(booleanString) ) {
					return true;
				} else if( "FALSE".equals(booleanString) ) {
					return false;
				}
			}
		} catch( Exception e ) {
			// if any exception occurs, just give 'false' back
			System.out.println("getBoolean() exception [returning false]: "+e);
			return false;
		}
		// If we get this far, the input could not be parsed as a boolean:
		System.out.println("getBoolean(): Could not parse input '"+booleanString+"' as boolean, returning false");
		return false;
	}
}

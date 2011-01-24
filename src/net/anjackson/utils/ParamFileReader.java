/**-------------------------------------------------------------
 * jBinLats - ParamFileReader.java
 * net.anjackson.utils.ParamFileReader
 * 
 * Created on 09-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.utils;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

/**
 * A parameter (.par) file is just the same as a normal java properties file,
 * but this interface simplifies reading the data in as doubles/floats etc.
 * 
 * Optional parameters can be checked with 'isDefined' and then specifying 
 * required = false in the getType method (no exception will be thrown then).
 * 
 * NOTE Could also specify default values, or do this instead of the required
 * flag (default == null implies required?).
 *
 * @author ajackso1
 * @version $Id: ParamFileReader.java 904 2006-07-11 17:27:10Z anj $
 *
 */
public class ParamFileReader extends Properties {
	// Whether to report the loaded parameters to stdout?
	boolean report_parameters = false;

	/**
	 * Getter for checking whether the parameters will be reported as they are 
	 * looked up.
	 * 
	 * @return boolean true if the parameters will be reported.
	 */
	public boolean getReportParameters() {
		return report_parameters;
	}
	
	/**
	 * Setter to define whether parameters will be reported as they are being
	 * looked up.
	 * 
	 * @param report true or false for reporting or not reporting the parameters.
	 */
	public void setReportParameters( boolean report ) {
		report_parameters = report;
	}
	
	/** 
	 * Simple routine to perform the basic file i/o.
	 * All you need to do is pass the filename.
	 * 
	 * @param filename of the properties file
	 */
	public void loadParamFile( String filename ) {
		try {
			this.load(new FileInputStream(filename));
		} catch( FileNotFoundException e ) {
			throw( new RuntimeException("ParamFileReader: File '"+filename+"' not found!", e));
		} catch( IOException e ) {
			throw( new RuntimeException("ParamFileReader: An I/O error occured while reading from '"+filename+"'.", e));
		}
	}
	
	/**
	 * Gets the parameter with the specified key and return it as a double.
	 * Throws an UndefinedParameterException if the parameter is null.
	 * 
	 * @param key the property to look up
	 * @return the value corresponding to key, as a double.
	 */
	public double getDouble( String key ) 
				throws UndefinedParameterException {
		String valueS = null;
		valueS = this.getPropertyAndCheck(key);
		double value = ParseUtils.getDouble(valueS);
		// Report to the user if required:
		if( report_parameters ) 
			System.out.println("# Parameter "+key+" has been given the double value "+value);
		return value;
	}

	/**
	 * Gets the parameter with the specified key and return it as a integer.
	 * Throws an UndefinedParameterException if the parameter is null.
	 * 
	 * @param key the property to look up
	 * @return the value corresponding to key, as a integer.
	 */
	public int getInteger( String key ) 
				throws UndefinedParameterException {
		String valueS = null;
		valueS = this.getPropertyAndCheck(key);
		int value = ParseUtils.getInt(valueS);
		// Report to the user if required:
		if( report_parameters ) 
			System.out.println("# Parameter "+key+" has been given the int value "+value);
		return value;
	}

	/**
	 * Gets the parameter with the specified key and return it as a long integer.
	 * Throws an UndefinedParameterException if the parameter is null.
	 * 
	 * @param key the property to look up
	 * @return the value corresponding to key, as a long integer.
	 */
	public long getLong( String key ) 
				throws UndefinedParameterException {
		String valueS = null;
		valueS = this.getPropertyAndCheck(key);
		long value = ParseUtils.getLong(valueS);
		// Report to the user if required:
		if( report_parameters ) 
			System.out.println("# Parameter "+key+" has been given the long value "+value);
		return value;
	}

	/**
	 * Gets the parameter with the specified key and return it as a boolean.
	 * Throws an UndefinedParameterException if the parameter is null.
	 * 
	 * @param key the property to look up
	 * @return the value corresponding to key, as a boolean.
	 */
	public boolean getBoolean( String key ) 
				throws UndefinedParameterException {
		String valueS = null;
		valueS = this.getPropertyAndCheck(key);
		boolean value = ParseUtils.getBoolean(valueS);
		// Report to the user if required:
		if( report_parameters ) 
			System.out.println("# Parameter "+key+" has been given the boolean value "+value);
		return value;
	}

	/**
	 * Gets the parameter with the specified key and return it as a String.
	 * Throws an UndefinedParameterException if the parameter is null.
	 * 
	 * @param key the property to look up
	 * @return the value corresponding to key, as a String.
	 */
	public String getString( String key ) 
				throws UndefinedParameterException {
		String value = null;
		value = this.getPropertyAndCheck(key);
		// Report to the user if required:
		if( report_parameters ) 
			System.out.println("# Parameter "+key+" has been given the String value "+value);
		return value;
	}

	/**
	 * Determine whether the user has specified a value for a particular key.
	 * Use to check whether to get optional keys:
	 * @param key the id of the property to test.
	 * @return true if the property has been defined, false otherwise.
	 */
	public boolean isDefined( String key ) {
		// If not defined, return false:
		if( this.getProperty(key) == null ) return false;
		// Otherwise:
		return true;
	}
	
	/**
	 * This class looks up the value and checks that it is set.
	 * @param key to look up
	 * @return the value as a string
	 * @throws UndefinedParameterException if the value for that key is null.
	 */
	private String getPropertyAndCheck( String key ) 
				throws UndefinedParameterException {
		String value = this.getProperty(key);
		if( value == null ) {
			throw new UndefinedParameterException("Parameter "+key+" is not defined!");
		}
		return value;
	}
	
	/**
	 * Override of the superclasses getProperty function.
	 * Does just the same except it trims any whitespace.
	 */
	public String getProperty( String key ) {
		// Call the superclass, but trim any whitespace:
		String value = super.getProperty(key);
		if( value != null ) return value.trim();
		return value;
	}
	
	/**
	 * Inner class defining the Exception where no matching parameter can be found.
	 *
	 * @author ajackso1
	 * 
	 */
	public class UndefinedParameterException extends Exception {
		public UndefinedParameterException( String error ) {
			super(error);
		}
	}
}

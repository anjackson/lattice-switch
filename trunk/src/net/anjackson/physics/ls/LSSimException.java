/**-------------------------------------------------------------
 * jBinLats - LSSimException.java
 * net.anjackson.physics.ls.LSSimException
 * 
 * Created on 17-Feb-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

/**
 * This is the RuntimeException for my simulation class.
 * 
 * Implemented as a RuntimeException as the complexity of declaring the kind of
 * errors that can be detected at various points throughout the code.
 * 
 * Should probably shift to compile-time Exceptions if possible.
 *
 * @author ajackso1
 * @version $Id:$
 *
 */
public class LSSimException extends RuntimeException {
	// The reason code for this exception:
	private int reason = UNKNOWN;

	// Add reason codes and string translations.
	/** If the reason for the exception is unknown: */
	public static final int UNKNOWN = -1;
	/** If the reason for the exception is that the code was requested to abort by some other process. */
	public static final int ABORTED = 0; 
	/** If the reason for the exception is that one of the crystals melted. */
	public static final int MELTED = 1;
	/** If the reason for the exception is that the simulation parameters were not acceptable. */
	public static final int PARAM = 2;
	/** If the reason for the exception is that the simulation failed a self-consistancy check. */
	public static final int CHECK = 3;
	
	/**
	 * Lookup the reason code for this exception.
	 * @return The reason code for the throws exception
	 */
	public int getReason() {
		return reason;
	}

	/**
	 * Look up the textual equivalent of the given reason code:
	 * @param _reason
	 * @return A string describing the reason.
	 */
	static public String getReasonString( int _reason ) {
		if( _reason == ABORTED ) {
			return "Simulation was requested to abort.";
		} else if( _reason == MELTED ) {
			return "One of the crystals melted.";
		} else {
			return "Unknown reason for run-time exception.";
		}
	}
	
	/**
	 * Look up the textual representation of the reason code for this exception:
	 * @return The string that describes the exception in detail.
	 */
	public String getReasonString() {
		return getReasonString(reason);
	}
	
	/**
	 * @param _reason The reason code for the exception, an int flag.
	 */
	public LSSimException( int _reason ) {
		super("LSSimException "+_reason+": "+getReasonString(_reason));
		reason = _reason;
	}
	
	/**
	 * @param _reason The reason code for the exception, an int flag.
	 * @param _cause The exception that caused the problem.
	 */
	public LSSimException( int _reason, Throwable _cause) {
		super(_cause);
		reason = _reason;
	}
	/**
	 * @param _reason The reason code for the exception, an int flag.
	 * @param _message String describing the problem.
	 */
	public LSSimException( int _reason, String _message) {
		super(_message);
		reason = _reason;
	}

	/**
	 * @param _reason The reason code for the exception, an int flag.
	 * @param _message String describing the problem.
	 * @param _cause The exception that caused the problem.
	 */
	public LSSimException( int _reason, String _message, Throwable _cause) {
		super(_message, _cause);
		reason = _reason;
	}


}

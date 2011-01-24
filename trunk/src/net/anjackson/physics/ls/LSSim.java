/**-------------------------------------------------------------
 * jBinLats - LSSim.java
 * net.anjackson.physics.ls.LSSim
 * 
 * Created on 23-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

import net.anjackson.utils.ParamFileReader;

/**
 * This is the main Lattice Switch Simulation class.
 * 
 * It should be run, and passed the relevant arguments from the command line.
 * It will then instanciate a simulation of the right type (via the Factory).
 * Then it will pass the arguments down so that the instance can be initialised.
 * It will then run the simulation.
 *
 *<p>
 * By passing it the -test flag, the code can be made to test itself.  It runs a known-good
 * simulation with a particular seed, and checks the answers are consistent.  This will fail 
 * if changes to code significantly alter the order of operation, or the order of random-number 
 * use.  However, it is a very good way of testing that simple changes like refactoring 
 * have not broken the code.
 *</p>
 *
 *<h3>To Do List</h3>
 *<p>
 *<ul>
 *  <li>IDEA Add a -test mode, see below, which runs tests with fixed seeds and checks for consistent output on some known systems.
 *  Add hs-fcchcp-nvt-n216-d0.7778-wf.dat and other best-known weight functions to the classpath? </li>
 * </ul>
 * @author ajackso1
 * @version $Id: LSSim.java 963 2006-08-14 21:02:10Z anj $
 *
 */
public class LSSim {
	// The simulate to be run:
	static LatticeSwitchSimulation lss = null;	

	/**
	 * Main routine - used to invoke this simulation.
	 * @param args The command line arguments
	 */
	public static final void main(String[] args) {
		// Turn the command-line arguments into a Parameter File instance.
		// Attempt to parse the arguments:
		if( args == null || args.length != 1 ) {
			System.err.println("No parameter file specified!");
			System.err.println("Syntax: [<parameter file>|-test]\n");
			System.exit(0);
		}
		
		// Are we in self-test mode?
		if( "-test".equals(args[0])) {
			boolean testResult = LSSim.testHSfcchcpNVT();
			if( testResult == false ) {
				System.out.print(" Test FAILED! ");
			}
		} else {
			// Take the argument to be the parameter file:
			String filename = args[0];
			
			// Load in the parameters properties file:
			ParamFileReader pars = new ParamFileReader();
			pars.loadParamFile(filename);
			
			// Run the simulation:
			LSSim.runLSSim(pars);
		}
	}

	/**
	 * Runs a LS simulation based on the supplied parameters.
	 * 
	 * @param pars
	 */
	private static void runLSSim( ParamFileReader pars ) {
		// Full debugging info:
		boolean debug = false;
		try {
			if( pars.isDefined("debug") ) debug = pars.getBoolean("debug");
		} catch( ParamFileReader.UndefinedParameterException e ) {}
		
		// Initialise the simulation:
		lss = LSSimFactory.createSimulation(pars);
		lss.init_simulation();

		// Attempt to run:
		try {
			// Run the simulation:
			lss.run_simulation();
		} catch( LSSimException e ) {
			System.out.println("Run-time exception from LatticeSwitchSimulation.run_simulation():\n"+e);
			if( debug ) {
				System.out.println("\n--- Stack Track ---");
				e.printStackTrace(System.out);
			}
		}
		
		// Finish off the simulation - final analysis and output:
		lss.finish_simulation();
		
	}

	/**
	 * IDEA Test code to be invoked here...
	 * 
	 * @return
	 */
	private static boolean testHSfcchcpNVT() {
		// Set up a ParamFileReader with the required parameters, and fire it up.  Check results are consistent.  Use a known good w/f.
		ParamFileReader pars = new ParamFileReader();
		
		// Run the simulation:
		LSSim.runLSSim(pars);
		
		// Test the results against known good values:
		
		return true;
	}
}

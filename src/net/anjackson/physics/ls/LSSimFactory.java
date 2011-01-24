/**-------------------------------------------------------------
 * jBinLats - LSSimFactory.java
 * net.anjackson.physics.ls.LSSimFactory
 * 
 * Created on 20-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

import net.anjackson.physics.ls.switches.AB2Fcc2SwitchMap;
import net.anjackson.physics.ls.switches.AB2FccSwitchMap;
import net.anjackson.physics.ls.switches.AB13Fcc2SwitchMap;
import net.anjackson.physics.ls.switches.AB13FccSwitchMap;
import net.anjackson.physics.ls.switches.AB13HcpSwitchMap;
import net.anjackson.physics.ls.switches.ABCsClFcc2SwitchMap;
import net.anjackson.physics.ls.switches.ABCsClFccSwitchMap;
import net.anjackson.physics.ls.switches.ABCsClHcpSwitchMap;
import net.anjackson.physics.ls.switches.FccHcpSwitchMap;
import net.anjackson.physics.ls.switches.RandomSwitchMap;
import net.anjackson.physics.ls.switches.SwitchMap;
import net.anjackson.utils.ParamFileReader;

/**
 * This class provides a factory that generates instances of the 
 * lattice-switch simulation. The options are, which potential to use, and 
 * which lattice mapping.
 * 
 * The factory instanciates the simulation classes, and puts the right lattices
 * etc into that simulation instance.
 *
 * @author ajackso1
 * @version $Id: LSSimFactory.java 1189 2006-12-14 15:18:43Z anj $
 *
 */
public class LSSimFactory {

	// The known potentials:
	public static final String POT_HS = "hs";
	public static final String POT_LJ = "lj";

	// The known lattice-switch mappings:
	public static final String LAT_FCC_HCP = "fcc-hcp";
	public static final String LAT_AB2_2FCC = "ab2-2fcc";
	public static final String LAT_AB2_FCC = "ab2-fcc";
	public static final String LAT_AB13_2FCC = "ab13-2fcc";
	public static final String LAT_AB13_FCC = "ab13-fcc";
	public static final String LAT_AB13_HCP = "ab13-hcp";
	public static final String LAT_ABCSCL_2FCC = "abCsCl-2fcc";	
	public static final String LAT_ABCSCL_FCC = "abCsCl-fcc";
	public static final String LAT_ABCSCL_HCP = "abCsCl-hcp";
	public static final String LAT_RANDOM = "random";

	// The known ensembles:
	public static final String NVT_SIM = "nvt";
	public static final String NPT_SIM = "npt";
	
	// Debugging log flag:
	public static boolean debug = false;
	
	/**
	 * This is a factory class for instanciating simulations.  Allows
	 * all types of simulation to be generated based on run-time configuration.
	 * 
	 * @param pars - the ParamFileReader containing the simulation properties.
	 * @return A LatticeSwitchSimulation instance, either loaded or new based on the supplied ParamFileReader.
	 */
	public static LatticeSwitchSimulation createSimulation( ParamFileReader pars ) {
		// Determine whether to load from an old file:
		boolean cold_start = true;
		try {
			// Attempt to see if debugging should be on.
			if( pars.isDefined("debug")) debug = pars.getBoolean("debug");
			// Find out if this is a cold start.
			cold_start = pars.getBoolean("cold_start");
		} catch( ParamFileReader.UndefinedParameterException e ) {
			// If any of them is not defined, cry and grumble.
			throw( new LSSimException(LSSimException.PARAM,e));
		}
		
		// Either run from cold start, or reload:
		if( cold_start ) {
			return createNewSimulation(pars);
		} else {
			return loadSimulation(pars);
		}
	}

	/** 
	 * A factory designed to load a serialized simulation from a file:
	 * @param pars the ParamFileReader
	 * @return the loaded simulation.
	 */
	private static LatticeSwitchSimulation loadSimulation( ParamFileReader pars ) {
		// Determine the filename:
		String filename = null;
		try {
			filename = pars.getString("load.config");
		} catch( ParamFileReader.UndefinedParameterException e ) {
			// If any of them is not defined, cry and grumble.
			throw( new LSSimException(LSSimException.PARAM,e));
		}
		// Load it and pass it back.
		LatticeSwitchSimulation lss = LatticeSwitchSimulation.load_config(filename);
		// Pass the debug flag down:
		lss.DEBUG = debug;
		return lss;
	}

	/**
	 * A factory to generate an entirely new simulation.
	 * 
	 * @param pars The ParamFileReader that describes the file.
	 * @return The new LatticeSwitchSimulation instance, created according to the ParamFileReader parameters.
	 */
	private static LatticeSwitchSimulation createNewSimulation( ParamFileReader pars ) {
		LatticeSwitchSimulation sim = null;
		SwitchMap switchmap = null;
		String ensm = null;
		

		// Load up the main simulation parameters:
		try {
			// Determine the ensemble:
			ensm = pars.getString("ensemble");
			
			// Determine the lattice-switch mapping to use:
			String lsmap = pars.getString("switchmap");
			if( LAT_FCC_HCP.equalsIgnoreCase(lsmap) ) {
				switchmap = new FccHcpSwitchMap();
			} else if( LAT_AB2_2FCC.equalsIgnoreCase(lsmap) ) {
				switchmap = new AB2Fcc2SwitchMap();
			} else if( LAT_AB13_2FCC.equalsIgnoreCase(lsmap) ) {
				switchmap = new AB13Fcc2SwitchMap();
			} else if( LAT_ABCSCL_2FCC.equalsIgnoreCase(lsmap) ) {
				switchmap = new ABCsClFcc2SwitchMap();
			} else if( LAT_AB2_FCC.equalsIgnoreCase(lsmap) ) {
				switchmap = new AB2FccSwitchMap();
			} else if( LAT_AB13_FCC.equalsIgnoreCase(lsmap) ) {
				switchmap = new AB13FccSwitchMap();
			} else if( LAT_AB13_HCP.equalsIgnoreCase(lsmap) ) {
				switchmap = new AB13HcpSwitchMap();
			} else if( LAT_ABCSCL_FCC.equalsIgnoreCase(lsmap) ) {
				switchmap = new ABCsClFccSwitchMap();
			} else if( LAT_ABCSCL_HCP.equalsIgnoreCase(lsmap) ) {
				switchmap = new ABCsClHcpSwitchMap();
			} else if( LAT_RANDOM.equalsIgnoreCase(lsmap) ) {
				switchmap = new RandomSwitchMap();
			} else {
				throw( new LSSimException(LSSimException.PARAM,
				"Unknown lattice switch mapping "+lsmap+". "+
				" Should be one of: "+LAT_FCC_HCP+", "+LAT_AB2_2FCC+", "+LAT_AB13_2FCC+", "+LAT_ABCSCL_2FCC+"...[see LSSimFactory code for full list]."));
			}
			
			// Set up the lattices, using the sizes specified in the parameter file:
            if( pars.isDefined("l0.nx")) {
            	// We are specifying both phases seperately:
            	int sizes[][] = new int[2][3];
            	sizes[0][0] = pars.getInteger("l0.nx");
            	sizes[0][1] = pars.getInteger("l0.ny");
            	sizes[0][2] = pars.getInteger("l0.nz");
            	sizes[1][0] = pars.getInteger("l1.nx");
            	sizes[1][1] = pars.getInteger("l1.ny");
            	sizes[1][2] = pars.getInteger("l1.nz");
            	// Initialise the switch map:
            	switchmap.setSizes(sizes);
            	switchmap.init();
            } else {
            	// Determine the system size:
        		int nx = 0, ny = nx, nz = nx;
            	nx = pars.getInteger("n.x");
            	ny = pars.getInteger("n.y");
            	nz = pars.getInteger("n.z");
            	
            	// Initialise the switch map:
            	switchmap.setSize(nx,ny,nz);
            	switchmap.init();
            }
			
			// Determine the potential:
			String pots = pars.getString("potential");
			// Create the simulation:
			if( POT_LJ.equalsIgnoreCase(pots)) {
				sim = new LJLS();
			} else if( POT_HS.equalsIgnoreCase(pots)) {
				// Are we dealing with a three-way binary crystal switch?
				if( LAT_AB2_2FCC.equalsIgnoreCase(lsmap) || 
						LAT_AB13_2FCC.equalsIgnoreCase(lsmap) || 
						LAT_ABCSCL_2FCC.equalsIgnoreCase(lsmap)) {
					sim = new Binary3PhaseHSLS();
				} else if( LAT_AB2_FCC.equalsIgnoreCase(lsmap) || 
						LAT_AB13_FCC.equalsIgnoreCase(lsmap) || 
						LAT_AB13_HCP.equalsIgnoreCase(lsmap) || 
						LAT_ABCSCL_FCC.equalsIgnoreCase(lsmap) ||
						LAT_ABCSCL_HCP.equalsIgnoreCase(lsmap)) {
						sim = new BinaryHSLS();
				} else {
					// This is a straight two-phase switch:
					sim = new HSLS();
				}
			} else {
				throw( new LSSimException(LSSimException.PARAM,
				"Could not parse system potential '"+pots+"'."));
			}
			
		} catch( ParamFileReader.UndefinedParameterException e ) {
			// If any of them is not defined, cry and grumble.
			throw( new LSSimException(LSSimException.PARAM,e));
		}
		
		// ----------------------------------------------------------------
        
		// Pass the SwitchMap to the simulation:
		sim.setSwitchMap(switchmap);

		// Pass the ensemble down:
		if( NPT_SIM.equalsIgnoreCase(ensm)) {
			sim.NPT_SIM = true;
		} else {
			sim.NPT_SIM = false;
		}
		// Now allow the simulation to read the parameters:
		sim.read_params_from_user(pars);

		// Pass the debug flag down:
		sim.DEBUG = debug;

		// Pass the initialised simulation class back:
		return sim;
	}
}

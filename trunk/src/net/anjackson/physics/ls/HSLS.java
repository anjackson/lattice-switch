/**-------------------------------------------------------------
 * jBinLats - HSLS.java
 * net.anjackson.physics.ls.HSLS
 * 
 * Created on 13-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

import net.anjackson.utils.BlockAnalysis;

/**
 * <h2>Introduction</h2>
 * <p>
 * This child of the lattice switch simulation class for monodisperse 
 * hard-sphere systems.  It overrides the necessary parts of the generic base
 * class to produce a valid hard-core simulation.
 * </p><p>
 * It actually simulates with a top-hat hard core of height 1, thus counting 
 * the number of overlaps.  The 'temperature' is set to be very low so that the
 * particles can <i>never</i> overlap.
 * </p>
 * 
 * @author ajackso1
 * @version $Id: HSLS.java 1011 2006-09-04 15:59:50Z anj $
 * 
 */
public class HSLS extends LatticeSwitchSimulation {
		
	/** Block analysis of the volume fraction */
	protected BlockAnalysis[] bVF = { new BlockAnalysis(), new BlockAnalysis() };

	/** 
	  * write_output_header:
	  * - write the definition of the simulation parameters to the 
	  *   standard output stream
	  **/
	protected void write_output_header() {
	    System.out.print(" H Hard-sphere lattice-switch monte carlo code: hcp v fcc:\n");
	    System.out.print(" H \n");
	    // Write out generic system headers
	    super.write_output_header();
	    // Write out HS specific stuff
	    System.out.print(" H -- Hard Spheres --\n");
	    System.out.println(" H dr = "+dr);
	    System.out.println(" H nn cut off = "+nn_cut_off);
	    System.out.println(" H there are ["+max_inn[0]+","+max_inn[1]+"] NNs given a cut-off of "+nn_cut_off);
	    System.out.println(" H The density is calculated to be "+calc_sys_density());
	    System.out.println(" H The volume fraction is calculated to be "+calcVolFrac());
	    System.out.println(" H ");
	}

		
	void setTemperature( double _t ) {
		// Running at a temperature of 0 (or close to) will prevent overlaps.
		super.setTemperature(1.0e-20);// Temperature should be very small, ideally 0.0.  Beware the tyranny of NaN.
	}
	
	void setPressure( double _pressure ) {
		// The factor of beta should be taken out of the pressure, and a factor relating to the overall scaling (sphere size) should be inserted.
		super.setPressure(_pressure);
		BetaPres = Pres;
	}


	protected void init_simulation() {
		// Initialize the simulation as normal.
		super.init_simulation();
		// As this is the hard-sphere case, any initial condition where e!=0 is WRONG!
		if( Math.abs(e) > E_TOL ) {
			// If this check fails, warn and die:
			throw( new LSSimException(LSSimException.PARAM,
					"ERROR: Hard-sphere system has overlaps in it's initial state: "+latt_name[c_lat]+" e = "+e));
		}
		if( !NPT_SIM && !zero_phase_ig && Math.abs(em) > E_TOL ) {
			// If this check fails, warn and die:
			throw( new LSSimException(LSSimException.PARAM,
					"ERROR: Hard-sphere system has overlaps in it's initial conjugate state: "+latt_name[m_lat]+" em = "+em));
		}
	}
	
	protected void report_simulation(int sweep) {
		// Do the superclass reporting:
		super.report_simulation(sweep);
		
		// For NPT, also output the volume fraction:
		if( NPT_SIM ) {
			System.out.println(" Vf "+calcVolFrac());
		}
	}
	
	protected void check_simulation(int sweep) {
		// Do all the usual checks:
		super.check_simulation(sweep);
		
		// TODO Should seperate checking and check-frequency measurements, so that equilibration times are always ignored.
		
		// Do additional volume-fraction analysis in NPT:
		if( NPT_SIM && sweep > this.equib_sweeps ) {
			// Add the current volume fraction to the block analyser:
			bVF[c_lat].add(sweep, calcVolFrac());
			// Output the current block analysis results:
			for( int ila=0; ila < 2; ila++ ) {
				if( bVF[ila].hasData() ) {
					System.out.println(" CV lat"+ila+" vf "+bVF[ila].getMean()+" +/- "+bVF[ila].getStdErr());
				}
			}
		}
		
	}
	
	/** 
	  * ij_inter:
	  * - calculates the pair potential as a function of dr^2.
	  * - a hard-core of height 1.
	  **/
	protected double ij_inter(double dr2, int i, int j, int il ) {
	    
	    if ( dr2 < diameter2 ) {
	    	return(1.0);
	    } else {
	    	return(0.0);
	    }

	}
	
	/**
	 * The actual pair-potential interaction in polydisperse systems:
	 * Calculates the pair potential as a function of dr^2 and pdiameter[il][n].
	 * 
	 * @param dr2 the seperation distance squared.
	 * @param i the identity of one particle
	 * @param j the identity of the other particle (usually i>j)
	 * @param il the lattice index - only required when performing potential-switches etc.
	 * @return The potential for the given distance^2
	 */
	protected double ij_inter_poly(double dr2, int i, int j, int il) {
		// Calculate the overlap seperation for these two particles:
		double sep = 0.5*(pdiameter[il][i]+pdiameter[il][j]);
		// If the given dr2 is less than this minimum seperation, return a 1.0;
		if ( dr2 < sep*sep ) {
	    	return(1.0);
	    } else {
	    	return(0.0);
	    }
	}
	

	/**
	 * Calculate the gound-state energy of a lattice.
	 * Zero for all hard-sphere systems.
	 */
	protected double calc_gs_by_struct(double densi, int latti) {
		// The Ground-state energy of the HS system is 0
		return 0;
	}



	/**
	 * Calculate the number of overlaps generated by this move:
	 */
	public double calc_dE_sphere_move(int ip, int t_lat) {
		// If we are moving in the current phase, we can assume the current local overlap count is zero:
		if( t_lat == c_lat && !zero_phase_ig) {
			return calc_local_energy(ip,t_lat,TRY_POS);
		} else {
			return calc_local_energy(ip,t_lat,TRY_POS) - calc_local_energy(ip,t_lat,CUR_POS);
		}
	}

	/**
	  * mc_basis_flip:
	  * - flip to other basis, if we can:
	  * Overriding the general formulation.
	  * This is just done as it is quicker than the more general routine.
	  **/
	protected void mc_basis_flip() {
		// Flip if m is zero, but note that m is now double, so the running average may wander:
	    if ( Math.abs(em) < E_TOL ) {
	    	c_lat = m_lat;
	    	m_lat = 1 - c_lat;
	    	// No need to swap energies as both are zero.
	    }
	}
	
	/**
	 * calculate the volume fraction for monodisperse hard spheres
	 */
	protected double calcVolFrac() {
		return n*sphere_volume(diameter)/(box[c_lat].x*box[c_lat].y*box[c_lat].z);
	}

}

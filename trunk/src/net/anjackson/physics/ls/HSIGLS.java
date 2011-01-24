/**-------------------------------------------------------------
 * LatticeSwitch - HSIGLS.java
 * net.anjackson.physics.ls.HSIGLS
 * 
 * Created on 25 Aug 2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

/**
 * Just like HSLS, but the potential is switch to an ideal-gas model
 * (i.e. no potential) in the LATT=0 phase.
 *
 * @author ajackso1
 * @version $Id$
 *
 */
public class HSIGLS extends HSLS {

	/** 
	  * ij_inter:
	  * - calculates the pair potential as a function of dr^2.
	  * - a hard-core of height 1.
	  **/
	protected double ij_inter(double dr2, int i, int j, int il ) {
		// Ideal gas:
	    if( il == 0 ) return 0.0;
	    // Normal:
	    return super.ij_inter(dr2, i, j, il);
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
		// Ideal gas:
	    if( il == 0 ) return 0.0;
	    // Normal:
	    return super.ij_inter_poly(dr2, i, j, il);
	}

	
}

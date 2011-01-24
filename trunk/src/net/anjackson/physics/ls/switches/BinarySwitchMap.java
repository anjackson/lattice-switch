/**-------------------------------------------------------------
 * jBinLats - ThreePhaseSwitchMap.java
 * net.anjackson.physics.ls.switches.ThreePhaseSwitchMap
 * 
 * Created on 26-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

/**
 * This is another switchmap but for the generalized case where a single
 * AB phase is being switched into a single A phase.  The A particles are listed
 * before the B particles, by convention.  
 * 
 * @author ajackso1
 * @version $Id$
 *
 */
public abstract class BinarySwitchMap extends SwitchMap {
	// The number of 'B' particles:
	int Ntp = -1;

	/**
	 * Where in the lattice array do to A particles end and the B's begin?
	 * i.e. How many A particles are there?
	 * 
	 * @return Index of the first particle of the B's, in the latt array.
	 */
	public int phaseIndex() {
		return n-Ntp;
	}
	
	/**
	 * Find out how many 'B' particles are there?
	 * @return The integer number of 'B' particles.
	 */
	public int getThirdN() {
		return Ntp;
	}
}

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

import net.anjackson.maths.VectorD3;
import 	net.anjackson.physics.lattices.Lattice;

/**
 * This is another switchmap but for the generalized case where a single
 * phase is being switched into two seperate phases.  The volumes and identities
 * of the first particles are by convention common to the first two phases.  The
 * 'extra' particles in the first phase are moved into the third-phase, which has
 * it's own volume and initial density.
 * 
 * i.e. the first two phases are treated just like the two phases of a lattice
 * switch, and the information required for the third phase is stored here.
 * 
 * @author ajackso1
 * @version $Id$
 *
 */
public abstract class ThreePhaseSwitchMap extends BinarySwitchMap {
	// The size of the initial box for the third-phase: 
	VectorD3 boxT = new VectorD3();
	// The lattice for the 'Third Phase'
	Lattice latT = null;
	
	/**
	 * This method gets the initial size of the third-phase box, for the overall
	 * density.
	 * @return A VectorD3 containing the box size in the x,y & z directions.
	 */
	public VectorD3 getThirdBox() {
		System.out.println(" DEBUG TPSM.getThirdBox() returning ["+boxT.x+","+boxT.y+","+boxT.z+"].");
		return boxT;
	}

}

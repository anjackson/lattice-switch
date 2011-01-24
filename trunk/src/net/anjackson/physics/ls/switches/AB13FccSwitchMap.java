/**-------------------------------------------------------------
 * jBinLats - AB2Fcc13SwitchMap.java
 * net.anjackson.physics.ls.switches.AB13Fcc2SwitchMap
 * 
 * Created on 25-Jun-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.lattices.AB13Lattice;
import net.anjackson.physics.lattices.FaceCentredCubicABCLattice;
import net.anjackson.physics.lattices.Lattice;

/**
 * This is the AB13 mapping where the AB13 crystal is mapped onto
 * two seperate phases, one of A particles (fcc) and one of 13n B particles (hcp).
 * 
 * 
 * @author ajackso1
 * @version $Id: AB2Fcc2SwitchMap.java 751 2006-05-19 12:08:15Z anj $
 *
 */
public class AB13FccSwitchMap extends BinarySwitchMap {

	/**
	 * Default constructor, must define npc
	 */
	public AB13FccSwitchMap() {
		npc[0] = 14;
		npc[1] = 1;
	}

	/**
	 * Constructor to create the three lattices for this binary switch.
	 * 
	 */
	public void generate() {
		
		// For this system, there are two particles per cell:
		int nc = size[0][0]*size[0][1]*size[0][2];
		n = 14*nc;
		
		// The AB2 Lattice:
		AB13Lattice ab13 = new AB13Lattice();
		ab13.setSizeInCells(size[0][0],size[0][1],size[0][2]);
		
		// Stick them in the latt member:
		latt = new VectorD3[2][];
		latt[0] = new VectorD3[n];		
		System.arraycopy(ab13.getLattice(0),0,latt[0], 0, nc);
		System.arraycopy(ab13.getLattice(1),0,latt[0], nc, 13*nc);
		latt_name[0] = "AB13";
		ab13.saveAsWRL("initAB13");
		
		// Build the A and B lattices, A first, then B second:

		// The seperate A and B2 lattices
		FaceCentredCubicABCLattice fcc = new FaceCentredCubicABCLattice();
		
		// Create the new system at the new size:
		Ntp = n-nc;
		fcc.setSizeInCells(size[1][0],size[1][1],size[1][2]);

		// Copy the two FCC lattices into the 1th phase:
		latt[1] = new VectorD3[n];
		latt_name[1] = "FCC";
		System.arraycopy(fcc.getLattice(), 0, latt[1], 0, n);

		/* Define distances between neighbouring atoms in the same stacking plane 
		 *  taking advantage of the fact that both lattices have the same Unit Cell: */
		double[] c_uc = new double[3];
		c_uc = ab13.getUnitCell();
		uc.x = c_uc[0];
		uc.y = c_uc[1];
		uc.z = c_uc[2];
		
		// Put the lattices in the class members, for later use:
		lats[0] = (Lattice)ab13;
		lats[1] = (Lattice)fcc;


	}

}

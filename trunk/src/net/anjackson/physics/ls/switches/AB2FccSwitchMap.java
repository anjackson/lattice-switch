/**-------------------------------------------------------------
 * jBinLats - AB2Fcc2SwitchMap.java
 * net.anjackson.physics.ls.switches.AB2Fcc2SwitchMap
 * 
 * Created on 20-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.lattices.AB2Lattice;
import net.anjackson.physics.lattices.FaceCentredCubicABCLattice;
import net.anjackson.physics.lattices.Lattice;

/**
 * This is the AB2 mapping where the AB2 crystal is mapped onto
 * two seperate phases, one of A particles (fcc) and one of 2n B particles (hcp).
 * 
 * @author ajackso1
 * @version $Id: AB2Fcc2SwitchMap.java 953 2006-08-08 16:38:23Z anj $
 *
 */
public class AB2FccSwitchMap extends BinarySwitchMap {

	/**
	 * Default constructor, must define npc
	 */
	public AB2FccSwitchMap() {
		npc[0] = 3;
		npc[1] = 1;
	}

	/**
	 * Constructor to create the two lattices for this switch.
	 * 
	 */
	public void generate() {
		
		// For this system, there are two particles per cell:
		int nc = size[0][0]*size[0][1]*size[0][2];
		n = 3*nc;
		
		// The AB2 Lattice:
		AB2Lattice ab2 = new AB2Lattice();
		ab2.setSizeInCells(size[0][0],size[0][1],size[0][2]);
		
		// Stick them in the latt member:
		latt = new VectorD3[2][];
		latt[0] = new VectorD3[n];		
		System.arraycopy(ab2.getLattice(0),0,latt[0], 0, nc);
		System.arraycopy(ab2.getLattice(1),0,latt[0], nc, 2*nc);
		latt_name[0] = "AB2";
		ab2.saveAsWRL("initAB2");
		
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
		c_uc = ab2.getUnitCell();
		uc.x = c_uc[0];
		uc.y = c_uc[1];
		uc.z = c_uc[2];
		
		// Put the lattices in the class members, for later use:
		lats[0] = (Lattice)ab2;
		lats[1] = (Lattice)fcc;
		
	}

}

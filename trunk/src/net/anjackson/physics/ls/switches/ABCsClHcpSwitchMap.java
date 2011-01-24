/**-------------------------------------------------------------
 * jBinLats - ABCsClFcc2SwitchMap.java
 * net.anjackson.physics.ls.switches.ABCsClFcc2SwitchMap
 * 
 * Created on 22-Mar-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.lattices.ABCsClLattice;
import net.anjackson.physics.lattices.HexagonalClosePackedLattice;
import net.anjackson.physics.lattices.Lattice;

/**
 * This class creates an AB CsCl to 2 * FCC lattice mapping, with
 * two FCC crystals each of nx*ny*nz particles being mapped onto a single 
 * CsCl structure.
 * 
 * @author ajackso1
 * @version $Id$
 *
 */
public class ABCsClHcpSwitchMap extends BinarySwitchMap {
	static final long serialVersionUID = 1L;

	/**
	 * Default constructor, must define npc
	 */
	public ABCsClHcpSwitchMap() {
		npc[0] = 2;
		npc[1] = 1;
	}
	
	/**
	 * 
	 */
	public void generate() {
		
		// For this system, there are two particles per cell:
		int nc = size[0][0]*size[0][1]*size[0][2];
		n = 2*nc;
		
		// The AB2 Lattice:
		ABCsClLattice ab = new ABCsClLattice();
		ab.setSizeInCells(size[0][0],size[0][1],size[0][2]);
		
		// Stick them in the latt member:
		latt = new VectorD3[2][];
		latt[0] = new VectorD3[n];		
		System.arraycopy(ab.getLattice(0),0,latt[0], 0, nc);
		System.arraycopy(ab.getLattice(1),0,latt[0], nc, nc);
		latt_name[0] = "ABCsCl";
		ab.saveAsWRL("initABCsCl");
		
		// Build the A and B lattices, A first, then B second:

		// The seperate A and B2 lattices
		HexagonalClosePackedLattice hcp = new HexagonalClosePackedLattice();
		
		// Create the new system at the new size:
		Ntp = n-nc;
		hcp.setSizeInCells(size[1][0],size[1][1],size[1][2]);

		// Copy the two FCC lattices into the 1th phase:
		latt[1] = new VectorD3[n];
		latt_name[1] = "HCP";
		System.arraycopy(hcp.getLattice(), 0, latt[1], 0, n);

		/* Define distances between neighbouring atoms in the same stacking plane 
		 *  taking advantage of the fact that both lattices have the same Unit Cell: */
		double[] c_uc = new double[3];
		c_uc = hcp.getUnitCell();
		System.out.println(" H DEBUG FCC UC [ "+c_uc[0]+" "+c_uc[1]+" "+c_uc[2]+" ]");
		c_uc = ab.getUnitCell();
		System.out.println(" H DEBUG ABCsCl UC [ "+c_uc[0]+" "+c_uc[1]+" "+c_uc[2]+" ]");
		uc.x = c_uc[0];
		uc.y = c_uc[1];
		uc.z = c_uc[2];
		
		// Put the lattices in the class members, for later use:
		lats[0] = (Lattice)ab;
		lats[1] = (Lattice)hcp;
	}

}

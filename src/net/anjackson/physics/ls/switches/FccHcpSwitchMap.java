/**-------------------------------------------------------------
 * jBinLats - FccHcpSwitchMap.java
 * net.anjackson.physics.ls.switches.FccHcpSwitchMap
 * 
 * Created on 18-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.physics.lattices.FaceCentredCubicABCLattice;
import net.anjackson.physics.lattices.HexagonalClosePackedLattice;

/**
 * This class defines the lattice-switch mapping corresponding to the 
 * original FCC-HCP slide switch.
 *
 * @author ajackso1
 * @version $Id: FccHcpSwitchMap.java 982 2006-08-21 14:51:25Z anj $
 *
 */
public class FccHcpSwitchMap extends SwitchMap {

	/**
	 * Default constructor, must define npc
	 */
	public FccHcpSwitchMap() {
		npc[0] = 1;
		npc[1] = 1;
	}
	

	/**
	 * Constructor to initialise this mapping, creating the 
	 * required lattices and defining the mapping between them.
	 * 
	 */
	public void generate() {
		
		int nx = size[0][0];
		int ny = size[0][0];
		int nz = size[0][0];
		int nc = nx*ny*nz;
		n = nc;
		
		// The HCP Lattice:
		latt_name[0] = "hcp";
		lats[0] = new HexagonalClosePackedLattice();
		
		// The FCC lattice
		lats[1] = new FaceCentredCubicABCLattice();
		latt_name[1] = "fcc";

		// Generate them and stick them in the latt member:
		for( int il = 0; il < 2; il++) {
			lats[il].setSizeInCells(nx,ny,nz);
			latt[il] = lats[il].getLattice();
		}
		
		/* Define distances between neighbouring atoms in the same stacking plane 
		 *  taking advantage of the fact that both lattices have the same Unit Cell: */
		double[] fcc_uc = new double[3];
		fcc_uc = lats[1].getUnitCell();
		uc.x = fcc_uc[0];
		uc.y = fcc_uc[1];
		uc.z = fcc_uc[2];
		
	}	
}

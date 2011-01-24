/**-------------------------------------------------------------
 * jBinLats - RandomSwitchMap.java
 * net.anjackson.physics.ls.switches.RandomSwitchMap
 * 
 * Created on 25-Aug-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.physics.lattices.RandomLattice;
import net.anjackson.maths.VectorD3;

/**
 * This class defines a random 'lattice' switchmap, although in this
 * case the same random lattice is used in both phases.
 *
 * @author ajackso1
 * @version $Id: FccHcpSwitchMap.java 982 2006-08-21 14:51:25Z anj $
 *
 */
public class RandomSwitchMap extends SwitchMap {

	/**
	 * Default constructor, must define npc
	 */
	public RandomSwitchMap() {
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
		
		// The Random Lattice:
		latt_name[0] = "Random";
		lats[0] = new RandomLattice();
		latt_name[1] = "Random";
		lats[1] = lats[0];

		// Generate them and stick them in the latt member:
		lats[0].setSizeInCells(nx,ny,nz);
		latt[0] = lats[0].getLattice();
		// Make latt[1] a clone of latt[0]:
		latt[1] = new VectorD3[n];
		for( int i=0; i<n; i++ ) {
			latt[1][i] = new VectorD3();
			latt[1][i].x = latt[0][i].x;
			latt[1][i].y = latt[0][i].y;
			latt[1][i].z = latt[0][i].z;
		}
		
		/* Taking advantage of the fact that both lattices have the same Unit Cell: */
		double luc[] = lats[0].getBasicUnitCell();
		uc.x = luc[0];
		uc.y = luc[1];
		uc.z = luc[2];
		
	}	
}

/**-------------------------------------------------------------
 * jBinLats - ABNaClLattice.java
 * net.anjackson.physics.lattices.ABNaClLattice
 * 
 * Created on 14-Feb-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.lattices;

import net.anjackson.maths.VectorD3;

/**
 * This class create a simple FCC lattice.
 * 
 * Here, we use two particles per cell, with two interlaced
 * simple cubic lattices being used to create FCC.
 * 
 * @author ajackso1
 * @version $Id: BCCLattice.java 993 2006-08-25 11:07:06Z anj $
 *
 */
public class BCCLattice extends Lattice {
	
	/**
	 * Constructor:
	 *   Initialised the presence of two lattice arrays:
	 *
	 */
	public BCCLattice() {
		// The unit cell, same as FCC in x and y, but out to 1.0 in Z:
		uc[0] = 1.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		//uc_sep = 1.0;
		// Three sites per unit cell:
		nspc = 2;
		name = "BCC";
	}

	/** 
	 * Generates a BCC Lattice.
	 * 
	 * (non-Javadoc)
	 * @see net.anjackson.physics.lattices.Lattice#generate()
	 */
	protected void generate() {
		// Exit if the system size has not been set.
    	if( !isInitialised() ) return;
    	
		// The number of atoms on each lattice:
		int nat = nx*ny*nz;
		// Total number of sites:
		na = nspc*nat;
		
		// There will be two lattices:
		lat = new VectorD3[na];
		
		// Instanciate an FCC lattice generator:
		SimpleCubicLattice sc = new SimpleCubicLattice();
		
		// Generate the cube lattice:
		sc.setSizeInCells(nx,ny,nz);
		sc.generate();
		System.arraycopy(sc.getLattice(), 0, lat, 0, nat);
		
		// Generate the face-centered lattice:
		sc.setOriginOffset( ucx*uc[0]*0.5, ucy*uc[1]*0.5, ucz*uc[2]*0.5 );
		sc.generate();
		System.arraycopy(sc.getLattice(), 0, lat, nat, nat);		
	}
	
	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public boolean checkSize(int nx, int ny, int nz) {
    	return true;
    }

}

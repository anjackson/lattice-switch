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
 * This class creates AB NaCl lattices.
 * 
 * An AB2 lattice can be built from three FCC lattices.  Relative to the A
 * lattice, the other two lattices (for the smaller particles) are placed at
 * the following offsets (shifted origins, in units of the unit cell):
 *  - A is at  (  x=0,  y=0,  z=0 )
 *  - B1 is at (   +0, +1/2, +1/2 )
 *  - B2 is at ( +1/2, +1/4, +1/2 )
 * 
 * @author ajackso1
 * @version $Id: ABNaClLattice.java 993 2006-08-25 11:07:06Z anj $
 *
 */
public class ABNaClLattice extends MultipleLattice {
	
	/**
	 * Constructor:
	 *   Initialised the presence of two lattice arrays:
	 *
	 */
	public ABNaClLattice() {
		// The unit cell, same as FCC in x and y, but out to 1.0 in Z:
		uc[0] = Math.sqrt(3.0)/2.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		uc[2] = Math.sqrt(2.0/3.0);
		//uc_sep = 1.0;
		// Three sites per unit cell:
		nspc = 2;
		name = "ABNaCl";
	}

	/** 
	 * Generates an AB2 Lattice.
	 * 
	 * (non-Javadoc)
	 * @see net.anjackson.physics.lattices.Lattice#generate()
	 */
	protected void generate() {
		// Exit if the system size has not been set.
    	if( !isInitialised() ) return;

		// There will be two lattices:
		lats = new VectorD3[2][];
		
		// The number of atoms on each lattice:
		int nat = nx*ny*nz;
		// Total number of sites:
		na = nspc*nat;
		
		// Instanciate an FCC lattice generator:
		FaceCentredCubicABCLattice fcc = new FaceCentredCubicABCLattice();
		
		// Generate the A lattice:
		fcc.setSizeInCells(nx,ny,nz);
		fcc.generate();
		lats[0] = fcc.getLattice();
		
		// Generate the B lattice:
		fcc.setOriginOffset( ucx*uc[0]*(1.0/3.0), ucy*uc[1]*0.5, ucz*uc[2]*0.5 );
		fcc.generate();
		lats[1]= fcc.getLattice();
		
	}
	
	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public boolean checkSize(int nx, int ny, int nz) {
    	if( nx%2 != 0 ) return false;
    	if( nz%3 != 0 ) return false;
    	return true;
    }

}

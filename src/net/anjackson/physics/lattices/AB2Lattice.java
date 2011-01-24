/**
 * 
 */
package net.anjackson.physics.lattices;

import net.anjackson.maths.*;

/**
 * This class creates AB2 lattices.
 * 
 * An AB2 lattice can be built from three AAA lattices.  Relative to the A
 * lattice, the other two lattices (for the smaller particles) are placed at
 * the following offsets (shifted origins, in units of the unit cell):
 *  - A is at  (  x=0,  y=0,  z=0 )
 *  - B1 is at (   +0, +1/2, +1/2 )
 *  - B2 is at ( +1/2, +1/4, +1/2 )
 * 
 * @author ajackso1
 * @version $Id: AB2Lattice.java 968 2006-08-15 16:54:11Z anj $
 *
 */
public class AB2Lattice extends MultipleLattice {
	
	/**
	 * Constructor:
	 *   Initialised the presence of two lattice arrays:
	 *
	 */
	public AB2Lattice() {
		// The unit cell, same as FCC in x and y, but out to 1.0 in Z:
		uc[0] = Math.sqrt(3.0)/2.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		//uc[2] = Math.sqrt(2.0/3.0);
		uc_sep = 1.0;
		// Three sites per unit cell:
		nspc = 3;
		name = "AB2";
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
		AAALattice aaa = new AAALattice();
		
		// Generate the A lattice:
		aaa.setSizeInCells(nx,ny,nz);
		aaa.generate();
		lats[0] = aaa.getLattice();
		
		// Generate the B1 lattice:
		aaa.setOriginOffset( ucx*uc[0]*(1.0/3.0), ucy*uc[1]*0.5, ucz*uc[2]*0.5 );
		aaa.generate();
		VectorD3[] B1 = aaa.getLattice();

		// Generate the B2 lattice:
		aaa.setOriginOffset( ucx*uc[0]*(2.0/3.0), ucy*uc[1]*0.0, ucz*uc[2]*0.5 );
		aaa.generate();
		VectorD3[] B2 = aaa.getLattice();
		
		// Copy the fcc B1 and B2 lattices into one array (lats[1]):
		lats[1] = new VectorD3[2*nat];
		System.arraycopy(B1, 0, lats[1], 0, nat);
		System.arraycopy(B2, 0, lats[1], nat, nat);
	}
	
	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public boolean checkSize(int nx, int ny, int nz) {
    	if( nx%2 != 0 ) return false;
    	return true;
    }

}

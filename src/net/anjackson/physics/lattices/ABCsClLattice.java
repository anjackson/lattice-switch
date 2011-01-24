/**-------------------------------------------------------------
 * jBinLats - ABCsClLattice.java
 * net.anjackson.physics.lattices.ABCsClLattice
 * 
 * Created on 21-Mar-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.lattices;

import net.anjackson.maths.VectorD3;

/**
 * Creates the AB CsCl lattice by combining two SC lattices.
 * The minor sites lie at (0.5,0.5,0.5) of the SC unit cell.
 *
 * @author ajackso1
 * @version $Id: ABCsClLattice.java 993 2006-08-25 11:07:06Z anj $
 *
 */
public class ABCsClLattice extends MultipleLattice {

	/**
	 * Set up the CsCl unit cell.
	 */
	public ABCsClLattice() {
		uc[0] = 1.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		uc_sep = 1.0;
		nspc = 2;
		name = "ABCsCl";
	}

	/** Generate a CsCl lattice by combining to SC lattices.
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
		SimpleCubicLattice sc = new SimpleCubicLattice();
		
		// Generate the A lattice:
		sc.setSizeInCells(nx,ny,nz);
		sc.generate();
		lats[0] = sc.getLattice();
		
		// Generate the B lattice:
		sc.setOriginOffset( ucx*uc[0]*0.5, ucy*uc[1]*0.5, ucz*uc[2]*0.5 );
		sc.generate();
		lats[1]= sc.getLattice();	
	}

}

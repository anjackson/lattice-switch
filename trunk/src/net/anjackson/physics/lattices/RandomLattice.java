/**-------------------------------------------------------------
 * jBinLats - RandomLattice.java
 * net.anjackson.physics.lattices.RandomLattice
 * 
 * Created on 25-Aug-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.lattices;

import edu.cornell.lassp.houle.RngPack.RandomSeedable;
import edu.cornell.lassp.houle.RngPack.Ranmar;
import net.anjackson.maths.VectorD3;

/**
 * This class create a simple Random 'lattice'.
 * 
 * @author ajackso1
 * @version $Id: RandomLattice.java 996 2006-08-25 13:03:14Z anj $
 *
 */
public class RandomLattice extends Lattice {
	
	/**
	 * Constructor:
	 *   With unit cell so that radius of 1 with one particle per cell gives
	 *   the same density as a close-packed crystal.
	 */
	public RandomLattice() {
		uc[0] = Math.pow(2.0,0.5)*Math.pow(0.25,1.0/3.0);
		uc[1] = uc[0];
		uc[2] = uc[0];
		nspc = 1;
		name = "Random";
	}

	/** 
	 * Generates a random 'lattice'.
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
		
		// Determine overall box size:
		double rbox[] = getSystemSize();
		
		// Initialise a random number generator:
		long iseed = RandomSeedable.ClockSeed();
		RandomSeedable rng = new Ranmar(iseed);
		
		// Place site at random locations in this box:
		for( int i=0; i< na; i++ ) {
			lat[i] = new VectorD3();
			lat[i].x = rng.raw()*rbox[0];
			if( i == 0 ) System.out.println("Generated "+i+" "+lat[i].x+" "+rbox[0]);
			lat[i].y = rng.raw()*rbox[1];
			lat[i].z = rng.raw()*rbox[2];
		}
		
	}
	
	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public boolean checkSize(int nx, int ny, int nz) {
    	return true;
    }

}

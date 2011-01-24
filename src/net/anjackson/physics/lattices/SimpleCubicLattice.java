/**-------------------------------------------------------------
 * jBinLats - SimpleCubicLattice.java
 * net.anjackson.physics.lattices.SimpleCubicLattice
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
 * A simple AAA packing lattice.  The x/y places contain normal hexagonally close-packed
 * layers, but the layers are placed directly on top of each other (in the A position).
 * Vertical spacing is therefore 1.
 *
 * @author ajackso1
 * @version $Id$
 *
 */
public class SimpleCubicLattice extends Lattice {

	/**
	 * The constructor, used to define the unit-seperation unit cell.
	 */
	public SimpleCubicLattice() {
		uc[0] = 1.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		uc_sep = 1.0;
	}

	/** Generate a SC lattice:
	 * @see net.anjackson.physics.lattices.Lattice#generate()
	 */
	protected void generate() {
		int ix,iy,iz,in;
		double xs,ys,zs;
		na = nx*ny*nz;
	
		/* Allocate memory for the lattice */
		lat = new VectorD3[na];
		
		/* Define distances between neighbouring atoms in the same stacking plane, 
		 * rescaled using the desired unit cell: */
		xs = uc[0] * ucx;
		ys = uc[1] * ucy;
		zs = uc[2] * ucz;
		
		/* Initialise particle counter */
		in = 0;
		/* Loop over the x-y planes */
		for ( iz = 1; iz <= nz; iz++ ) {
			/* Loop over x-dirn */
			for ( ix = 1; ix <= nx; ix++ ) {
				/* Loop over y-dirn */
				for ( iy = 1; iy <= ny; iy++ ) {
					lat[in] = new VectorD3();
					lat[in].x=loox + xs*((double)ix-1);
					lat[in].y=looy + ys*((double)iy-1);
					lat[in].z=looz + zs*((double)iz-1);
					in++;
				}
			}
		}
	}

}

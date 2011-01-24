/**-------------------------------------------------------------
 * jBinLats - AAALattice.java
 * net.anjackson.physics.lattices.AAALattice
 * 
 * Created on 02-Feb-2006 by ajackso1.
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
 * @version $Id: AAALattice.java 993 2006-08-25 11:07:06Z anj $
 *
 */
public class AAALattice extends Lattice {

	// The constructor, used to define the unit-seperation unit cell.
	public AAALattice() {
		uc[0] = Math.sqrt(3.0)/2.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		uc_sep = 1.0;
		name = "AAA";
	}
	
	protected void generate() {
		int ix,iy,iz,in;
		double xs,ys,zs,t;
		double xo,yo,zo,xoo,yoo;
		na = nx*ny*nz;
	
		/* Allocate memory for the lattice */
		lat = new VectorD3[na];
		
		
		/* Define offset in the x-dirn between A and B planes */
		t  = Math.sqrt(3.0)/3.0;
		/* Define distances between neighbouring atoms in the same stacking plane, 
		 * rescaled using the desired unit cell: */
		t  = t * ucx;
		xs = uc[0] * ucx;
		ys = uc[1] * ucy;
		zs = uc[2] * ucz;
		
		/* Initialise particle counter */
		in = 0;
		/* Loop over the x-y planes */
		for ( iz = 1; iz <= nz; iz++ ) {
			xo = 0.0;
			yo = 0.0;
			zo = zs * ( (double) (iz-1) );
			/* Loop over x-dirn */
			for ( ix = 1; ix <= nx; ix++ ) {
				xoo = xo + xs*( (double) (ix-1) );
				yoo = yo + ys*( (double) ((ix-1)%2) )/2.0;
				/* Loop over y-dirn */
				for ( iy = 1; iy <= ny; iy++ ) {
					lat[in] = new VectorD3();
					lat[in].x=loox + xoo;
					lat[in].y=looy + yoo+ys*((double)(iy-1));
					lat[in].z=looz + zo;
					in++;
				}
			}
		}	    		
	}

}

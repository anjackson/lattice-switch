package net.anjackson.physics.lattices;

import net.anjackson.maths.VectorD3;

/**
 * This is the basic CP lattice generator. 
 * 
 * For FCC-HCP switch, should be multiples of: nx = 6, ny = 6, nz = 6;
 *
 * @author ajackso1
 * @version $Id: ClosePackedLattice.java 516 2006-01-25 18:17:55Z anj $
 *
 */
public abstract class ClosePackedLattice extends Lattice {
    // Flags for which lattice to generate:
	public static final int FCC = 0;
	public static final int HCP = 1;
	
	// The constructor, used to define the unit-seperation unit cell.
	public ClosePackedLattice() {
		uc[0] = Math.sqrt(3.0)/2.0;
		uc[1] = 1.0;
		uc[2] = Math.sqrt(2.0/3.0);
		uc_sep = 1.0;
	}
	
	protected void cpGenerator( int nx, int ny, int nz, int latType ) {
		int ix,iy,iz,in;
		double xs,ys,zs,t;
		double xo,yo,zo,xoo,yoo;
		na = nx*ny*nz;
	
		/* Allocate memory for the lattices and displacements */
		double shift_x = 0.0;
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
			/* Define origin of the x-y planes */
			if ( iz%2 == 1 ) {
				xo = 0.0;
			} else {
				// ie if ( iz%2 == 0 ) 
				xo = +t;
			}
			if ( iz%6 == 1 || iz%6 == 2 ) shift_x = 0.0;
			if ( iz%6 == 3 || iz%6 == 4 ) shift_x = 2.0*t;
			if ( iz%6 == 5 || iz%6 == 0 ) shift_x = 1.0*t;
			yo = 0.0;
			zo = zs * ( (double) (iz-1) );
			/* Loop over x-dirn */
			for ( ix = 1; ix <= nx; ix++ ) {
				xoo = xo + xs*( (double) (ix-1) );
				yoo = yo + ys*( (double) ((ix-1)%2) )/2.0;
				/* Loop over y-dirn */
				for ( iy = 1; iy <= ny; iy++ ) {
					lat[in] = new VectorD3();
					if( latType == HCP ) {
						lat[in].x=loox + xoo;
						lat[in].y=looy + yoo+ys*((double)(iy-1));
						lat[in].z=looz + zo;
					} else if( latType == FCC ) {
						lat[in].x=loox + xoo+shift_x;
						lat[in].y=looy + yoo+ys*((double)(iy-1));
						lat[in].z=looz + zo;		    	
					} else {
						System.out.println("ClosePackedLattice.cpGenerator: Type of lattice unrecognized!  Placing site "+in+" at [0,0,0].");
						// Set the positions to zero:
						lat[in].x = 0;
						lat[in].y = 0;
						lat[in].z = 0;
					}
					in++;
				}
			}
		}	    		
	}

}

/**
 * 
 */
package net.anjackson.physics.lattices;

import net.anjackson.maths.*;

/**
 * This class creates AB13 lattices.
 * 
 * This lattice is composed of a 4A supercell, containing 4 A's arrange in SC.
 * In the middle of each cube there is an icosahedral group of 13 Bs with 
 * one right in the middle of the cell.
 * 
 * Neighbouring icosahedra are oriented oppositly.
 * 
 * See http://www.3dchem.com/inorganicmolecule.asp?id=484 for an image.
 * See various AB13 papers for more info.
 * 
 * @author ajackso1
 * @version $Id: AB13Lattice.java 993 2006-08-25 11:07:06Z anj $
 *
 */
public class AB13Lattice extends MultipleLattice {
	
	/**
	 * Constructor:
	 *   Initialised the presence of two lattice arrays:
	 *
	 */
	public AB13Lattice() {
		// The unit cell, same as SC:
		uc[0] = 1.0;
		uc[1] = 1.0;
		uc[2] = 1.0;
		uc_sep = 1.0;
		this.setUnitCellScaling(2.0,2.0,2.0);
		// Three sites per unit cell:
		nspc = 8*(1+12+1);
		name = "AB13";
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

		// The number of atoms on each lattice:
		int nat = nx*ny*nz/8;
		// Total number of sites:
		na = nspc*nat;
		
		// There will be two lattices:
		lats = new VectorD3[2][];
		lats[0] = new VectorD3[8*nat];
		lats[1] = new VectorD3[8*13*nat];
		
		// Instanciate an SC lattice generator:
		SimpleCubicLattice sc = new SimpleCubicLattice();
		
		// Loop over the 8 unit subcells of the supercell
		sc.setSizeInCells(nx/2,ny/2,nz/2);
		sc.setUnitCellScaling(2.0,2.0,2.0);
		int lcA = 0, lcB = 0;
		for( int ix=0; ix<2; ix++ ) {
			for( int iy=0; iy<2; iy++ ) {
				for( int iz=0; iz<2; iz++ ) {
					// Determine subcell origin and size:
					double cello[] = new double[3], cell[] = new double[3];
					cell[0] = ucx*uc[0]*0.5;
					cell[1] = ucy*uc[1]*0.5;
					cell[2] = ucz*uc[2]*0.5;
					cello[0] = cell[0]*ix;
					cello[1] = cell[1]*iy;
					cello[2] = cell[2]*iz;
					
					// Generate the A lattice per unit subcell:
					sc.setOriginOffset( cello[0], cello[1], cello[2] );
					sc.generate();
					System.arraycopy(sc.getLattice(), 0, lats[0], lcA*nat, nat);				
					lcA++;
					
					// Generate the 13 B lattices per unit subcell:
					// The center particle:
					sc.setOriginOffset( cello[0]+0.5*cell[0], cello[1]+0.5*cell[1], cello[2]+0.5*cell[2] );
					sc.generate();
					System.arraycopy(sc.getLattice(), 0, lats[1], lcB*nat, nat);
					lcB++;
					// Calculate if this is an odd or even site:
					boolean evens = false;
					if( (ix+iy+iz)%2 == 0 ) evens = true;
					// Loop over three directions and place rectangles of B's
					for( int d = 0; d < 3; d++ ) {
						// Loop over the corners of the rectangle:
						for( int ii = 0; ii < 2; ii++ ) {
							for( int jj = 0; jj < 2; jj++ ) {
								// Place the B for the current corner:
								double Bx, By, Bz;
								if( d == 0 ) {
									Bx = cello[0]+0.5*cell[0];
									if( evens ) {
										By = cello[1]+0.5*cell[1] - 4.0*cell[1]/12.0 + ii*8.0*cell[1]/12.0;
										Bz = cello[2]+0.5*cell[2] - 3.0*cell[2]/12.0 + jj*6.0*cell[1]/12.0;
									} else {
										By = cello[1]+0.5*cell[1] - 3.0*cell[1]/12.0 + ii*6.0*cell[1]/12.0;
										Bz = cello[2]+0.5*cell[2] - 4.0*cell[2]/12.0 + jj*8.0*cell[1]/12.0;										
									}
								} else if( d == 1 ) {
									By = cello[1]+0.5*cell[1];
									if( evens ) {
										Bx = cello[0]+0.5*cell[0] - 3.0*cell[0]/12.0 + ii*6.0*cell[0]/12.0;
										Bz = cello[2]+0.5*cell[2] - 4.0*cell[2]/12.0 + jj*8.0*cell[1]/12.0;
									} else {
										Bx = cello[0]+0.5*cell[0] - 4.0*cell[0]/12.0 + ii*8.0*cell[0]/12.0;
										Bz = cello[2]+0.5*cell[2] - 3.0*cell[2]/12.0 + jj*6.0*cell[1]/12.0;										
									}
								} else {
									Bz = cello[2]+0.5*cell[2];
									if( evens ) {
										Bx = cello[0]+0.5*cell[0] - 4.0*cell[0]/12.0 + ii*8.0*cell[0]/12.0;
										By = cello[1]+0.5*cell[1] - 3.0*cell[1]/12.0 + jj*6.0*cell[1]/12.0;
									} else {
										Bx = cello[0]+0.5*cell[0] - 3.0*cell[0]/12.0 + ii*6.0*cell[0]/12.0;
										By = cello[1]+0.5*cell[1] - 4.0*cell[1]/12.0 + jj*8.0*cell[1]/12.0;										
									}
								}
								// Add the lattice:
								sc.setOriginOffset( Bx, By, Bz );
								sc.generate();
								System.arraycopy(sc.getLattice(), 0, lats[1], lcB*nat, nat);
								lcB++;
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public  boolean checkSize(int nx, int ny, int nz) {
    	if( nx%2 != 0 ) return false;
    	if( ny%2 != 0 ) return false;
    	if( nz%2 != 0 ) return false;
    	return true;
    }

}

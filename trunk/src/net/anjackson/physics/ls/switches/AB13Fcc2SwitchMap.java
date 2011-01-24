/**-------------------------------------------------------------
 * jBinLats - AB2Fcc13SwitchMap.java
 * net.anjackson.physics.ls.switches.AB13Fcc2SwitchMap
 * 
 * Created on 25-Jun-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.lattices.AB13Lattice;
import net.anjackson.physics.lattices.FaceCentredCubicABCLattice;
import net.anjackson.physics.lattices.Lattice;

/**
 * This is the AB13 mapping where the AB13 crystal is mapped onto
 * two seperate phases, one of A particles (fcc) and one of 13n B particles (hcp).
 * 
 * NOTE The three-phase approach does not work well, and so this code has been neglected and may not function as it should.
 * 
 * @author ajackso1
 * @version $Id: AB2Fcc2SwitchMap.java 751 2006-05-19 12:08:15Z anj $
 *
 */
public class AB13Fcc2SwitchMap extends ThreePhaseSwitchMap {

	/**
	 * Default constructor, must define npc
	 * IDEA Properly split this into three phases.
	 */
	public AB13Fcc2SwitchMap() {
		npc[0] = 14;
		npc[1] = 14;
	}

	/**
	 * Constructor to create the three lattices for this binary switch.
	 * 
	 */
	public void generate() {
		
		// For this system, there are 8*14 particles per cell:
		int nx = size[0][0];
		int ny = size[0][0];
		int nz = size[0][0];
		int nc = nx*ny*nz;
		n = 8*(1+13)*nc;
		
		// The AB13 Lattice:
		AB13Lattice ab13 = new AB13Lattice();
		ab13.setSizeInCells(nx,ny,nz);
		
		// Stick them in the latt member:
		latt = new VectorD3[2][];
		latt[0] = new VectorD3[n];
		System.arraycopy(ab13.getLattice(0),0,latt[0], 0, 8*nx*ny*nz);
		System.arraycopy(ab13.getLattice(1),0,latt[0], 8*nx*ny*nz, 8*13*nx*ny*nz);
		latt_name[0] = "AB13";
		
		// Build the A and B lattices, A first, then B second:

		// The seperate A and B2 lattices
		FaceCentredCubicABCLattice fcca = new FaceCentredCubicABCLattice();
		fcca.setSizeInCells(2*nx,2*ny,2*nz);
		
		// B lattice requires more care:
		FaceCentredCubicABCLattice fccb = new FaceCentredCubicABCLattice();
		Ntp = 8*13*nx*ny*nz;
		int nbx = (int) Math.round(Math.pow(Ntp,1.0/3.0));
		int nby = nbx;
		int nbz = nbx;
		if( nbx*nby*nbz != Ntp ) {
			// Warn that we cannot create a sensible B-phase with this number of particles:
			System.out.println(" H WARNING: AB13Fcc2SwitchMap cannot create a cubic B-FCC phase with "+(2*nx*ny*nz)+" particles.");
			// Fall back to a system that is simply twice as long in one dimension (z):
			nbx = 2*nx;
			nby = 4*ny;
			nbz = 13*nz;
			System.out.println(" H WARNING: AB13Fcc2SwitchMap creating a phase that is longer in the z-direction: ["+nbx+","+nby+","+nbz+"].");
		}
		// Create the new system at the new size:
		fccb.setSizeInCells(nbx,nby,nbz);

		// Copy the two FCC lattices into the 1th phase:
		latt[1] = new VectorD3[n];
		latt_name[1] = "2.FCC";
		System.arraycopy(fcca.getLattice(), 0, latt[1], 0, 8*nx*ny*nz);
		System.arraycopy(fccb.getLattice(), 0, latt[1], 8*nx*ny*nz, Ntp);

		// Set up the dimensions of the third-phase box:
		boxT.x = fccb.getBasicUnitCell()[0]*nbx;
		boxT.y = fccb.getBasicUnitCell()[1]*nby;
		boxT.z = fccb.getBasicUnitCell()[2]*nbz;

		/* Define distances between neighbouring atoms in the same stacking plane 
		 *  taking advantage of the fact that both lattices have the same Unit Cell: */
		double[] c_uc = new double[3];
		c_uc = fcca.getUnitCell();
		System.out.println(" H DEBUG FCC UC [ "+c_uc[0]+" "+c_uc[1]+" "+c_uc[2]);
		c_uc = ab13.getUnitCell();
		System.out.println(" H DEBUG AB13 UC [ "+c_uc[0]+" "+c_uc[1]+" "+c_uc[2]);
		uc.x = c_uc[0];
		uc.y = c_uc[1];
		uc.z = c_uc[2];
		
		// Put the lattices in the class members, for later use:
		lats[0] = (Lattice)ab13;
		lats[1] = (Lattice)fcca;
		latT = (Lattice)fccb;

	}

	/**
	 * This overrides the mapping of the simulation boxes to take into account
	 * the fact that one of the phases is actually two seperate phases.
	 * 
     * This takes advantage of the fact that all three lattices are
     * the same dimensions, so the scaling and translation to a unit box is identical
     * for all three.
     * 
	 */
	protected void mapToSimulationCube() {
		/* Map the system into the simulation cube, by scaling, translation, and
		 by applying the periodic boundary conditions */
		double xo,yo,zo;
		VectorD3 mapBox;
		/* Loop over lattices: */
		for ( int ilatt = 0; ilatt < 2; ilatt++ ) {
			/* Loop over spheres: */
			for ( int in = 0; in < n; in++ ) {
				if( ilatt == 0 || (in < phaseIndex())) {
					mapBox = boxes[ilatt];
				} else {
					mapBox = boxT;
				}
				/* Copy vector from lattice array */
				xo = latt[ilatt][in].x;
				yo = latt[ilatt][in].y;
				zo = latt[ilatt][in].z;
				
				/*Scale and translate: */
				xo = ( xo/mapBox.x ) - 0.5;
				yo = ( yo/mapBox.y ) - 0.5;
				zo = ( zo/mapBox.z ) - 0.5;
				
				/*Apply periodic boundary conditions: */
				xo = xo - ( (int) ( xo+xo ) );
				yo = yo - ( (int) ( yo+yo ) );
				zo = zo - ( (int) ( zo+zo ) );
				
				/* Copy modified vector back into the lattice array */
				latt[ilatt][in].x = xo;
				latt[ilatt][in].y = yo;
				latt[ilatt][in].z = zo;
				
			}
		}		
	}

}

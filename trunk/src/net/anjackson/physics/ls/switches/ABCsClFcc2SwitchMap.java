/**-------------------------------------------------------------
 * jBinLats - ABCsClFcc2SwitchMap.java
 * net.anjackson.physics.ls.switches.ABCsClFcc2SwitchMap
 * 
 * Created on 22-Mar-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.lattices.ABCsClLattice;
import net.anjackson.physics.lattices.FaceCentredCubicABCLattice;
import net.anjackson.physics.lattices.Lattice;

/**
 * This class creates an AB CsCl to 2 * FCC lattice mapping, with
 * two FCC crystals each of nx*ny*nz particles being mapped onto a single 
 * CsCl structure.
 * 
 * NOTE The three-phase approach does not work well, and so this code has been neglected and may not function as it should.
 * 
 * @author ajackso1
 * @version $Id$
 *
 */
public class ABCsClFcc2SwitchMap extends ThreePhaseSwitchMap {
	
	/**
	 * Default constructor, must define npc
	 */
	public ABCsClFcc2SwitchMap() {
		npc[0] = 2;
		npc[1] = 2;
	}

	/**
	 * 
	 */
	public void generate() {
		
		// For this system, there are two particles per cell:
		int nx = size[0][0];
		int ny = size[0][0];
		int nz = size[0][0];
		int nc = nx*ny*nz;
		n = 2*nc;
		
		// The AB2 Lattice:
		ABCsClLattice ab = new ABCsClLattice();
		ab.setSizeInCells(nx,ny,nz);
		
		// Stick them in the latt member:
		latt = new VectorD3[2][];
		latt[0] = new VectorD3[n];		
		System.arraycopy(ab.getLattice(0),0,latt[0], 0, nx*ny*nz);
		System.arraycopy(ab.getLattice(1),0,latt[0], nx*ny*nz, nx*ny*nz);
		latt_name[0] = "ABCsCl";
		ab.saveAsXYZ("initABCsCl");
		ab.saveAsWRL("initABCsCl");
		
		// Build the A and B lattices, A first, then B second:

		// The seperate A and B2 lattices
		FaceCentredCubicABCLattice fcca = new FaceCentredCubicABCLattice();
		FaceCentredCubicABCLattice fccb = new FaceCentredCubicABCLattice();
		fcca.setSizeInCells(nx,ny,nz);
		
		// Create the new system at the new size:
		Ntp = nx*ny*nx;
		int nbx = nx;
		int nby = nbx;
		int nbz = nbx;
		fccb.setSizeInCells(nbx,nby,nbz);

		// Copy the two FCC lattices into the 1th phase:
		latt[1] = new VectorD3[n];
		latt_name[1] = "2.FCC";
		System.arraycopy(fcca.getLattice(), 0, latt[1], 0, nx*ny*nz);
		System.arraycopy(fccb.getLattice(), 0, latt[1], nx*ny*nz, Ntp);

		// Set up the dimensions of the third-phase box:
		double fccbBox[] = fccb.getSystemSize();
		boxT.x = fccb.getBasicUnitCell()[0]*nx;
		boxT.y = fccb.getBasicUnitCell()[1]*ny;
		boxT.z = fccb.getBasicUnitCell()[2]*nz;
		System.out.println(" DEBUG CsClSM() boxT = ["+fccbBox[0]+","+fccbBox[1]+","+fccbBox[2]+"].");
		System.out.println(" DEBUG CsClSM() boxT = ["+boxT.x+","+boxT.y+","+boxT.z+"].");

		/* Define distances between neighbouring atoms in the same stacking plane 
		 *  taking advantage of the fact that both lattices have the same Unit Cell: */
		double[] c_uc = new double[3];
		c_uc = fcca.getUnitCell();
		System.out.println(" H DEBUG FCC UC [ "+c_uc[0]+" "+c_uc[1]+" "+c_uc[2]);
		c_uc = ab.getUnitCell();
		System.out.println(" H DEBUG AB2 UC [ "+c_uc[0]+" "+c_uc[1]+" "+c_uc[2]);
		uc.x = c_uc[0];
		uc.y = c_uc[1];
		uc.z = c_uc[2];
		
		// Put the lattices in the class members, for later use:
		lats[0] = (Lattice)ab;
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

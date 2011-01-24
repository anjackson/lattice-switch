/**-------------------------------------------------------------
 * jBinLats - SwitchMap.java
 * net.anjackson.physics.ls.switches.SwitchMap
 * 
 * Created on 18-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls.switches;

import java.io.Serializable;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.lattices.Lattice;

/**
 * Base-class for the lattice-switch mappings.
 * All particular switches are children of this.
 * <p>
 * All lattices are stored at unit site seperation.
 * </p>
 * <ul>
 * <li>IDEA The lattice config testing and map generation are very messy and a different pattern of implementation should be used.  I don't have time to do this!</li>
 * </ul>
 * 
 * @author ajackso1
 * @version $Id: SwitchMap.java 996 2006-08-25 13:03:14Z anj $
 * 
 */
public abstract class SwitchMap implements Serializable {

    // Number of sites:
	protected int n;
	
	// Number of particles per unit cell for each lattice:
	// Initialised to values which will fail, to ensure subclasses override them.
	protected int npc[] = new int[] { -1, 2 };
	
	// Dimensions of lattices in unit cells:
	protected int size[][] = new int[2][3];
	
	// Other members:
	protected Lattice[] lats = new Lattice[2];
	protected VectorD3 latt[][] = new VectorD3[2][];
	protected String latt_name[] = new String[2];
	protected VectorD3 boxes[] = { new VectorD3(), new VectorD3() };
	protected VectorD3 uc = new VectorD3();
	protected boolean transformed = false;
	

	/**
	 * Define the size of the lattices (both the same size)
	 * 
	 * @param _nx Number of cells in the x-direction.
	 * @param _ny Number of cells in the y-direction.
	 * @param _nz Number of cells in the z-direction.
	 */
	public void setSize(int _nx, int _ny, int _nz ) {
	   for( int i = 0; i < 2; i++ ) {
		   size[i][0] = _nx;
		   size[i][1] = _ny;
		   size[i][2] = _nz;
	   }
	}
	
	/**
	 * Define the size of the lattices
	 * @param _sizes An [2][3] int array, giving the number of unit cells in each direction for the two lattices.
	 */
	public void setSizes( int[][] _sizes ) {
	   for( int i = 0; i < 2; i++ ) {
		   for( int d = 0; d < 3; d++ ) {
			   size[i][d] = _sizes[i][d];
		   }
	   }
	}
	
	/**
	 * Get the vectors describing the i'th lattice.
	 * @param i the index of the lattice of interest.
	 * @return the arrat of VectorD3 lattice vectors.
	 */
	public VectorD3[] getLatticeVectors(int i) {
		if( ! transformed ) this.mapToSimulationBox();
		return latt[i];
	}
	
	/**
	 * Get the textual name of the i'th lattice.
	 * @param i the index of the lattice of interest.
	 * @return The name of the lattice - a short text string.
	 */
	public String getLatticeName( int i ) {
		if( ! transformed ) this.mapToSimulationBox();
		return latt_name[i];
	}
	
	/**
	 * Get the overall size of the PBC simulation box:
	 * @return A VectorD3 describing the size of the simulation box.
	 */
	public VectorD3 getBox(int i) {
		if( ! transformed ) this.mapToSimulationBox();
		return boxes[i];
	}

	/**
	 * 
	 * This method sets up the simulation box, and initiates the mapping to 
	 * run-time coordinates.
	 * 
	 * Should be called automatically when a lattice is requested.
	 * 
	 * This has been made synchronised, to avoid a race-condition on 'transformed'.
	 * This is not likely to be a problem as the code is not intended for 
	 * multi-threaded use.
	 * 
	 */
	private synchronized void mapToSimulationBox() {
		// Do not do this twice!
		if( transformed ) return;
		transformed = true;
		
		/* Calc the cell size at diameter = 1.0 */
		for( int il = 0; il < 2; il++ ) {
			boxes[il].x = lats[il].getSystemSize()[0];
			boxes[il].y = lats[il].getSystemSize()[1];
			boxes[il].z = lats[il].getSystemSize()[2];
		}
		
		/* Resize the box to get the simulation cube */
		mapToSimulationCube();
		
	}
	
	/**
	 * This method takes the lattice defined in the constructor, and maps it onto
	 * a unit cell centred about 0,0 stretching from -0.5 to +0.5, for the sim.
	 * Also defines the relationship between the density and the true box size.
	 * 
	 */
	protected void mapToSimulationCube() {
		
		System.out.println("Got [0]x="+latt[0][0].x);
		System.out.println("Got [1]x="+latt[1][0].x);

		/* Map the system into the simulation cube, by scaling, translation, and
		 by applying the periodic boundary conditions */
		double xo,yo,zo;
		/* Loop over lattices: */
		for ( int ilatt = 0; ilatt < 2; ilatt++ ) {
			System.out.println(" TEST GOT boxes["+ilatt+"].x = "+boxes[ilatt].x);
			System.out.println(" TEST GOT boxes["+ilatt+"].y = "+boxes[ilatt].y);
			System.out.println(" TEST GOT boxes["+ilatt+"].z = "+boxes[ilatt].z);

			/* Loop over spheres: */
			for ( int in = 0; in < n; in++ ) {
				/* Copy vector from lattice array */
				xo = latt[ilatt][in].x;
				yo = latt[ilatt][in].y;
				zo = latt[ilatt][in].z;
				
				/*Scale and translate: */
				xo = ( xo/boxes[ilatt].x ) - 0.5;
				yo = ( yo/boxes[ilatt].y ) - 0.5;
				zo = ( zo/boxes[ilatt].z ) - 0.5;
				
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
		
		System.out.println("Generated [0]x="+latt[0][0].x);
		System.out.println("Generated [1]x="+latt[1][0].x);
		
	}
	
	/**
	 * Finds the largest 3 factors that make up a target N.
	 * 
	 * Useful in principle, but cannot cope with re-jigging the
	 * results to find the right arrangement to fit the PBC/UC
	 * requirements. i.e. if the system has to be multiples of
	 * particular amounts in particular directions (as for the 
	 * CP structures).
	 *
	 * @param target Desired product of 3 factors.
	 * @return An int[3] array such that i[0]*i[1]*i[2] = target.
	 */
	protected int[] findLargest3Factors( int target ) {
		int[] best = new int[3];
		
		// Brute-force search for the best factors:
		for( int a = 1; a <= n; a++ ) {
			for( int b = a; b <= n; b++ ) {				
				for( int c = b; c <= n; c++ ) {
					if( a*b*c == target ) {
						best[0] = a;
						best[1] = b;
						best[2] = c;
					}
				}
			}
		}
		// Return the best we have found:
		return best;
	}
	
	
	/**
	 * The actual routine that initializes and generates the switch map.
	 */
	public void init() {
		// Perform sanity checks:
		check();
		// Generate the actual lattice:
		generate();
	}
	
	/** 
	 * Check the requested sizes are sane:
	 */
	public void check() {
		int ns[] = getTotalNs();
		// Check N. dof match up:
		if( ns[0] != ns[1] ) {
			System.out.println("The lattices have different numbers of particles! "+ns[0]+" != "+ns[1]);
			System.exit(1);
		}
	}
	

	/**
	 * Calculate how many particles will be in each lattice.
	 * Should be overridden when complex cells are used?
	 * @return an integer array containing the total N for l0 and l1.
	 */
	public int[] getTotalNs() {
		int ns[] = new int[2];
		for( int i = 0; i < 2; i++ ) {
			ns[i] = npc[i]*size[i][0]*size[i][1]*size[i][2];
		}
		return ns;
	}
	

	
	/**
	 * Create the actual lattices and the mapping.
	 */
	public abstract void generate();
}

/**
 * 
 */
package net.anjackson.physics.lattices;

import net.anjackson.maths.*;
import java.io.*;

/**
 * IDEA Make this create a lattice based on CLI args, and output an XYZ file suitable for plotting in jMol/iMol.
 * 
 * @author ajackso1
 * @version $Id: GenerateXYZFile.java 981 2006-08-21 14:40:43Z anj $
 *
 */
public class GenerateXYZFile {

	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// Size - Old LS was 6 [216], 12 [1728] or 18 [5832].
		int nx = 12, ny = nx, nz = nx;
		// Create an HCP lattice
		HexagonalClosePackedLattice hcp = new HexagonalClosePackedLattice();
		hcp.setSizeInCells(nx,ny,nz);
		hcp.generate();
		// Output it to a file:
		outputXYZFile( hcp, "hcp.xyz");
		
		// Create an FCC lattice
		FaceCentredCubicABCLattice fcc = new FaceCentredCubicABCLattice();
		fcc.setSizeInCells(nx,ny,nz);
		fcc.generate();
		// Output it to a file:
		outputXYZFile( fcc, "fcc.xyz");
		
		// Create an AB2 lattice
		AB2Lattice ab2 = new AB2Lattice();
		ab2.setSizeInCells(nx,ny,nz);
		ab2.generate();
		// Output it to a file:
		outputXYZFile( ab2, "ab2.xyz");		
	}
	
	/**
	 * Take a lattice and output it as a simple XYZ file.
	 * 
	 * @param lat
	 * @param filename
	 */
	public static void outputXYZFile( Lattice lat, String filename ) {
		try {
		// Open the file:
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		// Output header:
		out.write(" "+lat.getNumberOfSites()+"\n");
		out.write(" \n");
		// Output each point:
		VectorD3[][] lv = null;
		if( lat.getNumberOfLattices() == 1 ) {
			lv = new VectorD3[1][];
			lv[0] = lat.getLattice();
		} else {
			MultipleLattice mlat = (MultipleLattice) lat;
			lv = mlat.getLattices();
		}
		// Output the lattices:
		char atom_char = 'X';
		for( int l = 0; l < lat.getNumberOfLattices(); l++ ) {
			for( int i = 0; i < lv[l].length; i++ ) {				
				if( l == 0 ) {
					atom_char = 'C';
				} else {
					atom_char = 'X';
				}
				out.write( atom_char+" "+lv[l][i].x+" "+lv[l][i].y+" "+lv[l][i].z+"\n");
			}
		}
		// Close the file:
		out.close();		
		} catch( IOException e ) {
			// Report this error!
			System.err.println("Writing to file "+filename+" failed!");
			System.err.println(e);
		}
	}

}

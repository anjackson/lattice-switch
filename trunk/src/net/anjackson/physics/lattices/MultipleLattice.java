/**-------------------------------------------------------------
 * jBinLats - MultipleLattice.java
 * net.anjackson.physics.lattices.MultipleLattice
 * 
 * Created on 07-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence:
 * -------------------------------------------------------------
 */
package net.anjackson.physics.lattices;

import java.io.FileOutputStream;
import java.io.PrintStream;

import net.anjackson.maths.*;

/**
 * This is the abstract base class for any lattice with more than one atom
 * per unit cell.  i.e. there are two or more lattices overlaid and spanning 
 * the same space.
 * 
 * 
 * @author ajackso1
 * @version $Id: MultipleLattice.java 893 2006-06-29 14:08:46Z anj $
 *
 */
public abstract class MultipleLattice extends Lattice {
    // Array of vectors for holding multiple lattices:
	VectorD3[][] lats = null;
	
	/**
	 * Get the complete lattice array
	 */
	public VectorD3[] getLattice() {
		// Sanity Check:
		if( lats == null ) generate();
		if( lats == null ) return null;
		if( lats.length == 0 ) return null;
		// Assemble:
		int na_tot = 0;
		// Count the total number of sites:
		for( int i = 0; i < lats.length; i++ ) na_tot += lats[i].length;
		// Create the array
		lat = new VectorD3[na_tot];
		// Copy the lattices into the overall lat:
		int na = 0;
		for( int i = 0; i < lats.length; i++ ) {
			System.arraycopy(lats[i],0, lat, na, lats[i].length);
			na += lats[i].length;
		}
		// Return:
		return super.getLattice();
	}
	
	/**
	 * How many single-particle lattices make up this multi-lattice?
	 * @return The number of lattices.
	 */
	public int getNumberOfLattices() {
		if( lats == null ) return 0;
		return lats.length;
	}
	/**
	 * Get the points for one of the lattices:
	 * @param i the index of the lattice to be recovered.
	 * @return The array of VectorD3 vectors that describes that lattice.
	 */
	public VectorD3[] getLattice( int i ) {
		if( lats == null ) generate();
		if( lats == null ) {
			System.out.println("MultipleLattice: ERROR: lats == null");
			return null;
		}
		if( i > lats.length ) {
			System.out.println("MultipleLattice: ERROR: There is no lattice "+i+" of ["+lats.length+"].");
			return null;
		}
		if( lats[i] == null ) {
			System.out.println("MultipleLattice: ERROR: Lattice "+i+" of ["+lats.length+"] is set to null!");
			return null;
		}
		// Lattice is okay - returning it:
		return lats[i];
	}

	/**
	 * Returns the array of lattice-site vectors.
	 * @return The array of lattices, each being an array of VectorD3s.
	 */
	public VectorD3[][] getLattices() {
		if( lats == null ) generate();
		return lats;
	}
	
	/**
	 * Save this lattice as an XYZ file:
	 * @param filename The name of the file to save into, as a String.
	 */
	public void saveAsXYZ(String filename) {
		char atomChar = 'C';
		
		/* open file and outpout number of atoms to come */
		String xyzfn = filename+".xyz";
		try {
			PrintStream out = new PrintStream(new FileOutputStream(xyzfn));
			out.print(" "+(na)+"\n\n");
			
			/* loop over all particles and output as carbons */
			for( int il = 0; il < lats.length; il++ ) {
				if( il == 0 ) {
					atomChar = 'C';
				} else {
					atomChar = 'H';
				}
				for( int i = 0; i < lats[il].length; i++ ) {
					out.println(atomChar+" "+lats[il][i].x+" "+lats[il][i].y+" "+lats[il][i].z);
				}
			}
			
			/* close the file */
			out.close();
		} catch ( Exception e ) {
			System.out.println("save_config_xyz("+filename+") Failed on exception: "+e);
		}
		
	}
	/**
	 * Save this lattice as a VRML (.wrl) file:
	 * @param filename File name to save the WRL to.
	 */
    public void saveAsWRL(String filename) {
    	double radius = 0.5;
    	double box[] = this.getSystemSize();
    	
		/* Output as a WRL file */
		String wrlfn = filename+".wrl";
		try {
			PrintStream out = new PrintStream(new FileOutputStream(wrlfn));
			
			/* Output header */
			out.print("#VRML V2.0 utf8\n");
			out.print("# By Andy Jackson (2006)\n");
			out.print("# META \"generator\" MCSim by Andrew Jackson\n");
			out.print("WorldInfo {\n");
			out.print("  title \"MCSim configuration\"\n");
			out.print("}\n");
			out.print("Viewpoint {\n");
			out.print("  description \"Side View\"\n");
			out.print("  position 0.0 0.0 30.0\n");
			out.print("}\n");
			
			/* Output the simulation box */
		    out.print("Transform { \n");
		    out.print("  translation "+0.5*box[0]+" "+0.5*box[1]+" "+0.5*box[2]+"\n");
		    out.print("  children Shape {\n");
		    out.print("    appearance Appearance {\n");
		    out.print("      material Material {\n");
		    out.print("        diffuseColor 0.8 0.8 0.8\n");
		    out.print("        transparency 0.75\n");
		    out.print("      }\n");
		    out.print("    }\n");
		    out.print("    geometry Box {\n");
		    out.print("      size "+box[0]+" "+box[1]+" "+box[2]+"\n");
		    out.print("    }\n");
		    out.print("  }\n");
		    out.print("}\n");
		    
			/* loop over all particles and output as carbons */
		    double sp[] = new double[3];
			for( int il = 0; il < lats.length; il++ ) {
				if( il == 0 ) {
					radius = 0.5;
				} else {
					radius = 0.5*5.0/12.0;
				}
				for( int i = 0; i < lats[il].length; i++ ) {
					// Calc sphere coords:
					sp[0] = lats[il][i].x;
					sp[1] = lats[il][i].y;
					sp[2] = lats[il][i].z;
					// Output Sphere:
					out.print("  Transform {\n");
					out.print("    translation "+sp[0]+" "+sp[1]+" "+sp[2]+"\n");
					out.print("    children Shape {\n");
					out.print("      appearance Appearance {\n");
					out.print("        material Material {\n");
					out.print("        diffuseColor 1.0 1.0 1.0\n");
					out.print("          transparency 0.0\n");
					out.print("        }\n");
					out.print("      }\n");
					out.print("      geometry Sphere {\n");
					out.print("        radius "+radius+"\n");
					out.print("      }\n");
					out.print("    }\n");
					out.print("  }\n");
				}
			}
						
			/* close the file */
			out.close();
		} catch ( Exception e ) {
			System.out.println("save_config_xyz("+wrlfn+") Failed on exception: "+e);
		}
		
    }
}

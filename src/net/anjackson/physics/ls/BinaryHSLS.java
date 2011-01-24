/**-------------------------------------------------------------
 * jBinLats - BinaryHSLS.java
 * net.anjackson.physics.ls.BinaryLS
 * 
 * Created on 08-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

import java.io.FileOutputStream;
import java.io.PrintStream;

import net.anjackson.maths.VectorD3;
import net.anjackson.physics.ls.switches.BinarySwitchMap;
import net.anjackson.utils.ParamFileReader;

/**
 * This is the binary hard-core simulation.  It overrides the monodisperse
 * case and overrides the methods needed to simulate the binary switch.
 * <p>
 * Note that there are three ways to do this.  You can have:
 * <ul>
 *  <li>ABx switching to Afcc and xBfcc seperate simulation boxes.</li>
 *  <li>ABx switching to Afcc and xBfcc coexistance (in the same box).</li>
 *  <li>ABx switching to (1+x)A pure phase (one-box to one-box).</li>
 * </ul>
 * This code implements the latter, one-box method.
 * </p>
 * <p>
 * Therefore, we combine a global transformation and a radii change into a single 
 * Hamiltonian switch.
 * </p>
 * <p>
 * This has the advantage of being essentially the same as the simple LS code. We work with
 * the same ensembles in the same way.  The ensembles are:
 * <ul>
 *  <li>NVT - in which case the boxes in the seperate phases are changed so that the densities match</li>
 *  <li>NPT - a single box sise is used, and the density is initialised in the initial phase</li>
 * </ul>
 * This behaviour should be inherited from the superclass.
 * </p>
 * <p>
 * By convention, phase zero is the composite phase (AB2) and phase 1 is made up of two phases.
 * </p>
 * <ul>
 * <li>IDEA TEST Set up FCC-to-2FCC transform, to check that part of things for 'contants.' But this is probably meaningless!</li>
 * </ul>
 * </p>
 * @author ajackso1
 * @version $Id: BinaryHSLS.java 980 2006-08-21 13:22:35Z anj $
 *
 */
public class BinaryHSLS extends HSLS {
	
	// Specializations are private as these should only be used in this class.	
	// The hard-sphere radii definitions:
	/** The ratio of diameters */
	private double rad_ratio = 1.0;
	/** The diameters, for the B particles */
	double diameterB = 1.0;
	/** The diameterB^2, for the B-B interactions */
	double diameterB2 = 1.0;
	/** The (diameterB/2+diameterA/2)^2, for the A-B interactions */
	double diameterAB2 = 1.0;	
	
	/** The point in the lattice arrays where we switch to B particles (A's come first) */
	private int iB = -1;
	/** The number of A particles: */
	private int nA = -1;
	/** The number of B particles: */
	private int nB = -1;
	/** The ratio of B particles to A particles: */
	private int nBoA = -1;
	
	/**
	 * read parameters, both general and model-specific:
	 */
	protected void read_params_from_user( ParamFileReader pars ) {
		// Read general parameters:
		super.read_params_from_user(pars);
		
		// Read model parameters for split-system:
		try {
			// Read in the binary ratio:
			rad_ratio = pars.getDouble("alpha");
		} catch( ParamFileReader.UndefinedParameterException e ) {
			// If any of them is not defined, cry and grumble.
			throw( new LSSimException( LSSimException.PARAM,e));
		}

	}

	/**
	 * Initialise both general and split-system variables:
	 */
	protected void init_variables() {		
		// Cast the switch map as a binary one:
	    BinarySwitchMap bsm = (BinarySwitchMap)lsm;
	    
		// Check where the B particles start from in the latt arrays:
		iB = bsm.phaseIndex();
		if( iB < 0 ) {
			throw( new LSSimException( LSSimException.PARAM,
					"No or invalid index supplied for the start of the B particle array: "+iB));
		}
		nA = iB;
		nB = bsm.getThirdN();
		nBoA = nB/nA;
		System.out.println(" H BinaryHSLS init: nB = "+nB+" start@ "+iB+"/"+(nA+nB)+" nBoA = "+nBoA);
				
		// Derived utility variables:
		setDiameterB(diameter*rad_ratio); // FCC: B Ratio is defined as rB/rA
		
		// Initialise general system variables:
		super.init_variables();

		// TODO Polyd stuff in here:  set up parts of diameter array.  B-particles in both phases need fixing.  densities need updating too?
				
	}
	
	protected void initBoxes() {
		// Do binary-system init:
		// Adjusts to fix the density in the mixed phase:
		/* Define the (initial) size of the simulation boxes:
		 * 
		 * The super-classes init_variables has already set up the density of A-FCC.
		 * 
		 * The AB phase box must be adjusted for the fixed particle diameter (cube-root of density in a monophase).
		 * 
		 * First, calculate the volume ratio V_1/V_0, derived from
		 *   \rho_0 = (N_a.\sigma_a^3 + N_b.\sigma_b^3)/(\sqrt(2).V_1)
		 *   \rho_1 = N.\sigma^3/(\sqrt(2).V_1)
		 * for \rho_0 = \rho_1 and \sigma = \sigma_a
		 * 
		 */
		double volumeRatio = (1.0+nBoA*Math.pow(rad_ratio, 3.0))/(1.0+nBoA);
		// Now adjust this scaling to take into account the relative volumes of the two systems:
		volumeRatio = volumeRatio/((box[0].x*box[0].y*box[0].z)/(box[1].x*box[1].y*box[1].z));
		// And find the linear factor corresponding to this volume scaling:
		double boxScaler = Math.pow( volumeRatio, 1.0/3.0 );
		// Alter the AB2 density by altering the box size.
		box[0] = box[0].times(boxScaler);
		setBoxes(box);
		
		// Do the vanilla init:
		super.initBoxes();
				
	}
	
	protected void setDiameterB( double _diameterB ) {
		diameterB = _diameterB;
		diameterB2 = diameterB*diameterB;
		// A-B interactions:
		diameterAB2 = 0.5*(diameter + diameterB);
		diameterAB2 = diameterAB2*diameterAB2; // AB2: A-B		
	}
	
	/**
	 * Write out the header, including binary-specific stuff.
	 */
	protected void write_output_header() {
		System.out.println(" H Binary hard sphere lattice-switch monte carlo code: AB2 v Afcc, Bfcc:");
		System.out.println(" H ");
		// Write out general system parameters:
		super.write_output_header();
		// Write out BinaryHS-specific parameters:
		System.out.println(" H ");
		System.out.println(" H The diameter ratio alpha = d(B)/d(A) = "+rad_ratio+" > == "+diameterB/diameter+" == "+Math.sqrt(diameterB2/diameter2));
		System.out.println(" H The diameters are "+diameter+" for A and "+diameterB+" for B particles.");
		System.out.println(" H The number of A particles is "+nA+" out of a total of "+n);
		System.out.println(" H The number of B particles is "+nB+" starting at "+iB);
		System.out.println(" H ");
		// Extra debugging information:
		System.out.println(" H The A diameters are ["+diameter+","+diameter2+"].");
		System.out.println(" H The B/AB diameters are [ "+diameterB+" "+diameterB2+" "+diameterAB2+" ].");
		System.out.println(" H The size of the box for the 0 phase is ["+box[0].x+","+box[0].y+","+box[0].z+"].");
		System.out.println(" H The size of the box for the 1 phase is ["+box[1].x+","+box[1].y+","+box[1].z+"].");
		System.out.println(" H VolFrac AB2 = "+calc_vol_frac(0));
		System.out.println(" H VolFrac A in AB2 = "+calc_partial_vol_frac(0));
		System.out.println(" H VolFrac B in AB2 = "+calc_partial_vol_frac(1));
		System.out.println(" H VolFrac A = "+calc_vol_frac(1));
		System.out.println(" H ");
				
	}

	
	
	/**
	 * ij_inter:
	 *   Calculate the overlap for a hard-sphere system, as function of dr^2.
	 *   If index&lt;iB then we have an A-particle, &gt; iB is a B-particle.
	 */
	protected double ij_inter(double dr2, int i, int j, int il ) {
		double testdiam2;

		// If A-A interaction:
		if( il == 1 || (i < iB && j < iB) ) {
			testdiam2 = diameter2;
		} else if( i >= iB && j >= iB ) {
			// B-B interaction:
			testdiam2 = diameterB2;
		} else {
			// A-B interaction:
			testdiam2 = diameterAB2;
		}
		
		// Check for an overlap:
	    if (dr2 < testdiam2 ) return(1.0);
	    
	    // Otherwise:
    	return(0.0);

	}
	
	/**
	 * The actual pair-potential interaction in polydisperse systems:
	 * Calculates the pair potential as a function of dr^2 and pdiameter[il][n].
	 */
	protected double ij_inter_poly(double dr2, int i, int j, int il) {
		throw new LSSimException(LSSimException.PARAM,"The BinaryHSLS class does not support polydisperse calculations at present!");
	}


	/**
	 * Output extra information pertaining to this binary system.
	 */
	protected void report_simulation(int sweep) {
		// Perform general reporting:
		super.report_simulation(sweep);
	}
	
	
	/**
	 * Perform extra checks and output for this system
	 */
	protected void check_simulation(int sweep) {
		// Perform all the usual checks:
		super.check_simulation(sweep);
	}

	/**
	 * Output extra information pertaining to this binary system.
	 */
	protected void output_sim_info_tail() {
	}
	
	/**
	 * Calculate the density of the AB and pure A phases:
	 */
	protected double calc_density() {
		if( c_lat == 0 ) {
			//System.out.println("Bin.calc_density(0) "+lsize[0][0]+" "+lsize[0][1]+" "+lsize[0][2]+" "+diameter+" "+diameterB+" "+box[c_lat].x+" "+box[c_lat].y+" "+box[c_lat].z);
			return(nA*(Math.pow(diameter,3.0) + nBoA*Math.pow(diameterB,3.0))/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));
		} else {
			//System.out.println("Bin.calc_density(1) "+lsize[1][0]+" "+lsize[1][1]+" "+lsize[1][2]+" "+diameter+" "+diameterB+" "+box[c_lat].x+" "+box[c_lat].y+" "+box[c_lat].z);
			return(n*Math.pow(diameter,3.0)/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));
			//return(nx*ny*nz*(Math.pow(diameter,3.0)+nBoA*Math.pow(diameter,3.0))/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));
		}
	}
	
     /**
      * Calculate the volume fraction of the phases
      * 
      * @param phase [0 = AB phase, 1 = A phase]
      * @return The volume fraction of the specified phase.
      */
	protected double calc_vol_frac(int phase) {
		if( phase == 0 ) {
			return ( (nA*sphere_volume(diameter) + nB*sphere_volume(diameterB)) )/(box[0].x*box[0].y*box[0].z);
		} else if( phase == 1 ) {
			return ( n*sphere_volume(diameter) )/(box[1].x*box[1].y*box[1].z);
		} else {
			// Fail with a runtime error:
			throw( new LSSimException(LSSimException.PARAM, "calc_vol_frac: Unknown phase number"+phase ));
		}
	}

	/**
	 * Calculate the partial volume fractions in the AB phase:
	 * @param phase [0 = A in AB, 1 = B in AB]
	 * @return The partial volume fraction such that Vf = Vf(A) + Vf(B).
	 */
	protected double calc_partial_vol_frac(int phase) {
		if( phase == 0 ) {
			return (nA*sphere_volume(diameter) )/(box[0].x*box[0].y*box[0].z);
		} else if( phase == 1 ) {
			return (nB*sphere_volume(diameterB) )/(box[0].x*box[0].y*box[0].z);
		} else {
			// Fail with a runtime error:
			throw( new LSSimException(LSSimException.PARAM, "calc_partial_vol_frac: Unknown phase number"+phase ));
		}
	}

	/**
	 * Overriding this function to alter the radii.
	 */
	protected void save_config_xyz(int il, VectorD3[] lattice, VectorD3 _box, int start, int number, String filename) {
		int i,ii,ij,ik;
		char atomChar = 'C';
		boolean output_xyz = false;
		
		if( output_xyz ) {
        /* open file and outpout number of atoms to come */
		String xyzfn = filename+".xyz";
		try {
			PrintStream out = new PrintStream(new FileOutputStream(xyzfn));
			out.print(" "+(number+8)+"\n\n");
			
			/* loop over all particles and output as carbons */
			for ( i = start; i < start+number; i++ ) {
				if( i<iB ) {
					atomChar = 'H';
				} else {
					atomChar = 'X';
				}
				out.println(atomChar+" "+(lattice[i].x+disp[i].x)*_box.x+" "+(lattice[i].y+disp[i].y)*_box.y+" "+(lattice[i].z+disp[i].z)*_box.z);
			}
			
			/* Loop over the corners of the pbc box, outputting as hydrogens */
			for ( ii = -1; ii <= 1; ii = ii + 2 ) {
				for ( ij = -1; ij <= 1; ij = ij + 2 ) {
					for ( ik = -1; ik <= 1; ik = ik + 2 ) {
						out.println("C "+_box.x*(ii*0.5)+" "+_box.y*(ij*0.5)+" "+_box.z*(ik*0.5));
					}
				}
			}
			
			/* close the file */
			out.close();
		} catch ( Exception e ) {
			System.out.println("save_config_xyz("+xyzfn+") Failed on exception: "+e);
		}
		}
	
	/****************************************************************/
	
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
	    out.print("  translation 0.0 0.0 0.0\n");
	    out.print("  children Shape {\n");
	    out.print("    appearance Appearance {\n");
	    out.print("      material Material {\n");
	    out.print("        diffuseColor 0.8 0.8 0.8\n");
	    out.print("        transparency 0.75\n");
	    out.print("      }\n");
	    out.print("    }\n");
	    out.print("    geometry Box {\n");
	    out.print("      size "+_box.x+" "+_box.y+" "+_box.z+"\n");
	    out.print("    }\n");
	    out.print("  }\n");
	    out.print("}\n");
	    
		/* loop over all particles and output as carbons */
	    double sp[] = new double[3];
		for ( i = start; i < start+number; i++ ) {
			// Calc sphere coords:
			sp[0] = (lattice[i].x+disp[i].x)*_box.x;
			sp[1] = (lattice[i].y+disp[i].y)*_box.y;
			sp[2] = (lattice[i].z+disp[i].z)*_box.z;
			// Output Sphere:
            out.print("  Transform {\n");
            out.print("    translation "+sp[0]+" "+sp[1]+" "+sp[2]+"\n");
            out.print("    children Shape {\n");
            out.print("      appearance Appearance {\n");
            out.print("        material Material {\n");
            if( i < iB ) {
              out.print("        diffuseColor 1.0 1.0 1.0\n");
            } else {
              out.print("        diffuseColor 1.0 0.0 0.0\n");            	
            }
            out.print("          transparency 0.0\n");
            out.print("        }\n");
            out.print("      }\n");
            out.print("      geometry Sphere {\n");
            // Choose the diameter based on the lattice ID and particle ID:
            if( il == 1 || i < iB ) {
              out.print("        radius "+0.5*diameter+"\n");
            } else {
              out.print("        radius "+0.5*diameterB+"\n");
            }
            out.print("      }\n");
            out.print("    }\n");
            out.print("  }\n");
		}
					
		/* close the file */
		out.close();
	} catch ( Exception e ) {
		System.out.println("save_config_xyz("+wrlfn+") Failed on exception: "+e);
	}
	}
	
	/**
	 * calculate the volume fraction for monodisperse hard spheres
	 */
	protected double calcVolFrac() {
		return calc_vol_frac(c_lat);
	}

}

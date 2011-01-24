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
import net.anjackson.physics.ls.switches.ThreePhaseSwitchMap;
import net.anjackson.utils.BlockAnalysis;
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
 * This code implements the former, three-box method.
 * </p>
 * <p>
 * These two different box rules require different overrides, and so the 
 * building of the NN lists and application of the box dimensions/radii are overridden.
 * </p>
 * <p>
 * So, the particle-move code is the same as before, with a little tweaking of the 
 * interations (two-radii) and the neighbour lists (doubled FCC phase simulations are independant
 * but live in the same periodic simulation box.  The original volume-fluctuation code is the same, and
 * controls the AB2 and A-FCC volumes.  The extra degrees of freedom concerning the B-FCC volume is implemented
 * as an extra MC volume move that is a non-interacting 'ghost' in the AB2 phase.
 * </p>
 * <p>
 * The radii are currently fixed to be the same in both phases, and the AB2 and A-FCC share the same volume.
 * This means the B particles will overlap when the initial density is set too high, and a fixed-density
 * simulation is mapping between rather arbitrary system densities.
 * The extra DOF of the B-FCC volume is also set to reproduce the same density as the other two.  
 * Of course, in an NPT simulation this is just an initial condition and the system is
 * free to explore other densities.  This also applies to the AB2 and A-FCC phases.
 * </p>
 * <p>
 * By convention, phase zero is the composite phase (AB2) and phase 1 is made up of two phases.
 * </p>
 * <p>
 * Due to the independent centre-of-mass motion of the two seperate phases, this switch is unstable unless
 * the displacements are locked into a fixed distance around the lattice sites.  Otherwise, the mappings drift.
 * <p>
 * <ul>
 * <li>IDEA Add radius and polydispersity to the third-phase for this binary system.</li>
 * <li>IDEA Could invert displacements of odd particles in the mixed phase to limit effects of CoM motion.</li>
 * </ul>
 * </p>
 * @author ajackso1
 * @version $Id: BinaryHSLS.java 907 2006-07-12 15:35:04Z anj $
 *
 */
public class Binary3PhaseHSLS extends HSLS {
	// The lattice switch-map, cast to three-phase form:
	ThreePhaseSwitchMap tpsm = null;
	
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
	/** The centre-of-mass of the A particles: */
	/*
	private VectorD3 com_A;
	*/
	/** The centre-of-mass of the B particles: */
	/*
	private VectorD3 com_B;
	*/
	
	// The second simulation box, for the B-particles, a ghost in AB2:
	private VectorD3 init_boxB = new VectorD3();
	private VectorD3 boxB = new VectorD3();
	//private VectorD3 box2B = new VectorD3();
	private double[] box2B = new double[3];
	// Block analysis of the volume for the B phase.
	private BlockAnalysis bDensblk = new BlockAnalysis();
	
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
		// Try to cast the SwitchMap as a three-phase switchmap:
		tpsm = (ThreePhaseSwitchMap) lsm;
		init_boxB = tpsm.getThirdBox();
		
		// Set B-box to the unit lattice while NN lists are constructed:
		setBoxB(init_boxB);
	
		// Check where the B particles start from in the latt arrays:
		iB = tpsm.phaseIndex();
		if( iB < 0 ) {
			throw( new LSSimException( LSSimException.PARAM,
					"No or invalid index supplied for the start of the B particle array: "+iB));
		}
		nA = iB;
		nB = tpsm.getThirdN();
		nBoA = nB/nA;
		System.out.println(" D The B-phase has been initialised at "+nB+" particles starting at "+iB);
		
		// Initialise the CoM recorders:
		/*
		com_A = new VectorD3(0,0,0);
		com_B = new VectorD3(0,0,0);
		*/
		
		// Initialise general system variables:
		super.init_variables();
		
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
		double volumeRatio = (1.0+nBoA*Math.pow(rad_ratio, 3.0));
		// Now adjust this scaling to take into account the relative volumes of the two systems:
		volumeRatio = volumeRatio/(box[0].x*box[0].y*box[0].z/(box[1].x*box[1].y*box[1].z));
		// And find the linear factor corresponding to this volume scaling:
		double boxScaler = Math.pow( volumeRatio, 1.0/3.0 );
		System.out.println(" DEBUG Using the box scalar = "+boxScaler+" nBoA = "+nBoA+" alpha = "+rad_ratio);
		// Alter the AB2 density by altering the box size.
		box[0] = box[0].times(boxScaler);
		setBoxes(box);
		
		// Make the initial /density/ of the B-crystal the same as the A-crystal, by shrinking it's box:
		setBoxB(init_boxB.times(rad_ratio));
		
		// Derived utility variables:
		setDiameterB(diameter*rad_ratio); // FCC: B Ratio is defined as rB/rA
		
		// Calculate the initial densities and check they are equal - only valid in NVT:
		if( ! NPT_SIM ) {
			int tmp_lat = c_lat, dens_errors = 0;
			double denses[] = new double[3];
			c_lat = 0;
			denses[0] = this.calc_density();
			c_lat = 1;
			denses[1] = this.calc_density();
			c_lat = tmp_lat;
			denses[2] = this.calc_B_density();
			if( Math.abs(denses[0]-denses[1]) > E_TOL ) {
				System.out.println(" ERROR: Densities in AB and A phase are different! AB@"+denses[0]+" A@"+denses[1]);
				dens_errors++;
			}
			if( Math.abs(denses[0]-denses[2]) > E_TOL ) {
				System.out.println(" ERROR: Densities in AB and B phase are different! AB@"+denses[0]+" B@"+denses[2]);
				dens_errors++;
			}
			if( Math.abs(denses[1]-denses[2]) > E_TOL ) {
				System.out.println(" ERROR: Densities in A and B phase are different! A@"+denses[1]+" B@"+denses[2]);
				dens_errors++;
			}
			if( dens_errors > 0 ) {
				System.exit(1);
			} else {
				System.out.println(" CHECK: All densities checked as equal: AB@"+denses[0]+" A@"+denses[1]+" B@"+denses[2]);
			}
		}
		// IDEA Polyd stuff in here:  set up parts of diameter array.  B-particles in both phases need fixing.  densities need updating too?
		
		// Save the initial A and B phase configurations:
		save_config_xyz(1, latt[1], box[1], 0, iB, 1+"-"+latt_name[1]+"-A-phase");
		save_config_xyz(1, latt[1], boxB, iB, nB, 1+"-"+latt_name[1]+"-B-phase");
		
	}
	
	protected void setDiameterB( double _diameterB ) {
		diameterB = _diameterB;
		diameterB2 = diameterB*diameterB;
		// A-B interactions:
		diameterAB2 = 0.5*(diameter + diameterB);
		diameterAB2 = diameterAB2*diameterAB2; // AB2: A-B		
	}
	
	protected void setBoxB( VectorD3 _boxb ) {
		boxB = _boxb;
		box2B[0] = boxB.x*boxB.x;
		box2B[1] = boxB.y*boxB.y;
		box2B[2] = boxB.z*boxB.z;		
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
		System.out.println(" H The size of the box for the B phase is ["+boxB.x+","+boxB.y+","+boxB.z+"].");
		System.out.println(" H The initial density of the B phase is "+calc_B_density());
		System.out.println(" H ");
		// Extra debugging information:
		System.out.println(" H The A diameters are ["+diameter+","+diameter2+"].");
		System.out.println(" H The B/AB diameters are [ "+diameterB+" "+diameterB2+" "+diameterAB2+" ].");
		System.out.println(" H The size of the box for the 0 phase is ["+box[0].x+","+box[0].y+","+box[0].z+"].");
		System.out.println(" H The size of the box for the 1 phase is ["+box[1].x+","+box[1].y+","+box[1].z+"].");
		System.out.println(" H The size of the box for the B phase is ["+boxB.x+","+boxB.y+","+boxB.z+"].");
		System.out.println(" H VolFrac AB2 = "+calc_vol_frac(0));
		System.out.println(" H VolFrac A in AB2 = "+calc_partial_vol_frac(0));
		System.out.println(" H VolFrac B in AB2 = "+calc_partial_vol_frac(1));
		System.out.println(" H VolFrac A = "+calc_vol_frac(1));
		System.out.println(" H VolFrac B = "+calc_vol_frac(2)+" density/cp = "+calc_B_density());
		System.out.println(" H ");
				
	}

	/**
	 * make_nn_lists.
	 *  - This does the main clever stuff.  For the 0 phase, it just runs as normal, but for the 
	 *  1 phase, it only looks for potential neighbours when both particles being tested are &lt; iB (in A) or
	 *  &gt; iB (in B).  This creates two seperate simulation running in the same space, with the same PBCs.
	 *  
	 */
	protected void make_nn_lists() {
		int il,i,j,inn,max_lnn;
		double dr2;

		/* allocate for the 3D NN array*/
		nns = new int[2][n][MAX_NNS];
		
		/* Initialise the neighbour calc: */
		max_dr2 = (nn_cut_off + NN_INC_TOL)*(nn_cut_off + NN_INC_TOL);
		
		/* Loop over both lattices */
		for (il=0; il<2; il++) {
			max_inn[il] = 0;
			/* Loop over all pairs */
			for (i=0; i<n; i++) {
				inn = 0;
				for (j=0; j<n; j++) {
					// This switch only adds neighbours if they are in the same phase - il=0 == AB2, il==1 is 2xFCC:
					if (i!=j && ( il==0 || ( (i<iB && j<iB) || (i>=iB && j>=iB) ) )) {
						/* Calculate pair seperation */
						dr2 = ij_sep2(i,j,il,CUR_POS);
						/* If within the cut-off, add to the list */
						if ( dr2 < max_dr2 ) {
							if( inn >= MAX_NNS ) { 
								System.err.println("make_nn_lists: Tried to add too many neighbours at dr2 = "+dr2);
							} else if( dr2 <= 0.0 ) {
								System.err.println("make_nn_lists: dr2<=0 ("+dr2+") for "+i+","+j+","+latt_name[il]+" "+iB);
							} else {
								nns[il][i][inn] = j;
								inn++;
							}
						}
					}
				}
				if( inn >= MAX_NNS ) { 
					throw( new LSSimException(LSSimException.PARAM, 
							"ERROR: BinaryHSLS.make_nn_lists: Too many neighbours! ["+inn+" >= "+MAX_NNS+"]"));
				}
				nns[il][i][inn] = NO_MORE_NNS;
				if (inn > max_inn[il]) max_inn[il] = inn;
			}
		}
		
		max_lnn = max_inn[0];
		if (max_inn[1] > max_lnn) max_lnn = max_inn[1];
		if (max_lnn > MAX_NNS) {
			throw( new LSSimException(LSSimException.PARAM, 
					"ljls: too many NNs; "+max_lnn+" > MAX_NNS("+MAX_NNS+"); exiting..."));
		}
		
	}

	
	/**
	 * ij_sep2:
	 * calculates the square of the seperation of particles i & j.
	 * uses nearest image pbcs.
	 * Overridden to make it work for two phases, by switching between two boxes
	 * depending on whether we are calculating a seperation in AB2/A or in B.
	 * The B phase is generally a smaller box with smaller particles.
	 * 
	 **/
	protected double ij_sep2(int i, int j, int il, int trial) {
		double[] d = new double[3], idisp = new double[3];
		
		/* use trial or current displacement: */
		if( trial == CUR_POS ) {
			idisp[0] = disp[i].x;
			idisp[1] = disp[i].y;
			idisp[2] = disp[i].z;
		} else if( trial == TRY_POS ) {
			idisp[0] = tdisp.x;
			idisp[1] = tdisp.y;
			idisp[2] = tdisp.z;
		}
				
		/* calculate ij vector */
		if( trial == ZERO_POS ) {
			d[0] = latt[il][j].x - latt[il][i].x;
			d[1] = latt[il][j].y - latt[il][i].y;
			d[2] = latt[il][j].z - latt[il][i].z;
		} else {
			d[0] = (latt[il][j].x + disp[j].x) - (latt[il][i].x + idisp[0]);
			d[1] = (latt[il][j].y + disp[j].y) - (latt[il][i].y + idisp[1]);
			d[2] = (latt[il][j].z + disp[j].z) - (latt[il][i].z + idisp[2]);
		}
		
		/* apply (EXACT) periodic boundary conditions */
		d[0] = d[0] - (int)(d[0]);
		d[0] = d[0] - (int)(d[0]+d[0]);
		d[1] = d[1] - (int)(d[1]);
		d[1] = d[1] - (int)(d[1]+d[1]);
		d[2] = d[2] - (int)(d[2]);
		d[2] = d[2] - (int)(d[2]+d[2]);
		
		/* return the square of the seperation in real space */
		// This changes the box depending on the LATTICE & NUMBER:
	    if( il == 0 || i < iB ) {
			// A-FCC phase and Mixed-phase uses normal boxes:
			return(d[0]*d[0]*box2[il][0] + d[1]*d[1]*box2[il][1] + d[2]*d[2]*box2[il][2]);
		} else {
			// But B-FCC phase uses this box:
			return(d[0]*d[0]*box2B[0] + d[1]*d[1]*box2B[1] + d[2]*d[2]*box2B[2]);			
		}
	}
	
	/**
	 * ij_inter:
	 *   Calculate the overlap for a hard-sphere system, as function of dr^2.
	 *   If index&lt;iB then we have an A-particle, &gt; iB is a B-particle.
	 */
	protected double ij_inter(double dr2, int i, int j, int il ) {
		double testdiam2;

		// If A-A interaction:
		if( i < iB && j < iB ) {
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
	 * Extra calcs to be done when a sphere is moved - used to update CoMs
	 */
	public void accepted_sphere_move(int _im, VectorD3 _dR, double _dE, double _dEM) {
		// Do the common stuff:
		super.accepted_sphere_move(_im, _dR, _dE, _dEM);
		// Update the A and B CoMs:
		/*
		if( _im < iB ) {
			com_A.x += _dR.x/nA;
			com_A.y += _dR.y/nA;
			com_A.z += _dR.z/nA;
		} else {
			com_B.x += _dR.x/nB;			
			com_B.y += _dR.y/nB;			
			com_B.z += _dR.z/nB;			
		}
		*/
	}

	/**
	 * Extra moves sweep.
	 * 
	 * In this case, we must also RW the second volume (the B-crystal).
	 * i.e. Calculate the energy change for a given fluctuation in the volume of 
	 * the B crystal. 
	 * 
	 */
	protected void simulation_sweep_extra( int ip ) {
		// Only implement additional B-Volume moves in an NPT simulation:
		if( NPT_SIM ) {
			// RW the B-volume. (looks rather like the main volume move code.
			double old_m, new_m, old_e, new_e, new_em;
			double nuVol, nuLinear, oldVol, cond;
			double old_ge, new_ge;
			int ss_decision;
			
			/* Only execute this routine one nth of the time */
			if ( rng.raw() >= one_over_n ) return;
			vol_trial+=1.0;
			
			/* store the old aspect ratio */
			VectorD3 old_box = boxB.copy();
			/* Calculate the old scales */
			oldVol = boxB.x*boxB.y*boxB.z;
			/* calc/store old orderparameter/energies */
			old_e = e;
			old_m = (2*c_lat-1)*(em-e);
			old_ge = calc_gs_by_struct(calc_sys_density(), c_lat);
			
			/* create a trial volume move uniform in V */
			nuVol = oldVol + dVol*( (rng.raw()) - 0.5 );
			nuLinear = Math.pow(nuVol/(init_boxB.x*init_boxB.y*init_boxB.z),1.0/3.0);
			setBoxB(init_boxB.times(nuLinear));
			
			/* calculate the order-parameter/energies */
			new_e = calc_e_from_scratch(c_lat);
			new_em = this.calc_e_from_scratch(m_lat);
			
			new_ge = calc_gs_by_struct(calc_sys_density(), c_lat);
			new_m = (2*c_lat-1)*(new_em-new_e);
			
			/* calc the acceptance probability */
			cond = -(beta*((new_ge-old_ge)+(new_e-old_e))+BetaPres*(nuVol-oldVol)-n*Math.log(nuVol/oldVol));
			
			/* accept? - update energies &c */
			ss_decision = ss.control( old_m, new_m, cond );
			if( ss_decision == 1 ) {
				/* accept */
				e = new_e;
				em = new_em;
				/* update acc counter */
				vol_acc+=1.0;
			} else {
				setBoxB( old_box );
			}
		}
	}

	/**
	 * Output extra information pertaining to this binary system.
	 */
	protected void report_simulation(int sweep) {
		// Perform general reporting:
		super.report_simulation(sweep);
		// Report CoM positions:
		/*
		System.out.print(" CoM A "+com_A.x+" "+com_A.y+" "+com_A.z);
		System.out.println(" B "+com_B.x+" "+com_B.y+" "+com_B.z);
		*/
		// Output the current B volume:
		if( NPT_SIM ) {
			System.out.println(" VB "+sweep+" "+(boxB.x*boxB.y*boxB.z)+" "+calc_B_density());
			// Add the current B density to the block analysis.
			bDensblk.add(sweep,calc_B_density());
		}
	}
	
	
	/**
	 * Perform extra checks and output for this system
	 */
	protected void check_simulation(int sweep) {
		// Perform all the usual checks:
		super.check_simulation(sweep);
		// If we are in the double-phase output the two simulated systems into seperate XYZ files.
		if( c_lat == 1 ) { // c_lat == 1 is the double-phase (e.g. 2*FCC).
			save_config_xyz(1, latt[1], box[1], 0, iB,"midsim-"+sweep+"-A-phase");
			save_config_xyz(1, latt[1], boxB, iB, nB,"midsim-"+sweep+"-B-phase");
		} else {
			save_config_xyz(1, latt[1], box[1], 0, iB,"midsim-conj-"+sweep+"-A-phase");
			save_config_xyz(1, latt[1], boxB, iB, nB,"midsim-conj-"+sweep+"-B-phase");
		}
	}

	/**
	 * Output extra information pertaining to this binary system.
	 */
	protected void output_sim_info_tail() {
		// Output information on the B-volume and the AB2 density.
		bDensblk.writeBlockAnalysisGraph("bdens.blk.dat");
		System.out.println(" Bdens "+bDensblk.getMean()+" +/- "+bDensblk.getStdErr());
	}
	
	/**
	 * Calculate the density of the B-phase:
	 */
	private double calc_B_density() {
		return(nBoA*lsize[1][0]*lsize[1][1]*lsize[1][2]*Math.pow(diameterB,3.0)/(Math.sqrt(2.0)*boxB.x*boxB.y*boxB.z));
	}
	
	/**
	 * Calculate the density of the AB and pure A phases:
	 */
	protected double calc_density() {
		System.out.println(" D DEBUG ALERT calc_density on BinaryHSLS!");
		if( c_lat == 0 ) {
			return(lsize[0][0]*lsize[0][1]*lsize[0][2]*(Math.pow(diameter,3.0)+nBoA*Math.pow(diameterB,3.0))/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));
		} else {
			return(lsize[1][0]*lsize[1][1]*lsize[1][2]*(Math.pow(diameter,3.0))/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));			
		}
	}
	
     /**
      * Calculate the volume fraction of the phases
      * 
      * @param phase [0 = AB phase, 1 = A phase, 2 = B phase]
      * @return The volume fraction of the specified phase.
      */
	protected double calc_vol_frac(int phase) {
		if( phase == 0 ) {
			return ( lsize[0][0]*lsize[0][1]*lsize[0][2]*(sphere_volume(diameter) + nBoA*sphere_volume(diameterB)) )/(box[0].x*box[0].y*box[0].z);
		} else if( phase == 1 ) {
			return (lsize[1][0]*lsize[1][1]*lsize[1][2]* sphere_volume(diameter) )/(box[1].x*box[1].y*box[1].z);
		} else if( phase == 2 ) {
			return (nBoA*lsize[1][0]*lsize[1][1]*lsize[1][2]*sphere_volume(diameterB) )/(boxB.x*boxB.y*boxB.z);		
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
			return (lsize[0][0]*lsize[0][1]*lsize[0][2]* sphere_volume(diameter) )/(box[0].x*box[0].y*box[0].z);
		} else if( phase == 1 ) {
			return (nBoA*lsize[1][0]*lsize[1][1]*lsize[1][2]* sphere_volume(diameterB) )/(box[0].x*box[0].y*box[0].z);
		} else {
			// Fail with a runtime error:
			throw( new LSSimException(LSSimException.PARAM, "calc_partial_vol_frac: Unknown phase number"+phase ));
		}
	}

	/**
	 * Calculate the volume of a sphere.
	 * @param _diameter The diameter of the sphere.
	 * @return The volume of that sphere.
	 */
	protected double sphere_volume( double _diameter ) {
		return 4.0*Math.PI*Math.pow(_diameter/2.0,3.0)/3.0;
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
            if( i < iB ) {
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
	
}

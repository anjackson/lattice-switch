/**-------------------------------------------------------------
 * jBinLats - LatticeSwitchSimulation.java
 * net.anjackson.physics.ls.LatticeSwitchSimulation
 * 
 * Created on 16-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.Date;

import edu.cornell.lassp.houle.RngPack.RanMT;
import edu.cornell.lassp.houle.RngPack.RandomSeedable;
import edu.cornell.lassp.houle.RngPack.Ranmar;
import net.anjackson.maths.VectorD3;
import net.anjackson.physics.ls.switches.SwitchMap;
import net.anjackson.physics.mc.StrongSampling;
import net.anjackson.utils.BlockAnalysis;
import net.anjackson.utils.ParamFileReader;
import net.anjackson.utils.Ran250_521;

/**
 * Abstract Lattice-Switch Monte Carlo Code
 * <p>
 * This is the super-class for all the Lattice-Switch Monte Carlo codes.
 * It defines the basic framework, structure and analysis for the 
 * simulations.  It does not specify the interaction or the crystal structures
 * which should be defined in sub-classes.
 * </p><p>
 * The sampling is controlled by the StrongSampling routines, which combine simple
 * importance sampling, MCMC and my own 'Strong' sampling routines into a single 
 * package.
 * </p><p>
 * The interaction determined through nearest-neighbour lists, and so is not 
 * this code is not directly amenable to fluid simulation.  It can fall back to
 * an O(N^2) calculation, but this will not be terribly efficient.
 * </p>
 * 
 * <h3>To Do List</h3>
 * <p>
 * <ul>
 * <li>FIXME Are the last measurements added twice?  Need proper measurement system.</li>
 * <li>FIXME NOTE: the differences in box size/shape make the switch impossible between v different structs (CsCl-FCC, NPT)?  Looks like it.</li>
 * <li>TODO Allow a ratio of volumes when switching between phases in NPT?  Stabilise with conseq for free-energy diff?</li>
 * <li>TODO Implement Ideal-Gas -to- Potential switch, used to refine Mansoori, may require cluster moves.</li>
 * <li>TODO Allow random initial conditions.  Requires Nlist updates periodically.</li>
 * <li>TODO Work out why I get a shoulder distribution 0.70,d0.75 wf-evol.</li>
 * <li>TODO Add test code to analyse the number of NNs and collect by #, output seperations for each.</li>
 * <li>TODO Talk to G and decide on paper strategies/titles.</li>
 * <li>TODO FUTURE Soft-potentials and LJ code need testing.</li>
 * <li>TODO FUTURE Refactor potentials into their own classes, LS just invokes from abstract base class.</li>
 * <li>IDEA Does u-scaling require proper boundary conditions so particles are always within the box? HOW to apply hard PBCs?</li>
 * <li>IDEA Add option to calculate CoMs and output them?</li>
 * <li>IDEA Re-scale everything to use normal unit cells (site-sep=1, not radius=1)?</li>
 * <li>IDEA Allow 'overlappiness' (sum of max{dr2-sep2,0}) as the order parameter?  Would this allow direct boxed-gas to crystal switching?</li>
 * <li>IDEA Use the strictfp (widefp/associativefp) modifiers?</li>
 * <li>IDEA Include the GPL with this software?</li>
 * <li>IDEA Allow radius moves to look at ordering phenomena? See Pronk reference concerning efficiency (p54). This requires an estimated chemical potential function, I think.</li>
 * <li>IDEA Allow particle identity swapping in one or both phases. Useful for studying A/B populations?  This requires an estimated chemical potential, I think.</li>
 * <li>IDEA Create run-time visualization?</li>
 * </p>
 *
 * <h3>Testing</h3>
 * <p>
 * The code has been constructed so that as much of the common algorithm as possible is shared between all the different potentials as lattices.
 * Therefore, by testing one system, we are testing the vast majority of the code for all systems.
 * <ul>
 * <li>Hard-sphere NVT code has been tested quite lot.
 *   <ul>
 *     <li>At N=216, density=0.7778 & using iseed=1140695816572 for 100,000,000 sweeps you get 'F PDFBlock T= 100000000 Df= 0.0010916226917780229 +/- 4.17029288363882E-5'.</li>
 *     <li>TODO At N=216 & density=0.7778, using 1,000,000 sweeps you get '?'.</li>
 *     <li>At N=1728 & density=0.7778, less sweeps are required, ... </li>
 *   </ul></li>
 * <li>Hard-sphere NPT code has not been tested yet.</li>
 * <li>Lennard-Jones code has not been tested yet, and may not be tested soon as non of that code overlaps with the binary hard-sphere crystal code</li>
 * </ul>
 * </p>
 * 
 * <h3>References</h3>
 * <p>
 * <ul>
 * <li>ANJ2001: The primary reference for this work is my thesis (Structural Phase Behaviour Via Monte Carlo Techniques, Andrew N Jackson, 2001)</li>
 * <li>...</li>
 * <li>NBW2002: N B Wilding and P Sollich. Grand canonical ensemble simulation studies of polydisperse fluids. Journal of Chemical Physics, 16:7116-7126, 2002.</li>
 * </ul>
 * </p>
 * 
 * <p>
 * $Id:LJLS.java 451 2005-12-15 16:14:42Z anj $
 * First version: 19th August 1999 : Andrew N Jackson
 * </p>
 * 
 * @author ajackso1
 * @version $Id: LatticeSwitchSimulation.java 1168 2006-11-30 16:58:05Z anj $
 *
 */
public abstract class LatticeSwitchSimulation implements Serializable {
	
	/** DEBUG flag: 0 means no debug info, higher numbers mean more output. */
	protected boolean DEBUG = true;
	
	/** SAVER FLAG: true if the simulation should save itself to a serialized .conf file on check/final */
	protected boolean SAVE_CONFIG = false;
	
	/** The SwitchMap that defines the LatticeSwitch we are simulating */
	protected SwitchMap lsm = null;
	
	/** False = NVT ensemble, True = NPT ensemble */
	protected boolean NPT_SIM = false;
	
	/** Is the NPT box allowed to flex - i.e. variable aspect ratio */
	protected boolean npt_fixed_aspect = true;
	
	/** Lattice-Switch off/on 0/1 */
	protected static boolean LAT_SWITCH = true;
	
	/** Initial phase (depends on sub-class) e.g. hcp/fcc 0/1 */
	protected static int INIT_PHASE = 0;
	
	/** Should interactions be switched off in the zeroth phase? */
	protected static boolean zero_phase_ig = false;
	
	/** Particle move type 1 = RW*/
	protected static int MOVE_RW = 1;
	/** Particle move type  2 = RW + cut-off */
	protected static int MOVE_RWCUT = 2;
	/** Particle move type 3 = Top-Hat */
	protected static int MOVE_TOPHAT = 3;
	/** Basic move default type */
	protected int move_type = MOVE_RW;
	
	/** frequency of basis-flip attempts:
		0 = once every move,
		1 = <once every sweep>,
		2 = as defined be bflipfreq */
	static int FLIP_FREQ = 1;

	/** Probability of a basis-flip/move for FLIP_FREQ=2 */
	protected double bflipfreq = 0.5;

	/** Max number of nearest neighbours to include */
	protected static int MAX_NNS = 400;
	
	/** NN distance-calc inclusion tolerance */
	protected static double NN_INC_TOL = 1e-6;
	
	/** Enable run-time visualization - NOT YET IMPLEMENTED */
	static boolean VISUALIZATION = false;
	
	/** Calc+ouput the CoM vector? */
	protected static int CALC_COM = 0;
	
	/** Calculate the Virial pressure */
	protected static int CALC_VIRP = 0;
	
	/** Automatically tweak the move step sizes */
	protected static int AUTO_STEPS = 1;
	/** Particle move acceptance rate:
	 *  Ideal value with low and high tolerances before tweaking */
	double DrARlo = 0.3;
	double DrARid = 0.35;
	double DrARhi = 0.4;
	double DrARtweak = 0.05;
	/** Volume move acceptance rate:
	 *  Ideal value with low and high tolerances before tweaking */	
	double DvARlo = 0.45;
	double DvARid = 0.5;
	double DvARhi = 0.55;
	double DvARtweak = 0.05;


	/** nnlists: flag used in list of NNs to indicate that there are no more NNs */
	protected static int NO_MORE_NNS = -1;
	/** Have the NNs been calculated? false = no, true = yes */
	protected boolean init_nncalc = false;

	/** verbosity flags: how much output: */
	protected static int LOUDLY = 0;
	static int QUIETLY = 1;
	
	/** ij_sep2 flag: chooses which seperation to calculate */
	protected static final int CUR_POS = 0;
	protected static final int TRY_POS = 1;
	protected static final int ZERO_POS = 2;
	
	/** Turn LS into a neighbour switch 0=LS, 1=N-Switch*/
	protected static int N_SWITCH = 0;
	/** Tolerance of energy checks */
	protected static double E_TOL = 1e-6;
	/** Use O(N^2) calculation instead of NN Lists */
	protected static int MAKEIT_NSQR = 0;

	/** Shall we evolve the weight-function, recalculating it from the TPM when the system is checked?*/
	protected boolean evolve_weights = false;
	/** Shall we reset the displacements to zero when the system is checked?*/
	protected boolean reset_displacements_on_check = false;
	
	/** System size in cells, [lattice][x/y/z-direction] */
	protected int[][] lsize = new int[2][3];
	/** System size - overall system size = npc*nx*ny*nz */
	protected int n;
	/** The value of 1/n (for conveinience, speed-up calc) */
	protected double one_over_n;
	
	/** The pressure (NPT ensemble) */
	protected double Pres;
	/** Beta*Pressure for convenience */
	protected double BetaPres;
	/** The volume-move step size */
	protected double dVol;
	/** density */
	protected double density;
	/** volume fraction */
	protected double volfrac;
	/** The temperature */
	protected double temp;
	/** The inverse temperature */
	protected double beta;
	/** The lattice dilation, c/a ratio, for each lattice.  Default to 1.0 */
	double cOa[] = { 1.0, 1.0 };

	/** The sphere-move step-size */
	protected double dr;
	/** The sphere-move maximum displacement */
	protected double dr_limit = -1;
	/** The sphere-move maximum displacement squared */
	protected double dr_limit2 = -1;
	/** Neighbour interaction cut-off distance */
	protected double nn_cut_off;
	/** Neighbour interaction cut-off distance squared */
	protected double max_dr2;
	/** interaction cut-off distance */
	protected double pot_cut_off;
	/** interaction cut-off distance squared */
	protected double pot_cut_off2;

	/** true = cold start, false = load old conf */
	protected boolean cold_start;
	/** total Monte Carlo Sweeps (MCS) to run */
	protected int total_sweeps;
	/** equilibration period in MCS */
	protected int equib_sweeps;
	/** output period in MCS */
	protected int output_period;
	/** system checking period in MCS */
	protected int check_period;
	/** random number generator seed */
	protected long iseed;
	/** The random number generator */
	protected RandomSeedable rng;
	/** Flag specified to use the Mersenne Twister RNG */
	protected static final int RANMT = 0;
	/** Flag specified to use the RanMar */
	protected static final int RANMAR = 1;
	/** Flag specified to use the 250/521 RNG, as used in the original HS calcs. */
	protected static final int RAN250_521 = 2;
	/** Flag specified to use the RanLux RNG */
	protected static final int RANLUX = 3;
	/** Flag specified to use the RanEcu RNG */
	protected static final int RANECU = 4;
	/** The type of RNG to use: */
	protected int rngType = RANMT;
	
	/** pointer to 2*lattice array */
	protected VectorD3[][] latt = new VectorD3[2][];
	/** textual names for the lattices */
	protected String[] latt_name = { "-", "-"};
	/** pointer to displacements array */
	protected VectorD3[] disp;
	/** change-in-displacement (DeltaR) holder */
	protected VectorD3 Ddisp = new VectorD3();
	/** trial displacement holder */
	protected VectorD3 tdisp = new VectorD3();
	/** nearest neighbour table */
	protected int[][][] nns;
	
	/** Energy in the currect phase */
	protected double e;
	/** Ground-state energy difference*/
	protected double dEgs;
	/** conjugate-phase energy */
	protected double em;
	/** the order parameter */
	protected double m;
	/* Accumilators */
	/** To calc average energy */
	protected double num_e;
	protected double ave_e;
	/** To calc average energy change/move */
	protected double num_De;
	protected double ave_De;
	/** To calc average energy cost of switch */
	protected double num_es;
	protected double ave_es;
	/* Self-check variables */
	/** The current system energy, to compare with running totals */
	double e_check;
	/** The conjugate system energy, to compare with running totals */
	double em_check;


	/** identity of the (c)urrent lattice */
	protected int c_lat;
	/** identity of the conjugate (m) order-param lattice */
	protected int m_lat;
	/** the MCMC weight function */
	double[] weights;
	/** range of the MCMC weight func in m[wlo_m,whi_m] */
	int wlo, whi;
	
	/** inital size of system box a,b,c */
	protected VectorD3[] init_boxes = { new VectorD3(), new VectorD3() };
	/** lattice box a,b,c */
	protected VectorD3[] box = { new VectorD3(), new VectorD3() };
	/** lattice box squared: a^2,b^2,c^2 */
	protected double[][] box2 = new double[2][3];
	/** the old box-size, for trial moves */
	protected double[] obox = new double[3];
	/** NPT box-size RW step for each lattice, in each dirn */
	protected double[][] drobox = new double[2][3];
	/** sphere diameter */
	protected double diameter;
	/** sphere diameter^2 */
	protected double diameter2;
	/** Ratio of the volfraction to the density wrt close packing.*/
	protected double VolFracOverDensity = Math.sqrt(2.0)*Math.PI/6.0;
	
	/** Polydispersity calculations - enabled? */
	protected boolean polyd = false;
	/** Diameters for the two phases */
	protected double polyd_d[] = { 1.0, 1.0 };
	/** Polydispersity as percentages of diameters for the two phases */
	protected double polyd_Dd[] = { 0.0, 0.0 };
	/** The diameter array, indexed by lattice and particle number */
	protected double pdiameter[][] = null;
	/** The diameter2 array, indexed by lattice and particle number */
	protected double pdiameter2[][] = null;
	
	/** the observed max no. of NNs */
	protected int[] max_inn = new int[2];
	/** The time spent in each phase */
	protected double[] time_in_phase = new double[2];
	/** The time spent in each phase, for block analysis */
	protected double[] btime_in_phase = new double[2];
	/** stores # of trial sphere moves */
	protected double trial_moves;
	/** stores # of accepted sphere moves */
	protected double acc_moves;
	/** the number of trial volume moves */
	protected double vol_trial;
	/** the number of accepted volume moves */
	protected double vol_acc;
	
	/** Local-displacement scaling factors - for mapping scaled displacements between lattices */
	protected double[] dispScaling = new double[] { 1.0, 1.0 };
	
	/* --- Block Analysis Members - Doubled up, one for each phase: --- */
	
	/** Block analysis of the energy */
	protected BlockAnalysis[] bE = { new BlockAnalysis(), new BlockAnalysis() };
	/** Block analysis of the Order Parameter */
	protected BlockAnalysis[] bM = { new BlockAnalysis(), new BlockAnalysis() };
	/** Block analysis of the density */
	protected BlockAnalysis[] bD = { new BlockAnalysis(), new BlockAnalysis() };
	/** Block analysis of the c/a ratio */
	protected BlockAnalysis[] bCoA = { new BlockAnalysis(), new BlockAnalysis() };
	/** Block analysis of the time-spent-in-phase */
	protected BlockAnalysis[] bTIP = { new BlockAnalysis(), new BlockAnalysis() };
	/** Block analysis of the PDF weight of phase */
	protected BlockAnalysis[] bPDFw = { new BlockAnalysis(), new BlockAnalysis() };	
	/** Block analysis of the free-energy difference */
	protected BlockAnalysis bDf = new BlockAnalysis();
	/** Block analysis of the free-energy difference calculated from visit time */
	protected BlockAnalysis bDfV = new BlockAnalysis();
	
	/* --- StrongSampling Parameters --- */
	
	/** The Strong Sampling MC controller */
	protected StrongSampling ss = new StrongSampling();	
	/** The sampling algorithm: strong, mcmc or imp */
	protected int samp_type;
	/** The order-parameter (OP) histogram: lower bound */
	protected double ss_lo;
	/** The OP histogram: upper bound */
	protected double ss_hi;
	/** clip weight function automatically */
	protected boolean wf_clip_auto = true;
	/** clip weight function up to here */
	protected double wf_clip_lo;
	/** clip weight function from here up */
	protected double wf_clip_hi;
	/** OP histogram: no. of bins */
	protected int ss_n;
	/** OP Transition Probability Matrix (TPM): max. +/- matrix diagonal size */
	protected int ss_nd;
	/** ss: period (in moves) for macrolock */
	protected int ss_lock;
	/** Use the TPM: true) empty the tpm, 0) use the old tpm */
	protected boolean samp_new_tpm;
	/** Use weight-function: true) mcmc uses "weights.in", 0) doesn't */
	protected boolean read_weights;
	/** exponential shift for no overflows, used in StrongSampling */
	static final double EXPO_SHIFT = 30.0;

	/** Filenames: The pre-simulation filename */
	protected String file_initial = "init";
	/** Filenames: The mid-simulation filename */
	protected String file_midsim = "midsim";
	/** Filenames: The final end-of-simulation filename */
	protected String file_final = "final";
	/** Filenames: The TPM data input file */
	protected String file_tpm_in = "TPM.in";
	/** Filenames: The TPM data output file */
	protected String file_tpm_out = "TPM.dat";
	/** Filenames: The PDF output file (via visited-states) */
	protected String file_pdf_vs = "pdf.vs.dat";
	/** Filenames: The PDF output file (via TPM analysis) */
	protected String file_pdf_tpm = "pdf.tpm.dat";
	/** Filenames: The weight function to be read in: */
	protected String file_weights_in = "weights.in";
	/** Filenames: The weight function used during this simulation: */
	protected String file_mcmc_weights_out = "mcmc.weights.dat";
	/** Filenames: The nearest-neighbour arrays */
	protected String file_nns = "nns";
	/** Filenames: The initial configuration to load (not a cold start): */
	protected String file_conf_in = "init.conf";
	
	
	/*------------------------------ METHODS -----------------------------*/
	
	/*----------------------
	 Simulation subroutines
	 ----------------------*/
	
	/**
	 * Setter to define the switch-map to be used.
	 * 
	 * @param _lsm The SwitchMap to simulate upon.
	 */
	public void setSwitchMap( SwitchMap _lsm ) {
		lsm = _lsm;
	}
	
	/**
	 * Load in the parameters for this simulation.
	 * 
	 * @param pars the ParameterFileReader to read from.
	 */
	protected void read_params_from_user( ParamFileReader pars ) {
		/* The diameter defaults to 1.0 */
		setDiameter(1.0);		
		
		// Now try to look up all the parameters: Required parameters first:
		try {
			// System size:
            if( pars.isDefined("l0.nx")) {
            	// We are specifying both phases seperately:
            	lsize[0][0] = pars.getInteger("l0.nx");
            	lsize[0][1] = pars.getInteger("l0.ny");
            	lsize[0][2] = pars.getInteger("l0.nz");
            	lsize[1][0] = pars.getInteger("l1.nx");
            	lsize[1][1] = pars.getInteger("l1.ny");
            	lsize[1][2] = pars.getInteger("l1.nz");
            } else {
            	// Determine the system size:
            	lsize[0][0] = pars.getInteger("n.x");
            	lsize[0][1] = pars.getInteger("n.y");
            	lsize[0][2] = pars.getInteger("n.z");
            	lsize[1][0] = pars.getInteger("n.x");
            	lsize[1][1] = pars.getInteger("n.y");
            	lsize[1][2] = pars.getInteger("n.z");
            }

			// General system parameters
			if( !pars.isDefined("polyd")) {
				if( pars.isDefined("density")) {
					setDensity(pars.getDouble("density"));
					if( pars.isDefined("volfrac ")) 
						throw( new LSSimException(LSSimException.PARAM, "You should only specify one of density and volfrac, not both of them!"));
				} else if( pars.isDefined("volfrac")) {
					setVolFrac(pars.getDouble("volfrac"));
				} else {
					throw( new LSSimException(LSSimException.PARAM,"You must specify either a density OR a volfrac!"));
				}
			} else {
				// Allow polydisperse calculations:
				polyd = pars.getBoolean("polyd");
				if( polyd ) setPolydispersity(
						pars.getDouble("polyd.d0"),
						pars.getDouble("polyd.d0.percent"),
						pars.getDouble("polyd.d1"),
						pars.getDouble("polyd.d1.percent"));
			}
			setTemperature(pars.getDouble("temperature"));
			// Get the move type, defaulting to RW, and get any required parameters:
			if( pars.isDefined("particle.move.type")) {
				if( "rw".equalsIgnoreCase(pars.getString("particle.move.type"))) {
					move_type = MOVE_RW;
				} else if( "rw-c".equalsIgnoreCase(pars.getString("particle.move.type"))) {
					move_type = MOVE_RWCUT;
				} else if( "tophat".equalsIgnoreCase(pars.getString("particle.move.type"))) {
					move_type = MOVE_TOPHAT;
				}
			}
			if( move_type == MOVE_RW || move_type == MOVE_RWCUT ) {
				setRWdr(pars.getDouble("dr"));				
			}
			if( move_type == MOVE_TOPHAT || move_type == MOVE_RWCUT ) {
				dr_limit = pars.getDouble("dr.max");
			}
			// Get NPT parameters:
			if( NPT_SIM ) {
				setPressure(pars.getDouble("pressure"));
				dVol = pars.getDouble("dV");
			}
			nn_cut_off = pars.getDouble("nn.incdist");			
			cold_start = ! pars.getBoolean("loadconf");
			// Allow the rng type and the seed to be specified.
			if( pars.isDefined("rng")) {
				String rngid = pars.getString("rng");
				if( "ran250_521".equalsIgnoreCase(rngid)) {
					rngType = RAN250_521;
					// IDEA Issue a warning about this RNG?
				} else if( "ranmt".equalsIgnoreCase(rngid)) {
					rngType = RANMT;
				} else if( "ranmar".equalsIgnoreCase(rngid)) {
					rngType = RANMAR;
				} else if( "ranlux".equalsIgnoreCase(rngid)) {
					rngType = RANLUX;
				} else if( "ranecu".equalsIgnoreCase(rngid)) {
					rngType = RANECU;
				} else {
					throw( new LSSimException( LSSimException.PARAM, "Unrecognised RNG type '"+rngid+"'!"));
				}
			}
			iseed = pars.getLong("rng.seed");

			// Phase parameters
			if( pars.isDefined("latticeswitch.enabled")) {
				LAT_SWITCH = pars.getBoolean("latticeswitch.enabled");
			}
			if( pars.isDefined("init.phase")) {
				INIT_PHASE = pars.getInteger("init.phase");
				if( INIT_PHASE < 0 || INIT_PHASE > 1 ) INIT_PHASE = 0;
			}
			
			// Time parameters
			total_sweeps = pars.getInteger("t.tot");
			equib_sweeps = pars.getInteger("t.equib");
			output_period = pars.getInteger("t.output");
			check_period = pars.getInteger("t.check");
			
			// Strong-sampling parameters:
			String sampmethod = pars.getString("samp.method");
			if( "strong".equals(sampmethod) ) {
				samp_type = StrongSampling.STRONG_SAMP;
			} else if( "mcmc".equals(sampmethod) ) {
				samp_type = StrongSampling.MCMC_SAMP;
			} else if( "imp".equals(sampmethod) ) {
				samp_type = StrongSampling.IMP_SAMP;
			} else {
				System.err.println("Could not understand sampling type "+sampmethod);
			}
			ss_lo = pars.getDouble("samp.histogram.low");
			ss_hi = pars.getDouble("samp.histogram.high");
			ss_n = pars.getInteger("samp.histogram.bins");
			ss_nd = pars.getInteger("samp.histogram.diagonal");
			ss_lock = pars.getInteger("samp.histogram.lock");
			samp_new_tpm = pars.getBoolean("samp.histogram.new_tpm");
			read_weights = pars.getBoolean("samp.weights.read_weights_file");
			// Optional WF clipping:
			if( pars.isDefined("samp.weights.clip.auto")) {
				wf_clip_auto = pars.getBoolean("samp.weights.clip.auto");
			}
			if( pars.isDefined("samp.weights.clip.low")) {
				wf_clip_lo = pars.getDouble("samp.weights.clip.low");
			} else {
				wf_clip_lo = ss_lo;
			}
			if( pars.isDefined("samp.weights.clip.high")) {
				wf_clip_hi = pars.getDouble("samp.weights.clip.high");
			} else {
				wf_clip_hi = ss_hi;
			}
			
			// Now look up optional parameters:
			
			// Allow the c/over/a ratio to be specified:
			if( pars.isDefined("cOa") ) {
				cOa[0] = pars.getDouble("cOa");
				cOa[1] = pars.getDouble("cOa");
			}
			if( pars.isDefined("l0.cOa") )
				cOa[0] = pars.getDouble("l0.cOa");
			if( pars.isDefined("l1.cOa") )
				cOa[1] = pars.getDouble("l1.cOa");			
			
			// Allow the displacements to be scaled differently for each phase:
			if( pars.isDefined("disp.scaling.factor")) 
				dispScaling[1] = dispScaling[0]*pars.getDouble("disp.scaling.factor");
			
			// Allow the weights to be evolved - recalculated every t.check sweeps
			if( pars.isDefined("evolve.weights")  )
				evolve_weights = pars.getBoolean("evolve.weights");
			// Force the displacements to be reset to zero every t.check sweeps
			if( pars.isDefined("reset.displacements.on.check")) 
				reset_displacements_on_check = pars.getBoolean("reset.displacements.on.check");
			
			// Allow fixed or variable aspect ratio
			if( pars.isDefined("npt.fixed.aspect"))
				npt_fixed_aspect = pars.getBoolean("npt.fixed.aspect");
			
			// Allow interactions to be switched off in the zero phase:
			if( pars.isDefined("zero.phase.ig"))
				zero_phase_ig = pars.getBoolean("zero.phase.ig");
			
		} catch( ParamFileReader.UndefinedParameterException e ) {
			// If any of them is not defined, cry and grumble.
			throw(new LSSimException( LSSimException.PARAM, e ) );
		}

	}
	
	/** 
	 * init_variables - initializes various variables and stuff, based 
	 * on the parameters and the SwitchMap.
	 */
	protected void init_variables() {
		// Get the lattices and put them in the latt field:
		for( int i = 0; i < 2; i++ ) {
			latt[i] = lsm.getLatticeVectors(i);
			latt_name[i] = lsm.getLatticeName(i);
			init_boxes[i] = lsm.getBox(i);
		}
		// Record the number of spheres:
		n = latt[0].length;
		
		// Also get the periodic box size:
		setBoxes(init_boxes);
		
		// Allocate space for the displacement array:
		disp = new VectorD3[n];
		
		// And initialise the displacements:
		for( int i = 0; i < n; i++ )
			disp[i] = new VectorD3(0,0,0);
		
		/* check that the NN cut-off is actually inside the box */
		for( int il = 0; il < 2; il++ ) {
			for( int ix = 0; ix < 3; ix++ ) {
				if (Math.sqrt(box2[il][ix]) < nn_cut_off) {
					throw( new LSSimException(LSSimException.PARAM,
							"ERROR: This cut-off ("+nn_cut_off+") does not fit in box ("+ix+","+Math.sqrt(box2[il][ix])+")"));
				}
			}
		}
	
		/** Construct the NN lists */
		make_nn_lists();
		
		// Set the box and derived parameters up:
		initBoxes();
		
		// If we are working in the NVT ensemble, check the densities are equal:
		if( !NPT_SIM ) {
			// Calculate the initial densities and check they are equal - only valid in NVT:
			int tmp_lat = c_lat, dens_errors = 0;
			double denses[] = new double[3];
			c_lat = 0;
			denses[0] = this.calc_density();
			c_lat = 1;
			denses[1] = this.calc_density();
			c_lat = tmp_lat;
			if( Math.abs(denses[0]-denses[1]) > E_TOL ) {
				System.out.println(" ERROR: Densities in 0 and 1 phase are different! 0@"+denses[0]+" q@"+denses[1]);
				dens_errors++;
			}
			if( Math.abs(denses[0]-density) > E_TOL ) {
				System.out.println(" ERROR: Densities in 0 does not match requested density! 0@"+denses[0]+" density@"+density);
				dens_errors++;
			}
			if( dens_errors > 0 ) {
				System.exit(1);
			} else {
				System.out.println(" D CHECK: All densities checked as equal: 0@"+denses[0]+" 1@"+denses[1]+" 2@"+denses[2]);
			}
		}
		
		
		
		/* initial value of m for 'cold' start */
		m = 0;
		
		/* initialise time spent in each phase */
		time_in_phase[0] = 0; time_in_phase[1] = 0;
		btime_in_phase[0] = 0; btime_in_phase[1] = 0;
		
		/* counters */
		trial_moves = 0.0; acc_moves = 0.0;
		num_e = 0.0; ave_e = 0.0;
		num_De = 0.0; ave_De = 0.0;
		num_es = 0.0; ave_es = 0.0;
		vol_acc = 0.0; vol_trial = 0.0;
		
		/* various */
		one_over_n = 1.0/((double)n);
		
		// initialise the random number generator:
		if( iseed <= 0 ) {
			// Use a clock-based seed if non is supplied:
			iseed = RandomSeedable.ClockSeed();
		}
		// Create the RNG:
		if( rngType == RANMT ) {
			rng = new RanMT(iseed);
		} else if( rngType == RANMAR ) {
			rng = new Ranmar(iseed);
		} else if( rngType == RAN250_521 ) {
			System.out.println(" H WARNING!  Using old Ran250/521 RNG!");
			rng = new Ran250_521(iseed);
		} else {
			throw( new LSSimException( LSSimException.PARAM, "Unknown RNG algorithm!"));
		}
		
		/* If this is a polydisperse calculation, build up the arrays */
		if( polyd ) {
			// Allocate some space:
			pdiameter = new double[2][n];
			pdiameter2 = new double[2][n];
			for( int il = 0; il < 2; il++ ) {
				for( int i = 0; i < n; i++ ) {
					// Use a Gaussian of stddev %tage D
					pdiameter[il][i] = polyd_d[il] + rng.gaussian(polyd_d[il]*polyd_Dd[il]/100.0);
					pdiameter2[il][i] = pdiameter[il][i]*pdiameter[il][i];
				}
			}
		}
		// Now recalc density and volfrac.
		System.out.println(" D density= "+density+" volfrac= "+volfrac+" sysdens= "+calc_sys_density() );
		// NOTE that this makes no sense for the binary system, where diameterB is not yet set.
		//density = calc_sys_density();
		volfrac = density*VolFracOverDensity;
		System.out.println(" D density= "+density+" volfrac= "+volfrac+" sysdens= "+calc_sys_density() );
		
		// Initialise the square of the RW cut-off/tophat.
		if( move_type == MOVE_RWCUT || move_type == MOVE_TOPHAT ) {
			dr_limit2 = dr_limit*dr_limit;
		}

	}

	/**
	 * Box dimensions are initialized here, so that they can be overridden if necessary.
	 *
	 */
	protected void initBoxes() {
		// Apply the required c/a ratio, usually 1.0
		for( int il = 0; il < 2; il++ ) {
		  if( cOa[il] > 0.0 ) {
			// Fix the z/y ratio and rescale to fix the volume:
			double o_V = box[il].x*box[il].y*box[il].z;
			box[il].z = cOa[il] * box[il].y;
			double n_V = box[il].x*box[il].y*box[il].z;
		    box[il] = box[il].times( Math.pow(o_V/n_V, 1.0/3.0) );
			/*
			// TODO Determine ratio w.r.t. the current aspect ratio:
			double cOa_r = cOa[il]*(box[il].z/box[il].y);
			box[il].z *= cOa_r;
			box[il].y /= cOa_r;
			box[il].x /= cOa_r;
			*/
		  }
		}
		// Set the box and derived parameters up:
		setBoxes(box);
		
		// If we are working in the NPT ensemble, set both boxes to be the same:
		if( NPT_SIM ) {
			setBoxes(box[INIT_PHASE]);
		}
	}
	
	/**
	 * Related variables setter - set temperature and beta
	 * @param _t The desired dimensionless temperature
	 */
	void setTemperature( double _t ) {
		temp = _t;
		beta = 1/_t;
	}

	/**
	 * Related variables setter - allow the volume fraction and set the density:
	 * @param _volfrac The desired volume fraction:
	 */
    void setVolFrac( double _volfrac ) {
		volfrac = _volfrac;
		setDensity(volfrac/VolFracOverDensity);
	}

    /**
     * Related variables setter - given the density, specify related properties.
     * @param _density The desired density.
     */
	void setDensity( double _density ) {
		/** Whether to fix the density by shrinking the radius, or scaling the box.  Should be equivelant. */
		final boolean density_by_radius = true;
		
		// Update the density & volume fraction:
		density = _density;
		volfrac = density*VolFracOverDensity;
		
		// Fix the diameters:
		if( density_by_radius ) {
			setDiameter(Math.pow(density,1.0/3.0));
		} else {
			setDiameter(1.0);
		}
		
		/* Rescale the simulation box to get the desired density: */
		if( !density_by_radius ) {
			VectorD3[] new_boxes = { new VectorD3(), new VectorD3() };
			for( int il = 0; il < 2; il++ ) {
				new_boxes[il].x = Math.pow(density,-1.0/3.0)*init_boxes[il].x;
				new_boxes[il].y = Math.pow(density,-1.0/3.0)*init_boxes[il].y;
				new_boxes[il].z = Math.pow(density,-1.0/3.0)*init_boxes[il].z;
			}
			setBoxes( new_boxes );
		}
		
		System.out.println(" D Set density to "+density+" c.f. sysdens "+calc_sys_density());
		System.out.println(" D Set volfrac to "+volfrac);		
	}

	/**
	 * Related variables setter - specify the pressure.
	 * @param _pressure The desired dimensionless pressure.
	 */
	void setPressure( double _pressure ) {
		Pres = _pressure/Math.pow(diameter,3.0);
		BetaPres = beta*Pres;
	}

	/**
	 * Related variables setter - specify the same box size for both boxes:
	 * @param _box The desired box-size as a VectorD3.
	 */
	void setBox( VectorD3 _boxes) {
		setBoxes( new VectorD3[]{ _boxes, _boxes });
	}
	
	/**
	 * Related variables setter - specify the box sizes:
	 * @param _box The desired box-sizes as a VectorD3 array.
	 */
	void setBoxes( VectorD3[] _boxes) {
		/* store box side length ^ 2 */
		for( int il = 0; il < 2; il ++ ) {
			box[il] = _boxes[il];
			box2[il][0] = box[il].x*box[il].x;
			box2[il][1] = box[il].y*box[il].y;
			box2[il][2] = box[il].z*box[il].z;
		}
		
		/* This also affect the RW move size (indirectly) - reset using current value */
		setRWdr(dr);
	}
	
	/**
	 * Related variables setter - specify the box sizes using one box:
	 * @param _box The desired box-size as a VectorD3 array.
	 */
	void setBoxes( VectorD3 _box) {
		setBoxes( new VectorD3[]{_box, _box} );
	}
	
	/**
	 * Related variables setter - specify the diameter and the diameter2
	 * @param _diameter The desired diameter.
	 */
	void setDiameter ( double _diameter ) {
		diameter = _diameter;
		diameter2 = diameter*diameter;
	}
	
	/**
	 * Set up a polydisperse calculation, along with associated vars:
	 * @param d0 The diameter for the 0 phase.
	 * @param poly0 The polydispersity for the 0 phase (as %tage of d0).
	 * @param d1 The diameter for the 1 phase.
	 * @param poly1 The polydispersity for the 1 phase (as %tage of d1).
	 */
	void setPolydispersity(double d0, double poly0, double d1, double poly1) {
		// Identify the simulation as polydisperse:
		polyd = true;
		// Record the base diameters and polydispersities
		polyd_d[0] = d0;
		polyd_Dd[0] = poly0;
		polyd_d[1] = d1;
		polyd_Dd[1] = poly1;
	}

	/**
	 * Related variables setter - Set the Random Walk step size.
	 * @param _dr The desired random-walk step size.
	 */
	void setRWdr( double _dr ) {
		dr = _dr;
		/* calc RW step in each dirn */
		for( int il=0; il<2; il++ ) {
			drobox[il][0] = dr/box[il].x;
			drobox[il][1] = dr/box[il].y;
			drobox[il][2] = dr/box[il].z;
		}
	}
	
	/** ---------------------------------------------------------
	 read_weights_file
	 ---------------------------------------------------------*/	
	void read_weights_file(String filename) {
		int wi;
		
		/* Define size etc of the weight array based on n */
		if ( n == 216 ) {
			wlo =  0;
			whi =  200 - 1;
		} else if ( n == 1728 ) {
			wlo =  0;
			whi =  1000 - 1;
		} else if ( n == 5832 ) {
			wlo =  0;
			whi =  10000 - 1;
		} else {
			System.err.print("ERROR: read_weights_file cannot cope with n = "+n);
		}
		/* Calculate where M=0 lies in the [0,whi] array*/
		//wzero = whi/2; // This is not used anywhere in the code.
		
		/* Allocate memory for the weights array */
		weights = new double[whi-wlo+1];
		
		/* If there is a weights file, load it */
		try {
			// The input file:
			BufferedReader in = new BufferedReader(new FileReader(filename));
			for ( wi = wlo; wi <= whi; wi++ ) {
				weights[wi] = Double.parseDouble(in.readLine());
			}
			/* close the weights file */
		} catch( Exception e ) {
			System.out.println("Could not load weights file '"+filename+"' - using flat weight function.");
			/* If no weights file, define a flat weight function */
			for ( wi = wlo; wi <= whi; wi++ ) {
				weights[wi] = 1.0;
			}
		}
		
	}

	/**
	 * Updates the weight function based on the current SS data.
	 * Optionally loads them in from a file.
	 * Also performs clipping if required.
	 */
	public void updateWeightFunction(boolean load_weights) {
		ss.getMCMCWeightsFromTPM();
		if ( load_weights ) ss.loadMCMCWeights(file_weights_in);
		if( wf_clip_auto ) {
			ss.clipMCMCWeights();				
		} else {
			ss.clipMCMCWeights(wf_clip_lo,wf_clip_hi);				
		}
	}
	
	/** -------------------------------------------------------
	 tweak_for_acc_rates
	 Changes the step-sizes to make the acc_rates right.
	 ------------------------------------------------------- */
	protected void tweak_for_acc_rates() {
		double vAR, rAR;
		
		/* If position acc_rate is out of bounds, tweak it */
		rAR = ((double)acc_moves)/((double)trial_moves);
		if( rAR < DrARlo || rAR > DrARhi ) {
			if( rAR > 0.0 && rAR < 1.0 ) {
				setRWdr( dr*Math.log(DrARid)/Math.log(rAR) );
			} else {
				if( rAR < DrARlo ) setRWdr(dr*(1.0 + DrARtweak));
				if( rAR > DrARhi ) setRWdr(dr*(1.0 - DrARtweak));
			}
		}
		
		if( NPT_SIM ) {
			/* If volume acc_rate is out of bounds, tweak it */
			vAR = ((double)vol_acc)/((double)vol_trial);
			if( vAR < DvARlo || vAR > DvARhi ) {
				if( vAR > 0.0 && vAR < 1.0 ) {
					dVol = dVol*Math.log(DvARid)/Math.log(vAR);
				} else {
					if( vAR < DvARlo ) dVol = dVol*(1.0 + DvARtweak);
					if( vAR > DvARhi ) dVol = dVol*(1.0 - DvARtweak);
				}
			}
			/* Notify */
			System.out.println(" AR dr = "+dr+" dVol = "+dVol);
		} else {
			/* Notify */
			System.out.println(" AR dr = "+dr);
		}
		
	}
	
	/**
	 * TODO Describe this new method on LatticeSwitchSimulation.	
	 ---------------------------------------------------------
	 save_config
	 Save configuration to a file
	 - in wrl (vrml) format for use with VRMLview and others.
	 - in xyz format for use with xmol
	 ---------------------------------------------------------
	 *
	 * @param il lattice identity
	 * @param lattice
	 * @param _box
	 * @param start
	 * @param number
	 * @param filename
	 */
	protected void save_config_xyz(int il, VectorD3[] lattice, VectorD3 _box, int start, int number, String filename) {
		int i,ii,ij,ik;
		boolean output_xyz = false;

		/* Output as an XYZ file: */
		if( output_xyz ) {
		/* open file and outpout number of atoms to come */
		String xyzfn = filename+".xyz";
		try {
			PrintStream out = new PrintStream(new FileOutputStream(xyzfn));
			out.print(" "+(number+8)+"\n\n");
			
			/* loop over all particles and output as carbons */
			for ( i = start; i < start+number; i++ ) {
				out.println("C "+(lattice[i].x+disp[i].x)*_box.x+" "+(lattice[i].y+disp[i].y)*_box.y+" "+(lattice[i].z+disp[i].z)*_box.z);
			}
			
			/* Loop over the corners of the pbc box, outputting as hydrogens */
			for ( ii = -1; ii <= 1; ii = ii + 2 ) {
				for ( ij = -1; ij <= 1; ij = ij + 2 ) {
					for ( ik = -1; ik <= 1; ik = ik + 2 ) {
						out.println("X "+_box.x*(ii*0.5)+" "+_box.y*(ij*0.5)+" "+_box.z*(ik*0.5));
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
                out.print("        diffuseColor 1.0 1.0 1.0\n");
	            out.print("          transparency 0.0\n");
	            out.print("        }\n");
	            out.print("      }\n");
	            out.print("      geometry Sphere {\n");
	            out.print("        radius "+0.5*diameter+"\n");
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
	 * Save the current configuration to a file:
	 * 
	 * @param filename The filename (leafname) to use.
	 */
	protected void save_config_xyz(String filename) {
		save_config_xyz(c_lat, latt[c_lat], box[c_lat], 0, n, filename);
	}
	
	/** 
	 * The simulation loader, actually a simulation class factory based
	 * on standard java object serialization.
	 * 
	 * See http://java.sun.com/docs/books/tutorial/essential/io/providing.html
	 * and http://java.sun.com/j2se/1.4.2/docs/api/java/io/ObjectInputStream.html
	 *  
	 */
	public static LatticeSwitchSimulation load_config(String filename) {
		/* The simulation to return */
		LatticeSwitchSimulation ls_sim = null;
		if( ls_sim.DEBUG ) 
			System.out.println("load_config: Attempting to load LS simulation class from file "+filename);
		
		// Attempt to read the simulation state from a file */
		try {
			// Open the object serialization output stream:
			FileInputStream fis = new FileInputStream(filename);
			ObjectInputStream ois = new ObjectInputStream(fis);

			// Read in the date and the object (the entire simulation) */
			Date tstamp = (Date) ois.readObject();
			Object sim_obj = ois.readObject();
			
			// Notify the user:
			if( ls_sim.DEBUG ) {
				System.out.println("load_config: Loader factory loaded a simulation class '"
					+ sim_obj.getClass().getName() + "' dated " + tstamp);
			}
			
			// Cast the object as a generic LS simulation:
			ls_sim = (LatticeSwitchSimulation)sim_obj; 

			// Close the input stream:
			ois.close();
			
		} catch( FileNotFoundException e ) {
			System.err.println("save_config: FileNotFoundException while saving this class! "+e);
		} catch( IOException e ) {
			System.err.println("save_config: IOException while saving this class! "+e);			
		} catch( ClassNotFoundException e ) {
			
		}
		
		// Pass the new simulation object back to the user:
		return ls_sim;
		
	}

	/**
	 * Saves the entire simulation state via Java object serialization.
	 * See http://java.sun.com/docs/books/tutorial/essential/io/providing.html
	 * and http://java.sun.com/j2se/1.4.2/docs/api/java/io/ObjectOutputStream.html
	 *  
	 */
	protected void save_config(String filename) {
		if( DEBUG ) 
			System.out.println("save_config: Attempting to save this LS simulation class into a file called "+filename);
		
		// Attempt to save this simulation as a file:
		try {
			// Open the object serialization output stream:
			FileOutputStream fos = new FileOutputStream(filename);
			ObjectOutputStream oos = new ObjectOutputStream(fos);

			// Write out the date and this object (the entire simulation) */
			Date tstamp = new Date();
			oos.writeObject(tstamp);
			oos.writeObject(this);
			
			// Notify the user:
			if( DEBUG ) {
				System.out.println("save_config: saved this simulation class into file '"
					+ filename + "' including date " + tstamp);
			}
			
			// Close the output stream:
			oos.close();
			
		} catch( FileNotFoundException e ) {
			System.err.println("save_config: FileNotFoundException while saving this class! "+e);
		} catch( IOException e ) {
			System.err.println("save_config: IOException while saving this class! "+e);			
		}
	}

	/** -------------------------------------------------------------
	 ij_sep2:
	 - calculates the square of the seperation of particles i & j.
	 - uses nearest image pbcs.
	 -------------------------------------------------------------*/
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
			d[0] = (latt[il][j].x + dispScaling[il]*disp[j].x) - (latt[il][i].x + dispScaling[il]*idisp[0]);
			d[1] = (latt[il][j].y + dispScaling[il]*disp[j].y) - (latt[il][i].y + dispScaling[il]*idisp[1]);
			d[2] = (latt[il][j].z + dispScaling[il]*disp[j].z) - (latt[il][i].z + dispScaling[il]*idisp[2]);
		}
		
		/* apply (EXACT) periodic boundary conditions, via truncate towards zero (int) */
		for( int di = 0; di < 3; di++ ) {
			d[di] = d[di] - (double)((int)(d[di]));
			d[di] = d[di] - (double)((int)(d[di]+d[di]));
		}
		/* return the square of the seperation in real space */
		return(d[0]*d[0]*box2[il][0] + d[1]*d[1]*box2[il][1] + d[2]*d[2]*box2[il][2]);
		
	}

	/** ---------------------------------------------------------
	 make_nn_lists:
	 - builds the static neighbour lists:
	 ---------------------------------------------------------*/
	protected void make_nn_lists() {
		int il,i,j,inn,max_lnn;
		double dr2;
		
		/* allocate for the 3D NN array*/
		//nns = imatrix3d(0,1,0,n,0,MAX_NNS);
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
					if (i!=j) {
						/* Calculate pair seperation */
						dr2 = ij_sep2(i,j,il,CUR_POS);
						/*
						 if ( il == 0 && i == 0 ) System.out.printf(" NN %i %i %i %lf %lf\n",il,i,j,dr2,sqrt(dr2));
						 */
						/* If within the cut-off, add to the list */
						if ( dr2 < max_dr2 ) {
							nns[il][i][inn] = j;
							inn++;
						}
					}
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

	/** -------------------------------------------------------
	 - Writes out NN list
	 ------------------------------------------------------- */
	protected void write_out_nn_list(int ill, String fname) {
		int nni, nnj, nnn;
		
		try {
		PrintStream out = new PrintStream(new FileOutputStream(fname));
		for( nni = 0; nni < n; nni++ ) {
			nnn = 0; nnj = nns[ill][nni][nnn];
			while( nnj != NO_MORE_NNS ) {
				out.println(nni+" "+nnj);
				nnn++; nnj = nns[ill][nni][nnn];
			}
		}
		out.close();
		} catch( Exception e ) {
			System.err.println("write_out_nn_list failed while writing file '"+fname+"' with exception: "+e);			
		}
	}

	/**
	 * Calculate the current density, in appropriate reduced units.
	 * 
	 * Provides a function to map between the densities and radii and boxes.
	 * 
	 * Calculate the reduced density (wrt close packing) for the hard-sphere system.
	 * <h2>Definitions</h2>
     * <p>
     * For the hard-sphere system, the density is defined for a
     * system of N hard-spheres (in a volume V ) in terms of the reduced
     * number density,
     * </p><p>
     * $\tilde{\rho} = \frac{1}{\tilde{v}} = \frac{\rho}{\rho_{cp}} = \frac{N/V}{\sqrt{2}/\sigma^3}$
     * </p><p>
     * where $\rho_{cp}$ is the number density at close-packing, 
     * $\sigma$ is the diameter of the spheres, and $\tilde{v}$ is the volume
     * per particle (THESIS:ANJ eqn. 4.1). For the NPT ensemble, the pressure is 
     * defined in reduced units such that,
     * </p><p>
     * $\tilde{p} = \frac{p \simga^3}{kT}
     * </p><p>
     * (THESIS:ANJ eqn. 4.2) A factor of $kT$ has been folded into the pressure and so 
     * the temperature only affects the system indirectly, in that the value of the 
     * reduced pressure is dependent upon it.
 	 *</p>
 	 * 
	 * @return The current density, wrt volfrac of close-packed monodisperse hard spheres.
	 */
	protected double calc_density() {
		//System.out.println(" D ALERT calc_density on LS! **************************************");
		//System.out.println(" D diamter= "+diameter+" box[c_lat].z= "+box[c_lat].z+" => "+(n*Math.pow(diameter,3.0)/(box[c_lat].x*box[c_lat].y*box[c_lat].z)));
		return(n*Math.pow(diameter,3.0)/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));
	}
	
	/**
	 * Calculate the current density for the polydisperse system, in appropriate reduced units.
	 * @return The current density, wrt volfrac of close-packed monodisperse hard spheres.
	 */
	protected double calc_polyd_density() {
		double pvol = 0.0;
		for( int i = 0; i < n; i++) {
			pvol += Math.pow(pdiameter[c_lat][i],3.0);
		}
		return(pvol/(Math.sqrt(2.0)*box[c_lat].x*box[c_lat].y*box[c_lat].z));
	}
	
	/**
	 * Calculate the current density.  Switches between mono and polydisperse cases:
	 * @return The current density, in reduced units appropriate to the model.
	 **/
	protected double calc_sys_density() {
		if( polyd ) {
			return calc_polyd_density();
		} else {
			return calc_density();			
		}
	}


	/** ---------------------------------------------------------
	 output_nn_dist:
	 ---------------------------------------------------------*/
	protected void output_nn_dist( boolean showAllAbout0 ) {
		int i,max_steps = 250;
		int[][] t_nn = new int[2][250];
		
		if( showAllAbout0 ) {
			for (i=0; i<max_steps; i++) {
				nn_cut_off = 0.999 + 2.0*(double)i/(double)max_steps;
				make_nn_lists();
				t_nn[0][i] = max_inn[0]; t_nn[1][i] = max_inn[1];
				/*	System.out.printf(" %f %i %i\n",nn_cut_off,max_inn[0],max_inn[1]);*/
				if (i>0) System.out.println(nn_cut_off+" "+(t_nn[0][i]-t_nn[0][i-1])+" "+(t_nn[1][i]-t_nn[1][i-1]));
			}
			
			// This is a test routine, and the code should exit after reporting.
			throw( new LSSimException(LSSimException.ABORTED));
			
		} else {
			// Loop over the nn sets, calculate the seperation as normal, and output that.
			int il,jc,j;
			double dr2;
			String filename;
			for(il = 0; il < 2; il++ ) {
				// Create a sensible filename:
				filename = il+"-"+latt_name[il]+"-NNs.dat";
				try {
					PrintStream out = new PrintStream(new FileOutputStream(filename));
					for(i=0; i<n; i++) {
						jc = 0; j = nns[il][i][jc];
						while ( j!=NO_MORE_NNS ) {
							if ( j < i ) {
								dr2 = ij_sep2(i,j,il,0);
								// Write the seperation out to the file:
								out.println(" "+i+" "+j+" "+Math.sqrt(dr2));
							}
							jc++; j = nns[il][i][jc];
						}
					}
					/* close the file */
					out.close();
				} catch ( Exception e ) {
					System.out.println("output_nn_dist(false) using "+filename+". Failed on exception: "+e);
				}
		    }
			
		}
	}


	/**
	 * Calculates the value of the MCMC order parameter given the lattice energies:
	 * @param energy The energy in the current lattice.
	 * @param Cenergy The energy in the conjugate lattice.
	 * @return The order parameter.
	 */
	protected double calcOrderPar(double energy, double Cenergy ) {
		return (2*c_lat-1)*(Cenergy-energy);
	}
	
	/**
	 * Calculates the current value of the MCMC order parameter.
	 * @return The order parameter.
	 */
	protected double calcOrderPar() {
		return (2*c_lat-1)*(em-e);
	}
	/** -------------------------------------------------------
	 calc_com_vect( VectorD3 comvect ):
	 calculate the com vector and put it in the given array:
	 ------------------------------------------------------- */
	protected VectorD3 calc_com_vect() {
		int i;
		VectorD3 comvect = new VectorD3();
		
		comvect.x = 0.0;
		comvect.y = 0.0;
		comvect.z = 0.0;
		
		for( i=0; i<n; i++ ) {
			/*
			 comvect->x += latt[c_lat][i].x + disp[i].x;
			 comvect->y += latt[c_lat][i].y + disp[i].y;
			 comvect->z += latt[c_lat][i].z + disp[i].z;
			 */
			comvect.x += disp[i].x;
			comvect.y += disp[i].y;
			comvect.z += disp[i].z;
		}
		
		comvect.x /= (double)n;
		comvect.y /= (double)n;
		comvect.z /= (double)n;
	
		return comvect;
	}
	
	/** ---------------------------------------------------------
	 mc_sphere_move:
	 ---------------------------------------------------------*/
	protected void mc_sphere_move(int i) {
		double cond, dr2;
		double dE,dEM,old_m,new_m;
		int im, ss_decision;
		
		/* select sphere at random */
		if( rng instanceof Ran250_521 ) {
			// This range-checking is not required when a well-behaved random number generator, but 250/521 is not well behaved.
			do {
				im = (int)( ((double)(rng.raw())) * ((double)n) );
				if( im < 0 || im >=n )
					System.out.println("ERROR! Chosen particle no. "+im+" out of bounds! Getting another rnd#...");
			} while( im < 0 || im >= n );
		} else {
			im = (int)( ((double)(rng.raw())) * ((double)n) );
		}
		
		/* generate trial move */
		if( move_type == MOVE_RW ) {
			Ddisp.x = (rng.raw()-0.5)*drobox[c_lat][0];
			Ddisp.y = (rng.raw()-0.5)*drobox[c_lat][1];
			Ddisp.z = (rng.raw()-0.5)*drobox[c_lat][2];
			tdisp.x = disp[im].x + Ddisp.x;
			tdisp.y = disp[im].y + Ddisp.y;
			tdisp.z = disp[im].z + Ddisp.z;
		} else if( move_type == MOVE_RWCUT ) {
			do {
				Ddisp.x = (rng.raw()-0.5)*drobox[c_lat][0];
				Ddisp.y = (rng.raw()-0.5)*drobox[c_lat][1];
				Ddisp.z = (rng.raw()-0.5)*drobox[c_lat][2];
				tdisp.x = disp[im].x + Ddisp.x;
				tdisp.y = disp[im].y + Ddisp.y;
				tdisp.z = disp[im].z + Ddisp.z;
				dr2 = tdisp.x*tdisp.x*box[c_lat].x*box[c_lat].x +
				tdisp.y*tdisp.y*box[c_lat].y*box[c_lat].y + 
				tdisp.z*tdisp.z*box[c_lat].z*box[c_lat].z;
			} while( dr2 > dr_limit2 );
		} else if( move_type == MOVE_TOPHAT ) {
			do {
				tdisp.x = (rng.raw()-0.5)*dr_limit/box[c_lat].x;
				tdisp.y = (rng.raw()-0.5)*dr_limit/box[c_lat].y;
				tdisp.z = (rng.raw()-0.5)*dr_limit/box[c_lat].z;
				Ddisp.x = tdisp.x - disp[im].x;
				Ddisp.y = tdisp.y - disp[im].y;
				Ddisp.z = tdisp.z - disp[im].z;
				dr2 = tdisp.x*tdisp.x*box[c_lat].x*box[c_lat].x +
					tdisp.y*tdisp.y*box[c_lat].y*box[c_lat].y + 
					tdisp.z*tdisp.z*box[c_lat].z*box[c_lat].z;
			} while( dr2 > dr_limit2 );
		} else {
			throw( new LSSimException(LSSimException.PARAM,"The requested particle move type is not supported."));
		}
		// Accrue the trail moves counter.
		trial_moves++;
		
		
		dE = calc_dE_sphere_move(im, c_lat);
		cond = -dE*beta;
		
		dEM = calc_dE_sphere_move(im, m_lat);
		old_m = calcOrderPar();
		new_m = calcOrderPar(e+dE, em+dEM);
		
		/* - accept? */
		ss_decision = ss.control( old_m, new_m, cond );
		if( ss_decision == 1 ) {
			e = e + dE;
			em = em + dEM;
			disp[im].x = tdisp.x;
			disp[im].y = tdisp.y;
			disp[im].z = tdisp.z;
			acc_moves++;
			num_De++;
			ave_De += Math.abs(dE);
			// Perform any additional actions required for this type of simulation:
			accepted_sphere_move(im, Ddisp, dE, dEM );
		} else if( ss_decision == -1 ) {
			throw( new LSSimException(LSSimException.PARAM,
					"mc_move_sphere: Move ["+old_m+","+new_m+"] rejected as the move exceeds the diagonal maximum move!"));
		} else if( ss_decision == -2 ) {
			throw( new LSSimException(LSSimException.PARAM,
					"mc_move_sphere: Move ["+old_m+","+new_m+"] rejected as the move exceeds the histogram bounds!"));
		}
		num_e++;
		ave_e += Math.abs(e);
	}
	
	/**
	 * Override this to perform extra computations when a particle move is accepted.
	 * 
	 * @param _im The number of the particle that just got moved·
	 * @param _dR The VectorD3 is was moved by.
	 * @param _dE The energy change (in c_lat) associated with that move.
	 * @param _dEM The conjugate energy change (in m_lat) associated with that move.
	 */
	public void accepted_sphere_move(int _im, VectorD3 _dR, double _dE, double _dEM ) {
	}

	/** ---------------------------------------------------------
	 mc_volume_move():
	 - try to perform a volume move
	 The same box is used for BOTH phases in the NPT case.
	 ---------------------------------------------------------*/
	protected void mc_volume_move() {
		double old_m, new_m, old_e, new_e, new_em;
		double nuVol, nuLinear, oldVol, cond, oldLinear;
		double old_ge, new_ge;
		int ss_decision;
		
		/* Only execute this routine one nth of the time */
		if ( rng.raw() >= one_over_n ) return;
		vol_trial+=1.0;
		
		/* store the old aspect ratio */
		VectorD3 old_box = box[c_lat].copy();
		/* Calculate the old scales */
		oldVol = box[c_lat].x*box[c_lat].y*box[c_lat].z;
		oldLinear = Math.pow(oldVol/(init_boxes[c_lat].x*init_boxes[c_lat].y*init_boxes[c_lat].z),1.0/3.0);
		/* calc/store old orderparameter/energies */
		old_e = e;
		old_m = calcOrderPar();
		old_ge = calc_gs_by_struct(calc_sys_density(), c_lat);

		// Decide whether to use fixed aspect ratio or flexible box moves:
		if( npt_fixed_aspect ) {
			/* create a trial volume move uniform in V */
			nuVol = oldVol + dVol*( (rng.raw()) - 0.5 );
			nuLinear = Math.pow(nuVol/(init_boxes[c_lat].x*init_boxes[c_lat].y*init_boxes[c_lat].z),1.0/3.0);
			setBoxes( init_boxes[c_lat].times(nuLinear) );
		} else {
			// Random-walk each dimension of the box:
			box[c_lat].x += dVol*( (rng.raw()) - 0.5 );
			box[c_lat].y += dVol*( (rng.raw()) - 0.5 );
			box[c_lat].z += dVol*( (rng.raw()) - 0.5 );
			nuVol=box[c_lat].x*box[c_lat].y*box[c_lat].z;
			nuLinear = 0.0; // Not relevant in this case, as the LJ short-cut will not work.
			setBoxes( box[c_lat] );
		}
		/* calculate the order-parameter/energies - capable of using running-average methods: */
		new_e = calc_e_volume_move(c_lat, oldLinear, nuLinear, oldVol, nuVol);
		new_em = calc_e_volume_move(m_lat, oldLinear, nuLinear, oldVol, nuVol);
		
		new_ge = calc_gs_by_struct(calc_sys_density(), c_lat);
		new_m = calcOrderPar(new_e,new_em);
		
		/* calc the acceptance probability */
		cond = -(beta*((new_ge-old_ge)+(new_e-old_e))+BetaPres*(nuVol-oldVol)-n*Math.log(nuVol/oldVol));
		
		/* TODO Remove this debugging stuff 
		System.out.println("mc_vol_move beta="+beta+"+ BetaPres="+BetaPres);
		System.out.println("mc_vol_move "+oldVol+"->"+nuVol);
		System.out.println("mc_vol_move "+(new_ge-old_ge)+"->"+(new_e-old_e));
		*/
		/* accept? - update energies &c */
		ss_decision = ss.control( old_m, new_m, cond );
		if( ss_decision == 1 ) {
			/* accept */
			e = new_e;
			em = new_em;
			/* update acc counter */
			vol_acc+=1.0;
			// Do model-specific accept step:
			accepted_volume_move(oldLinear,nuLinear,oldVol,nuVol);
		} else {
			// Set both boxes to be as the old box:
			setBoxes( old_box );
		}
	}
	protected void accepted_volume_move( double oldLinear, double nuLinear, double oldVol, double nuVol ) {
	}
	

	/** ---------------------------------------------------------
	 mc_basis_flip:
	 - flip to other basis, if we can:
	 ---------------------------------------------------------*/
	protected void mc_basis_flip() {
		double cond, old_m;
		int ss_decision;
		
		if( FLIP_FREQ == 1 ) {
			/* Only execute this routine one nth of the time */
			if ( rng.raw() >= one_over_n ) return;
		} else if( FLIP_FREQ == 2 ) {
			/* Only execute this routine bflipfreq of the time */
			if ( rng.raw() >= bflipfreq ) return;
		}
		
		cond = -(em-e)*beta;
		// OP does not change upon lattice switch.
		old_m = calcOrderPar();
		
		if( NPT_SIM ) {
			/* include the gs energy diff if in the NPT ensemble */
			cond += -(calc_gs_by_struct(calc_sys_density(), m_lat)
					- calc_gs_by_struct(calc_sys_density(), c_lat));
		}
		
		/* pass this to the SS routines */
		ss_decision = ss.control( old_m, old_m, cond );
		if( ss_decision == 1 ) {
			/* accept */
			c_lat = m_lat;
			m_lat = 1 - c_lat;
			//System.out.println("Swapping "+e+" "+em);
			// This does not work - presumably primitives are passed by value:
			//swap_doubles(e, em);
			//System.out.println("Swapped "+e+" "+em);
			double temp = e; e = em; em = temp;
			//System.out.println("Really Swapped "+e+" "+em);
			// Do model-specific accept step:
			accepted_basis_flip();
		} else if( ss_decision == -1 ) {
			/* reject badly! */
			System.out.println("mc_basis_flip: Move rejected, not enuff diagonal space!");
			System.out.println(c_lat+"-"+m_lat+": "+em+" "+e+" ("+(em-e)+")");
		}
		
	}	
	protected void accepted_basis_flip() {
	}

	public void run_simulation() {
		int sweep,i;
				
		/* Main time-sweep loop... */
		for ( sweep = 1; sweep <= total_sweeps; sweep++ ) {
			
			/* If the simulation has equilibrated: */
			if( sweep > equib_sweeps ) {
				/* Count the time (sweeps) spent in each phase: */
				if( samp_type == StrongSampling.MCMC_SAMP ) {
					time_in_phase[c_lat] += ss.getWeightAt(calcOrderPar());
					btime_in_phase[c_lat] += ss.getWeightAt(calcOrderPar());
				} else {
					time_in_phase[c_lat] += 1.0;
					btime_in_phase[c_lat] += 1.0;
				}
			
				/* Periodically output the system status */
				if ( sweep%output_period == 0 ) report_simulation(sweep);
			}
			
		    /* Periodically check the system */
			if ( sweep%check_period == 0) check_simulation(sweep);
			
			/* Sphere loop, where the Monte Carlo actually happens... */
			for ( i = 0; i < n; i++ ) {
				
				/* Volume move */
				if( NPT_SIM ) mc_volume_move();
				
				/* Sphere move */
				mc_sphere_move(i);
				
				/* Basis-flip move */				
				if( LAT_SWITCH ) mc_basis_flip();
				
				/* Perform any additional per-sphere actions for this model */
				simulation_sweep_extra(i);
			}
			
		}
		
	}
	
	
	/**
	 * Extra moves in the sweep:
	 * 
	 * If a particular subclass required more operations at the
	 * per-sphere level, then can go in here:
	 * 
	 * @param i the current particle.
	 */
	protected void simulation_sweep_extra(int i) {
	}
	
	protected void init_simulation() {
		// Initialise the system variables:
		init_variables();
		
		/* output the two structures and then initialise c_lat */
		for( c_lat = 0; c_lat < 2; c_lat++ ) {
			m_lat = 1 - c_lat;
			save_config_xyz(c_lat+"-"+latt_name[c_lat]);	
		}
		c_lat = INIT_PHASE; m_lat = 1 - c_lat;
		
		// Also output the neighbour distributions:
		output_nn_dist(false);
				
		/* calc initial energies... */
		e = calc_e_from_scratch(c_lat);
		em = calc_e_from_scratch(m_lat);
		System.out.println(" Einit "+e+" "+c_lat+" "+em+" "+calcOrderPar());
		
		if( N_SWITCH == 1 ) {
			// Change the hcp-fcc LS into a single-phase NN-switch:
			alter_for_neighbour_switch();
		}
		if( MAKEIT_NSQR == 1 ) {
			make_it_order_n_squared(0);
			make_it_order_n_squared(1);
		}
		write_out_nn_list(0,file_nns+".0.dat");
		write_out_nn_list(1,file_nns+".1.dat");
		/*calc_e_v_dens();*/
		if ( !cold_start ) load_config(file_conf_in);
		write_output_header();
		/* sampling init - including RNG of our choice */
		/* NOTE that the ss_lock is only used when using the Strong Sampling algorithm. */
		ss.initSS(ss_lo,ss_hi,ss_n,ss_nd,ss_lock,samp_type,rng);
		// TODO let the SS set this automatically via the overall max? ss.setExpShift(EXPO_SHIFT);
		ss.loadTPM(file_tpm_in); // load the old TPM in from a standard file: 
		// generate a wf from the tpm:
		if ( samp_type == StrongSampling.MCMC_SAMP ) {
			updateWeightFunction(read_weights);
			ss.outputMCMCWeights(file_mcmc_weights_out);
		}
		
		// reset the TPM to zero, if required:
		if( samp_new_tpm ) ss.setTPMToZero();
		
		ss.printSamplingInfo();
		
	}
	
	protected void write_output_header() {
		System.out.println(" H --Code definitions--");
		System.out.print(" H ensemble = ");
		if ( !NPT_SIM ) System.out.print("NVT\n");
		if ( NPT_SIM ) System.out.print("NPT\n");
		System.out.print(" H sphere move mechanism = ");
		if ( move_type == MOVE_RW ) System.out.print("random walk\n");
		if ( move_type == MOVE_RWCUT ) System.out.print("random walk + cut-off\n");
		if ( move_type == MOVE_TOPHAT ) System.out.print("top-hat\n");
		System.out.print(" H \n");
		System.out.print(" H --User definitions--\n");
		System.out.print(" H Lattice 0 is ("+lsize[0][0]+","+lsize[0][1]+","+lsize[0][2]+") cells = "+n+" spheres\n");
		System.out.print(" H Lattice 1 is ("+lsize[1][0]+","+lsize[1][1]+","+lsize[1][2]+") cells = "+n+" spheres\n");
		System.out.print(" H density = "+density+"\n");
		System.out.print(" H volfrac = "+volfrac+"\n");
		if( polyd ) {
			System.out.print(" H Polydispersity [0]="+polyd_d[0]+" +/- "+polyd_Dd[0]+"%.");
			System.out.print(" H Polydispersity [1]="+polyd_d[1]+" +/- "+polyd_Dd[1]+"%.");
		} else {
			System.out.print(" H  ...diameter = "+diameter+"\n");
		}
		System.out.print(" H beta = "+beta+"\n");
		System.out.print(" H dr = "+dr+"\n");
		if( NPT_SIM ) {
			System.out.print(" H NPT: Pres = "+Pres+", dVol = "+dVol+"\n");
		}
		System.out.print(" H interaction cut off = "+pot_cut_off+"\n");
		System.out.print(" H nn cut off = "+nn_cut_off+"\n");
		if ( !cold_start ) System.out.print(" H loading initial configuration from "+file_conf_in+"\n");
		if ( cold_start ) System.out.print(" H cold start\n");
		System.out.println(" H total_sweeps = "+total_sweeps);
		System.out.println(" H equib_sweeps = "+equib_sweeps);
		System.out.println(" H output_period = "+output_period);
		System.out.println(" H check_period = "+check_period);
		System.out.println(" H iseed = "+iseed);
		System.out.println(" H ");
		System.out.print(" H --Notes--\n");
		System.out.print(" H lat 0: up to "+max_inn[0]+" NNs given a cut-off of "+nn_cut_off+".\n");
		System.out.print(" H lat 1: up to "+max_inn[1]+" NNs given a cut-off of "+nn_cut_off+".\n");
		System.out.println(" H dEgs(fcc-hcp): = "+dEgs+" "+dEgs/((double)n));
		System.out.print(" H There now follows an easily machine-readable list of parameter values\n");
		if ( !NPT_SIM ) System.out.println(" P NVT");
		if ( NPT_SIM ) System.out.println(" P NPT");
		System.out.println(" P 0 "+lsize[0][0]+" "+lsize[0][1]+" "+lsize[0][2]+" "+n);
		System.out.println(" P 1 "+lsize[1][0]+" "+lsize[1][1]+" "+lsize[1][2]+" "+n);
		if ( !NPT_SIM ) System.out.println(" P "+density+" "+beta);
		if ( NPT_SIM ) System.out.println(" P "+density+" "+Pres+" "+beta);
		System.out.println(" P "+0.0 +" "+0.0 +" "+ 0.0/((double)n));
		System.out.println(" P "+dr+" "+pot_cut_off+" "+nn_cut_off);
		System.out.println(" P "+cold_start+" "+total_sweeps+" "+
				equib_sweeps+" "+output_period+" "+check_period);
		System.out.println(" P "+iseed);
		System.out.println(" P "+ss_lo+" "+ss_hi+" "+ss_n+" "+ss_nd+" "+ss_lock+" "+
				samp_type+" "+wf_clip_lo+" "+wf_clip_hi);
	}

	/** ---------------------------------------------------------
	 calc_local_energy:
	 ---------------------------------------------------------*/
	protected double calc_local_energy(int i, int il, int trial) {
		int jc,j;
		double dr2,newe;
		
		newe = 0.0;
		jc = 0; j = nns[il][i][jc];
		while ( j!=NO_MORE_NNS ) {
			dr2 = ij_sep2(i,j,il,trial);
			newe += ij_inter_act(dr2, i, j, il);
			jc++; j = nns[il][i][jc];
		}
		
		return(newe);
	}

	protected void report_simulation(int sweep) {
		// Vector to hold the Centre of Mass:
		VectorD3 comv = new VectorD3();
		
		num_es++; ave_es += em-e;
		double op = calcOrderPar();
		System.out.println(" E "+sweep+" "+e+" "+c_lat+" "+op);// Extra op/n +" "+op/n);
		// Add into block analysers:
		bE[c_lat].add(sweep,e);
		bM[c_lat].add(sweep,op);
		// For NPT simulations:
		if( NPT_SIM ) {
			System.out.println(" V  "+sweep+" "+(box[c_lat].x*box[c_lat].y*box[c_lat].z)+" "+calc_sys_density()+" "+(box[c_lat].z/box[c_lat].y));
			//System.out.println(" Vc "+sweep+" "+(box[m_lat].x*box[m_lat].y*box[m_lat].z));
		}
		/* NOTE: Visualization goes here:
		 if SGL_PLOTS == 1
		 sgl_plot_line((double) sweep, 1.0, (double)(sweep+1), calc_density() );
		 //sgl_plot_fcircle((double)(sweep+1), calc_density(), 0.1);
		  sgl_update();
		  endif
		  */
		if( CALC_COM == 1 ) {
			comv = calc_com_vect();
			System.out.println(" CoM "+comv.x+" "+comv.y+" "+comv.z);
		}
		if( CALC_VIRP == 1 ) {
			System.out.println(" VP "+calc_virial_pressure() );
		}
		
		// In the ideal-gas switch case, update the NN lists:
		if( zero_phase_ig ) {
			make_nn_lists();
		}
	}
	
	protected void check_simulation(int sweep) {
		// Output some status info:
		System.out.println(" C Checking the simulation at sweep "+sweep);
		if( SAVE_CONFIG ) save_config(file_midsim+".conf");
		//save_config_xyz(file_midsim+"-"+sweep);
		//save_config_xyz(m_lat, latt[m_lat], box[m_lat], 0, n, file_midsim+"-conj-"+sweep);
		ss.outputTPM(file_tpm_out);
		ss.outputPDF_TPM(file_pdf_tpm);
		ss.outputPDF_VS(file_pdf_vs);

		// Check the self-consistency of the system:
		check_the_system(LOUDLY);
		
		// Perform density and c/a analysis:
		if( NPT_SIM && sweep > this.equib_sweeps ) {
			// Add to block analysis
			bD[c_lat].add(sweep,calc_sys_density());
			// FIXME This should be a ratio w.r.t. original z/y???
			bCoA[c_lat].add(sweep,box[c_lat].z/box[c_lat].y);
			
			// Output info on densities:
			for( int ila=0; ila < 2; ila++ ) {
				if( bD[ila].hasData() ) {
					System.out.println(" CV lat"+ila+" dens "+bD[ila].getMean()+" +/- "+bD[ila].getStdErr());
				}
				if( bCoA[ila].hasData() ) {
					System.out.println(" CV lat"+ila+" c/a "+bCoA[ila].getMean()+" +/- "+bCoA[ila].getStdErr());
				}
			}
		}
		
		/* Optional: Update using the latest estimate for the w/f */
		if ( evolve_weights && samp_type == StrongSampling.MCMC_SAMP ) {
			updateWeightFunction(false);
			ss.outputMCMCWeights("sweep-"+sweep+"."+file_mcmc_weights_out);
		}

		/* Optional: Zero the displacements, so that the centre of the w/f is sampled rapidly */
		if( reset_displacements_on_check ) {
			for( int i = 0; i < n; i++ ) {
				disp[i].x = 0;
				disp[i].y = 0;
				disp[i].z = 0;
			}
			// Must also zero any incremental counters
			e = calc_e_from_scratch(c_lat);
			em = calc_e_from_scratch(m_lat);
		}
	
		/* Also show the current free-energy difference calc: */
		if( sweep > equib_sweeps && sweep < total_sweeps )
			output_free_energy_difference(sweep);
		
		// Acceptance rate calculations:
		if( AUTO_STEPS == 1 && !zero_phase_ig ) {
			tweak_for_acc_rates();
		}
		/* re-zero the acceptance rate counters */
		acc_moves = 0; trial_moves = 0;
		if( NPT_SIM ) {
			vol_acc = 0; vol_trial = 0;
		}
	}
	
	
	protected void check_the_system(int noise) {
		
		/* Output the acceptance rate */
		if( noise == LOUDLY && trial_moves > 0 )
			System.out.println(" C acc_rate = "+acc_moves/trial_moves);
		if( NPT_SIM && noise == LOUDLY && vol_trial > 0)
			System.out.println(" C dV_acc_rate = "+vol_acc/vol_trial);
		
		
		/* Output the average energy and energy step */
		if (noise == LOUDLY) {
			if( num_e > 0) System.out.println(" C <e> = "+ave_e/num_e);
			if( num_De > 0 ) System.out.println(" C <De> = "+ave_De/num_De);
			if( num_es > 0 ) System.out.println(" C <eS> = "+ave_es/num_es);
		}
		
		
		/* Check the energy running totals*/
		e_check = calc_e_from_scratch(c_lat);
		em_check = calc_e_from_scratch(m_lat);
		if (noise == LOUDLY) {
			System.out.println(" C e = "+e+", e_check = "+e_check);
			System.out.println(" C m = "+em+", m_check = "+em_check);
		}
		
		if ( e_check > 0.0 && Math.abs((e-e_check)/e_check) > E_TOL ) {
			throw( new LSSimException(LSSimException.CHECK,
					"fatal error: (e="+e+") != (e_check="+e_check+") by "+Math.abs(e-e_check)));
		} else {
			e = e_check;
		}
		
		if ( em_check > 0.0 && Math.abs((em-em_check)/em_check) > E_TOL ) {
			throw( new LSSimException(LSSimException.CHECK,
					"fatal error: (em="+em+") != (em_check="+em_check+") by "+Math.abs(em-em_check)));
		} else {
			em = em_check;
		}
		
		// The melt detector - if a particle wanders more that nn_dist away from it's site (inc. CoM motion).
		if( !zero_phase_ig ) {
			VectorD3 comv = calc_com_vect();
			int far_spheres = 0;
			double displacement;
			for( int i = 0; i < n; i++ ) {
				displacement = (disp[i].x-comv.x)*(disp[i].x-comv.x) + (disp[i].y-comv.y)*(disp[i].y-comv.y) + (disp[i].z-comv.z)*(disp[i].z-comv.z);
				if( displacement*displacement > nn_cut_off*nn_cut_off ) far_spheres++;
			}
			// Die with an exception if any spheres are now free of their neighbours.
			if( far_spheres > 0 ) {
				throw( new LSSimException( LSSimException.MELTED, "Crystal has melted! "+far_spheres+" are now > nn_cut_off away from their site(s)!"));
			}
		}		
	}

	protected void finish_simulation() {
		// Check the energy:
		em = calc_e_from_scratch(m_lat);
		check_the_system(LOUDLY);
		// Same the configuration:
		if( SAVE_CONFIG ) save_config(file_final+".conf");
		save_config_xyz(file_final);
		// Finish off the block analysis:
		output_free_energy_difference(total_sweeps);		
		for( int il = 0; il < 2; il++ ) {
			System.out.println(" F Writing out block-analysis results for "+latt_name[il]);
			bE[il].writeBlockAnalysisGraph(latt_name[il]+".E.blk.dat");
			bM[il].writeBlockAnalysisGraph(latt_name[il]+".M.blk.dat");
			if( NPT_SIM ) {
				bD[il].writeBlockAnalysisGraph(latt_name[il]+".dens.blk.dat");
				bCoA[il].writeBlockAnalysisGraph(latt_name[il]+".cOa.blk.dat");
			}
		}
		
		// Output the tail of the analysis:
		output_sim_info_tail();
		
		// Let the SS routines finish up:
		ss.finish();
	}
	
	
	/***********************************************************************
	 * Below are the actual interface methods.  They must be implemented
	 * to create a valid simulation as a sub-class.
	 * 
	 */
	
	protected void output_sim_info_tail() {		
		/* acceptance rate */
		if( trial_moves > 0 )
			System.out.println(" C dr = "+dr+" :: acc_rate = "+acc_moves/trial_moves);
		
		/* average energy and energy step */
		if( num_e != 0.0 )  System.out.println(" C <e> = "+ave_e/num_e);
		if( num_De != 0.0 ) System.out.println(" C <De> = "+ave_De/num_De);
		if( num_es != 0.0 ) System.out.println(" C <eS> = "+ave_es/num_es);
		
		/* free-energy difference calc: */
		output_free_energy_difference(total_sweeps);
	}

	/**
	 * Output the free-energy difference:
	 */
	private void output_free_energy_difference( int sweep ) {
		
		/* Calculate BLOCK Df from Time-In-Phase Data
		 * -------------------------------------------------------- */
		// Block-analyse the time spent in the phases:
		bTIP[0].add(sweep, btime_in_phase[0]);
		bTIP[1].add(sweep, btime_in_phase[1]);
		// Df from block time-spent-in-phase:
		double dfeTIP = Math.log((double)btime_in_phase[1]/(double)btime_in_phase[0])/(double)n;
		// Add this estimate for Df to the block analyser:
		bDfV.add(sweep, dfeTIP);
		bDfV.writeRawDatafile("DfV.obs.dat");
		bDfV.writeBlockAnalysisGraph("DfV.blk.dat");
		// Reset the counters
		btime_in_phase[0] = 0.0; btime_in_phase[1] = 0.0;
		
		/* Calculate BLOCK Df from Peak-Weight Data
		 * -------------------------------------------------------- */
		
		// Also attempt to calculate Df from the StrongSampling block-pdf
		double phw[] = ss.getBlockPDFTotals(ss.getBlockPDF());
		double dfePDF = Math.log(phw[1]/phw[0])/(double)n;
		// Add these estimates for P(0/1) to the block-analysers:
		bPDFw[0].add(sweep, phw[0]);
		bPDFw[1].add(sweep, phw[1]);		
		// Add this estimate for Df to the block-analyser
		bDf.add(sweep, dfePDF);
		// Write out the current block analysis
		bDf.writeRawDatafile("Df.obs.dat");
		bDf.writeBlockAnalysisGraph("Df.blk.dat");
		
		// Reset the PDF so that a new one is collected:
		ss.outputBlockPDF("BlockPDF.dat");
		ss.resetBlockHistogram();
		
		
		/* Output the Dfs via the different calculations, for comparison 
		 * ---------------------------------------------------------------*/
		System.out.println(" F --- Df from whole simulation --- ");
		/* Calculate OVERALL Df from Time-In-Phase Data*/
		System.out.println(" F TIP-VS Df/N= "+Math.log((double)time_in_phase[1]/(double)time_in_phase[0])/(double)n);
		/* Calculate TOTAL Df from Peak-Weight Data */		
		double phwVS[] = ss.getBlockPDFTotals(ss.getPDF_VS());
		System.out.println(" F PDF-VS Df/N= "+(Math.log(phwVS[1]/phwVS[0])/(double)n));

		// Now output data from blocking the simulation:
		System.out.println(" F --- Df from simulation blocks --- ");
		// Output the block-based Df - calculated each way:
		System.out.println(" T BTIP Df/N= "+dfeTIP);
		System.out.println(" F BPDF Df/N= "+dfePDF);
		// Also print out the difference between the TIP and PDF dfe for this block:
		System.out.println(" F B(TIP-PDF) Difference = "+(dfeTIP-dfePDF));
		
		// Estimate the uncertainty in Df so far:
		System.out.println(" F BTIP BA(Df) Df/N= "+bDfV.getMean()+" +/- "+bDfV.getStdErr());		
		System.out.println(" F BPDF BA(Df) Df/N= "+bDf.getMean()+" +/- "+bDf.getStdErr());
		
		// Calculate the values from the averaged weight of PDFs:
		System.out.println(" F BTIP BA(BTIP) Df/N= "+Math.log(bTIP[1].getMean()/bTIP[0].getMean())/(double)n);
		System.out.println(" F BPDF BA(BPDF) Df/N= "+Math.log(bPDFw[1].getMean()/bPDFw[0].getMean())/(double)n);
		System.out.println(" F BPDF BA(BPDF)-A Df/N= "+Math.log(bPDFw[1].getMean()/(1.0-bPDFw[1].getMean()))/(double)n);
		System.out.println(" F BPDF BA(BPDF)-B Df/N= "+Math.log((1.0-bPDFw[0].getMean())/bPDFw[0].getMean())/(double)n);
		
		// Calculate some bounds on the error for BTIP:
		System.out.println(" F BTIP BA(BTIP) Df[lo]/N= "+Math.log((bTIP[1].getMean()-bTIP[1].getStdErr())/(bTIP[0].getMean()+bTIP[0].getStdErr()))/(double)n);
		System.out.println(" F BTIP BA(BTIP) Df[hi]/N= "+Math.log((bTIP[1].getMean()+bTIP[1].getStdErr())/(bTIP[0].getMean()-bTIP[0].getStdErr()))/(double)n);
		
		
	}
	
	/**
	 * The actual pair-potential interaction:
	 * Calculates the pair potential as a function of dr^2.
	 * 
	 * @param dr2 the seperation distance squared.
	 * @param i the identity of one particle
	 * @param j the identity of the other particle (usually i>j)
	 * @param il the lattice index - only required when performing potential-switches etc.
	 * @return The potential for the given distance^2
	 */
	protected abstract double ij_inter(double dr2, int i, int j, int il);
	
	/**
	 * The actual pair-potential interaction in polydisperse systems:
	 * Calculates the pair potential as a function of dr^2 and pdiameter[il][n].
	 * 
	 * @param dr2 the seperation distance squared.
	 * @param i the identity of one particle
	 * @param j the identity of the other particle (usually i>j)
	 * @param il the lattice index - only required when performing potential-switches etc.
	 * @return The potential for the given distance^2
	 */
	protected abstract double ij_inter_poly(double dr2, int i, int j, int il);
		
	/**
	 * The actual pair-potential interaction:
	 * Decides whether this is a polydisperse calculation or not.
	 * 
	 * @param dr2 the seperation distance squared.
	 * @param i the identity of one particle
	 * @param j the identity of the other particle (usually i>j)
	 * @param il the lattice index - only required when performing potential-switches etc.
	 * @return The potential for the given distance^2
	 */
	protected double ij_inter_act(double dr2, int i, int j, int il) {
		// Switch off interactions in the zero phase if desired.
		if( zero_phase_ig && il == 0 ) return 0.0;
		
		// Otherwise, if polydisperse:
		if( polyd ) return ij_inter_poly(dr2,i,j,il);
		
		// Normal mono-disperse case:
		return ij_inter(dr2,i,j,il);
	}	
	
	/** ---------------------------------------------------------
	 calc_e_from_scratch:
	 ---------------------------------------------------------*/
	protected double calc_e_from_scratch(int il) {
	    int i,jc,j,olaps;
	    double dr2;
	
	    olaps = 0;
	    for(i=0; i<n; i++) {
		jc = 0; j = nns[il][i][jc];
		while ( j!=NO_MORE_NNS ) {
		    if ( j < i ) {
			dr2 = ij_sep2(i,j,il,0);
			olaps += (int)ij_inter_act(dr2,i,j,il);
		    }
		    jc++; j = nns[il][i][jc];
		}
	    }
	
	    return olaps;
	}
	
	/** ---------------------------------------------------------
	 calc_gs_by_struct()
	 ---------------------------------------------------------*/
	protected abstract double calc_gs_by_struct(double densi, int latti );

	/**
	 * Calculate the number of overlaps generated by this move:
	 */
	public double calc_dE_sphere_move(int ip, int t_lat) {
		return calc_local_energy(ip,t_lat,TRY_POS) - calc_local_energy(ip,t_lat,CUR_POS);
	}


	/**
	 * Basic routine to calculate E for a given lattice for a new volume.
	 * Can be overridden for systems with short-cuts, e.g. the running-average system for
	 * Lennard Jones.
	 * 
	 * @param t_lat  The lattice id to calculate the energy for.
	 * @param oldLinear The old linear box dimension.
	 * @param nuLinear The new linear box dimension.
	 * @param oldVol The old volume.
	 * @param nuVol The new volume.
	 * @return the total energy at the given volume.
	 */
	protected double calc_e_volume_move( int t_lat, double oldLinear, double nuLinear, double oldVol, double nuVol ) {
		// Default to simple, non running-average calculation:
		return calc_e_from_scratch(t_lat);
	}
	
	/**********************************************************************
	 * Optional extras - if these are not required.
	 */
	
	protected double calc_virial_pressure() {return 0.0;}
	protected void alter_for_neighbour_switch() {}
	protected void make_it_order_n_squared(int ilatt) {}
	protected void virtual_n_squared() {}

	/**
	 * Calculate the volume of a sphere.
	 * @param _diameter The diameter of the sphere.
	 * @return The volume of that sphere.
	 */
	protected double sphere_volume(double _diameter) {
		return 4.0*Math.PI*Math.pow(_diameter/2.0,3.0)/3.0;
	}



}

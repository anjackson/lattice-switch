/**-------------------------------------------------------------
 * jBinLats - LJLS.java
 * net.anjackson.physics.ls.LJLS
 * 
 * Created on 13-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: The GNU General Public License (GPL), as provided in LICENCE.txt.
 * 
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with 
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple 
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.ls;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;

import net.anjackson.maths.*;
import net.anjackson.utils.*;
import net.anjackson.physics.mc.*;
import edu.cornell.lassp.houle.RngPack.*;

/**
 * Lennard-Jones Lattice-Switch Monte Carlo Code
 * 
 * Notes:
 * -FIX ME Split lines parsing correctly? ie index 0/1 depending on leading space.
 * -TA SK: Include the GPL with this software?
 * -TA SK change to sure arrays instead of structs?
 * -TA SK more dynamic allocation of multidimensional arrays?
 * -TA SK change weight file format to XY.
 * -TA SK change to NAG-style Nd arrays?
 * -TA SK NPT to be implemented!
 * -TA SK Load-save to be implemented!
 *  
 * First version: 19th August 1999 : Andrew N Jackson
 * 
 * @author ajackso1
 * @version $Id: LJLSOriginal.java 530 2006-01-31 18:13:19Z anj $
 *
 */
public class LJLSOriginal {
	
	/*--------------------------------
	 Global definitions and variables
	 --------------------------------*/
	
	/* "Hard-coded" system parameters */
	static int NPT_SIM     =1;                              /* 0=NVT ensemble, 1=NPT ensemble */
	static int LAT_SWITCH  =1;                              /* Lattice-Switch off/on 0/1 */
	static int INIT_PHASE  =1;                              /* Initial phase hcp/fcc 0/1 */
		/* Init NNs in NNmode, 0=some, 1=all */
	
	static int VIRT_NSWC   =0;                        /* Use the virtual NN-switch move*/
	static int N_SWITCH    =0;                        /* Turn LS into a neighbour-switch */
	static int INIT_NPHASE =0;                        /* init(+always) lattice-id for NN calc */
	static double EN_SHIFT =0.0;                      /* shifting the excitation energies */
		/* 14.0 for 216, 80 for 1728 */
	
	static int MOVE_TYPE   =1;                        /* 1 = RW, 2 = RW + cut-off */
	static int FLIP_FREQ   =0;                        /* frequency of basis-flip attempts:
		0 = once every move,
		1 = <once every sweep>,
		2 = as defined be bflipfreq */
	static int MAX_NNS       =50;                     /* Max number of nearest neighbours */
	static double E_TOL      =1e-6;                   /* Tolerance of energy checks */
	static double NN_INC_TOL =1e-6;                   /* NN inclusion tolerance */
	static double EXPO_SHIFT =0.0;                    /* exponential shift for no overflows */
	static int SGL_PLOTS     =0;                      /* use SGL graphical output? */
	static int CALC_COM      =0;                      /* Calc+ouput the CoM vector? */
	static int CALC_VIRP     =0;                      /* Calculate the Virial pressure */
	static int AUTO_STEPS    =1;                      /* Automatically tweak the stepsizes */
	static int MAKEIT_NSQR   =0;                      /* Use O(N^2) calculation */
	
	/* useful static flag values */
	static int NO_MORE_NNS =-1;                       /* nnlists: no more nns flag */
	static int LOUDLY      =0;                        /* verbosity: how much output: */
	static int QUIETLY     =1;
	static int CUR_POS     =0;                        /* ij_sep2: calc which seperation: */
	static int TRY_POS     =1;
	static int ZERO_POS    =2;
	
	/* User defined system parameters */
	int     nx,ny,nz;                                 /* dimensions in nos of spheres */
	double  Pres,BetaPres;                            /* pressure, Beta x pressure */
	double  dVol;
	double  density;                                  /* density */
	double  temp;                                     /* inverse temperature */
	double  dr;                                       /* random-walk length parameter */
	double  pot_cut_off,pot_cut_off2;                 /* interaction cut-off distance */
	double  nn_cut_off, max_dr2;                      /* NN inclusion distance */
	boolean cold_start;                               /* true = cold start, false = load old conf */
	int     total_sweeps;                             /* total MCS to run */
	int     equib_sweeps;                             /* equilibration period in MCS */
	int     output_period;                            /* output period in MCS */
	int     check_period;                             /* system checking period in MCS */
	long    iseed;                                    /* random number generator seed */
	
	/* System storage */
	/* Pointers to be malloc'ed */
	VectorD3[][] latt = new VectorD3[2][];            /* pointer to 2*lattice array */
	VectorD3[] disp;                                  /* pointer to displacements array */
	VectorD3 tdisp = new VectorD3();                  /* trial displacements */
	double[] weights;                                 /* pointer to weight function */
	int[][][] nns;                                    /* nearest neighbour table */
	
	/* General simulation variables */
	int     n;                                        /* total number of spheres */
	double  beta;                                     /* inverse temperature */
	double  lj_eta;                                   /* lj leading (energy) coefficient */
	double  e;                                        /* current energy of system */
	double  dEgs;                                     /* ground state energy diff. */
	double  em,m;                                     /* conjugate energy and order param. */
	int     c_lat, m_lat;                             /* identity of the (c)urrent and (m) order-param lattices [=1/2]*/
	double  one_over_n;                               /* value of 1/n (speed-up calc) */
	int     wlo, whi;                                 /* range of the MCMC weight func in m[wlo_m,whi_m] */
	int       wzero;
	VectorD3  init_box = new VectorD3();              /* inital size of system box a,b,c */
	VectorD3  box; double[] box2 = new double[3];     /* current system box a,b,c & ^2 */
	double[]  obox = new double[3];                   /* old box-size, for trial moves */
	double[]  drobox = new double[3];                 /* RW step in each dirn */
	double  diameter,diameter2;                       /* sphere diameter & diameter^2 */
	int[]   max_inn = new int[2];                     /* the observed max no. of NNs */
	int[]   time_in_phase = new int[2];               /* stores time spent in each phase */
	double  bflipfreq = 0.5;                          /* probability of a basis-flip/move */
	double  e12, e6;                                  /* two components of the current energy */
	double  em12, em6;                                /* two components of the conj. energy */
	
	/* things to measure */
	double trial_moves, acc_moves;                    /* stores # of trial/accepted moves */
	double num_e,ave_e;                               /* to calc average energy */
	double num_De,ave_De;                             /* to calc average energy change/move */
	double num_es,ave_es;                             /* to calc average energy cost of switch */
	double vol_acc, vol_trial;
	
	/* sampling routine parameters */
	int     samp_type;                                /* sampling type: strong, mcmc or imp */
	double  ss_lo;                                    /* histogram: lower bound */
	double  ss_hi;                                    /* histogram: upper bound */
	double  wf_clip_lo;                               /* clip weight function up to here */
	double  wf_clip_hi;                               /* clip weight function from here up */
	int     ss_n;                                     /* histogram: no. of bins */
	int     ss_nd;                                    /* tpm: max. +/- matrix diagonal size */
	int     ss_lock;                                  /* ss: period (in moves) for macrolock */
	int     samp_new_tpm;                             /* 1) empty the tpm, 0) use the old tpm */
	int     read_weights;                             /* 1) mcmc uses "weights.in", 0) doesn't */
	RandomSeedable rng;                               /* The random number generator */
	StrongSampling ss = new StrongSampling();         /* The Strong Sampling MC controller */
	int init_nncalc = 0;                              /* Have the NNs been calculated? 0 = no, 1 = yes */

	
	/* automatic step-size tweaking parameters */
	double DrARlo = 0.3, DrARid = 0.35, DrARhi = 0.4, DrARtweak = 0.05;
	double DvARlo = 0.45, DvARid = 0.5, DvARhi = 0.55, DvARtweak = 0.05;
		
	/*----------------------
	 Simulation subroutines
	 ----------------------*/
	
	/* ---------------------------------------------------------
	 read_params_from_user:
	 - read the simulation defining parameters from the user
	 - NB file format is quite strict
	 ---------------------------------------------------------*/
	void read_params_from_user() {
		
		// The system size:
		String dummy = StdinReader.getString();
		int[] nsize = StdinReader.getIntegerArray();
		nx = nsize[0];
		ny = nsize[1];
		nz = nsize[2];
		n = nx*ny*nz;
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		density = StdinReader.getDouble();
		//if NPT_SIM == 0
		diameter = Math.pow(density,1.0/3.0);
		//else
		//diameter = 1.0;
		//endif
		diameter2 = diameter*diameter;
		//System.out.println("Got density = "+density);
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		temp = StdinReader.getDouble();
		beta = 1.0/temp;
		lj_eta = 4.0*beta;
		//System.out.println("Got temp = "+temp);
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		dr = StdinReader.getDouble();
		//System.out.println("Got dr = "+dr);
		
		if( NPT_SIM == 1 ) {
			dummy = StdinReader.getString();
			dummy = StdinReader.getString();
			String ssp[] = StdinReader.getString().split("\\s+");
			Pres = Double.parseDouble(ssp[1]);
			dVol = Double.parseDouble(ssp[2]);
			BetaPres = beta*Pres/density;
		} else {
			dummy = StdinReader.getString();
			dummy = StdinReader.getString();
			pot_cut_off = StdinReader.getDouble();
			pot_cut_off2 = pot_cut_off*pot_cut_off;
		}
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		nn_cut_off = StdinReader.getDouble();
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		cold_start = StdinReader.getBoolean();
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		total_sweeps = StdinReader.getInt();
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		equib_sweeps = StdinReader.getInt();
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		output_period = StdinReader.getInt();
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		check_period = StdinReader.getInt();
		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		iseed = StdinReader.getLong();
		
		/* scanf("%s %lf %lf %i %i %i %i %i %lf %lf\n\n",*/		
		dummy = StdinReader.getString();
		dummy = StdinReader.getString();
		String ssp[] = StdinReader.getString().split("\\s+");
		//for( int i = 0; i < ssp.length; i++ ) System.out.println(" GOT "+i+" "+ssp[i]);
		dummy = ssp[1];		
		if ( "strong".equals(dummy) ) samp_type = StrongSampling.STRONG_SAMP;
		if ( "mcmc".equals(dummy) ) samp_type = StrongSampling.MCMC_SAMP;
		if ( "imp".equals(dummy) ) samp_type = StrongSampling.IMP_SAMP;
		ss_lo        = Double.parseDouble(ssp[2]);
		ss_hi        = Double.parseDouble(ssp[3]);
		ss_n         = Integer.parseInt(ssp[4]);
		ss_nd        = Integer.parseInt(ssp[5]);
		ss_lock      = Integer.parseInt(ssp[6]);
		samp_new_tpm = Integer.parseInt(ssp[7]);
		read_weights = Integer.parseInt(ssp[8]);
		wf_clip_lo   = Double.parseDouble(ssp[9]);
		wf_clip_hi   = Double.parseDouble(ssp[10]);
		
	}
	
	
	/* ---------------------------------------------------------------------
	 read_weights_file:
	 - read-in or construct the multi-canonical weight function
	 - if 'filename' exists, then read it, else use a flat weight function
	 - also, only build to cope with N=216, N=1728 or N=5832
	 - should change to XY format.
	 ---------------------------------------------------------------------*/
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
		wzero = whi/2;
		
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
	
	
	/* ---------------------------------------------------------
	 write_output_header:
	 - write the definition of the simulation parameters to the
	 standard output stream
	 ---------------------------------------------------------*/
	
	void write_output_header() {
		
		System.out.println(" H Lennard-Jones lattice-switch monte carlo code: hcp v fcc:");
		System.out.println(" H ");
		System.out.println(" H --Code definitions--");
		System.out.print(" H ensemble = ");
		if ( NPT_SIM == 0 ) System.out.print("NVT\n");
		if ( NPT_SIM == 1 ) System.out.print("NPT\n");
		System.out.print(" H sphere move mechanism = ");
		if ( MOVE_TYPE == 1 ) System.out.print("random walk\n");
		if ( MOVE_TYPE == 2 ) System.out.print("random walk + cut-off\n");
		if ( MOVE_TYPE == 3 ) System.out.print("top-hat\n");
		System.out.print(" H \n");
		System.out.print(" H --User definitions--\n");
		System.out.print(" H system of ("+nx+","+ny+","+nz+") = "+n+" spheres\n");
		System.out.print(" H density = "+density+"\n");
		System.out.print(" H  ...diameter = "+diameter+"\n");
		System.out.print(" H beta = "+beta+"\n");
		System.out.print(" H lj_eta = 4*beta = +"+lj_eta+"\n");
		System.out.print(" H dr = "+dr+"\n");
		if( NPT_SIM == 1 ) {
			System.out.print(" H NPT: Pres = "+Pres+", dVol = "+dVol+"\n");
		}
		System.out.print(" H interaction cut off = "+pot_cut_off+"\n");
		System.out.print(" H nn cut off = "+nn_cut_off+"\n");
		if ( !cold_start ) System.out.print(" H loading initial configuration from init.conf\n");
		if ( cold_start ) System.out.print(" H cold start\n");
		System.out.println(" H total_sweeps = "+total_sweeps);
		System.out.println(" H equib_sweeps = "+equib_sweeps);
		System.out.println(" H output_period = "+output_period);
		System.out.println(" H check_period = "+check_period);
		System.out.println(" H iseed = "+iseed);
		System.out.println(" H ");
		if( VIRT_NSWC == 1 ) {
			System.out.println(" H VIRTUAL n-switch, structure "+INIT_PHASE+", en_shift = "+EN_SHIFT);
		}
		System.out.print(" H --Notes--\n");
		System.out.print(" H lat 0: up to "+max_inn[0]+" NNs given a cut-off of "+nn_cut_off+".\n");
		System.out.print(" H lat 1: up to "+max_inn[1]+" NNs given a cut-off of "+nn_cut_off+".\n");
		System.out.println(" H dEgs(fcc-hcp): = "+dEgs+" "+dEgs/((double)n));
		System.out.print(" H There now follows an easily machine-readable list of parameter values\n");
		if ( NPT_SIM == 0 ) System.out.println(" P NVT");
		if ( NPT_SIM == 1 ) System.out.println(" P NPT");
		System.out.println(" P "+nx+" "+ny+" "+nz+" "+n);
		if ( NPT_SIM == 0 ) System.out.println(" P "+density+" "+beta);
		if ( NPT_SIM == 1 ) System.out.println(" P "+density+" "+Pres+" "+beta);
		System.out.println(" P "+0.0 +" "+0.0 +" "+ 0.0/((double)n));
		System.out.println(" P "+dr+" "+pot_cut_off+" "+nn_cut_off);
		System.out.println(" P "+cold_start+" "+total_sweeps+" "+
				equib_sweeps+" "+output_period+" "+check_period);
		System.out.println(" P "+iseed);
		System.out.println(" P "+ss_lo+" "+ss_hi+" "+ss_n+" "+ss_nd+" "+ss_lock+" "+
				samp_type+" "+wf_clip_lo+" "+wf_clip_hi);
	}
	
	/*----------------------------------------------------------------
	 make_lattices:
	 - builds the fcc and hcp lattices, incorporating the LS mapping.
	 - requires 6n planes in the z direction for normal pbcs.
	 -----------------------------------------------------------------*/
	void make_lattices() {
		int ilatt,ix,iy,iz,in;
		double xs,ys,zs,t;
		double xo = 0.0,yo,zo,xoo,yoo;
		VectorD3[] shifts = new VectorD3[2];
		shifts[0] = new VectorD3();
		shifts[1] = new VectorD3();
		
		/* Allocate memory for the lattices and displacements */
		latt[0] = new VectorD3[n];
		latt[1] = new VectorD3[n];
		disp = new VectorD3[n];
		
		/* Define offset in the x-dirn between A and B planes */
		t  = Math.sqrt(3.0)/3.0;
		/* Define distances between neighbouring atoms in the same stacking plane */
		xs = Math.sqrt(3.0)/2.0;
		ys = 1.0;
		zs = Math.sqrt(2.0/3.0);
		
		/* Initialise the shift arrays for both lattices */
		for ( ilatt=0; ilatt<2; ilatt++ ) {
			shifts[ilatt].x = 0.0;
			shifts[ilatt].y = 0.0;
			shifts[ilatt].z = 0.0;
		}
		
		/* Initialise particle counter */
		in = 0;
		/* Loop over the x-y planes */
		for ( iz = 1; iz <= nz; iz++ ) {
			/* Define origin of the x-y planes */
			if ( iz%2 == 1 ) xo = 0.0;
			if ( iz%2 == 0 ) xo = +t;
			if ( iz%6 == 1 || iz%6 == 2 ) shifts[1].x = 0.0;
			if ( iz%6 == 3 || iz%6 == 4 ) shifts[1].x = 2.0*t;
			if ( iz%6 == 5 || iz%6 == 0 ) shifts[1].x = 1.0*t;
			yo = 0.0;
			zo = zs * ( (double) (iz-1) );
			/* Loop over x-dirn */
			for ( ix = 1; ix <= nx; ix++ ) {
				xoo = xo + xs*( (double) (ix-1) );
				yoo = yo + ys*( (double) ((ix-1)%2) )/2.0;
				/* Loop over y-dirn */
				for ( iy = 1; iy <= ny; iy++ ) {
					latt[0][in] = new VectorD3();
					latt[0][in].x=xoo;
					latt[0][in].y=yoo+ys*((double)(iy-1));
					latt[0][in].z=zo;
					latt[1][in] = new VectorD3();
					latt[1][in].x=xoo+shifts[1].x;
					latt[1][in].y=yoo+ys*((double)(iy-1));
					latt[1][in].z=zo;
					in++;
				}
			}
		}
		
		/* Calc the cell size at diameter = 1.0 */
		init_box.x = xs*( (double) nx );
		init_box.y = ys*( (double) ny );
		init_box.z = zs*( (double) nz );
		box = new VectorD3();
		if ( NPT_SIM == 1 ) {
			/* if this is an NPT simulation, resize the box to get the desired density */
			box.x = Math.pow(density,-1.0/3.0)*init_box.x;
			box.y = Math.pow(density,-1.0/3.0)*init_box.y;
			box.z = Math.pow(density,-1.0/3.0)*init_box.z;
		} else {
			/* use standard box-size:*/
			box.x = init_box.x;
			box.y = init_box.y;
			box.z = init_box.z;
		}
		/* store box side length ^ 2 */
		box2 = new double[3];
		box2[0] = box.x*box.x;
		box2[1] = box.y*box.y;
		box2[2] = box.z*box.z;
		/* calc RW step in each dirn */
		drobox = new double[3];
		drobox[0] = dr/box.x;
		drobox[1] = dr/box.y;
		drobox[2] = dr/box.z;
		
		/* check that the NN cut-off is actually inside the box */
		for( ix = 0; ix < 3; ix++ ) {
			if (Math.sqrt(box2[ix]) < nn_cut_off) {
				System.err.print("ERROR: This cut-off ("+nn_cut_off+") does not fit in box ("+ix+","+Math.sqrt(box2[ix])+")\n");
				System.exit(1);
			}
		}
		
		/* Map the system into the simulation cube, by scaling, translation, and
		 by applying the periodic boundary conditions */
		/* Loop over lattices: */
		for ( ilatt = 0; ilatt < 2; ilatt++ ) {
			/* Loop over spheres: */
			for ( in = 0; in < n; in++ ) {
				/* Copy vector from lattice array */
				xo = latt[ilatt][in].x;
				yo = latt[ilatt][in].y;
				zo = latt[ilatt][in].z;
				
				/*Scale and translate: */
				xo = ( xo/box.x ) - 0.5;
				yo = ( yo/box.y ) - 0.5;
				zo = ( zo/box.z ) - 0.5;
				
				/*Apply periodic boundary conditions: */
				xo = xo - ( (int) ( xo+xo ) );
				yo = yo - ( (int) ( yo+yo ) );
				zo = zo - ( (int) ( zo+zo ) );
				
				/* Copy modified vector back into the lattice array */
				latt[ilatt][in].x = xo;
				latt[ilatt][in].y = yo;
				latt[ilatt][in].z = zo;
				
				/* Initialise the displacement array */
				disp[in] = new VectorD3();
				disp[in].x = 0.0;
				disp[in].y = 0.0;
				disp[in].z = 0.0;
			}
		}
		
	}
	
	/* -------------------------------------------------------------
	 ij_sep2:
	 - calculates the square of the seperation of particles i & j.
	 - uses nearest image pbcs.
	 -------------------------------------------------------------*/
	double ij_sep2(int i, int j,int il,int trial) {
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
		return(d[0]*d[0]*box2[0] + d[1]*d[1]*box2[1] + d[2]*d[2]*box2[2]);
		
	}
	
	/* ---------------------------------------------------------
	 make_nn_lists:
	 - builds the static neighbour lists:
	 ---------------------------------------------------------*/
	void make_nn_lists() {
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
			System.out.print("ljls: too many NNs; "+max_lnn+" > MAX_NNS("+MAX_NNS+"); exiting...\n");
			System.exit(1);
		}
		
	}
	
	/* ---------------------------------------------------------
	 ---------------------------------------------------------*/
	void load_config(String filename) {
		System.out.println("l_c:  "+filename);
	}
	
	/* ---------------------------------------------------------
	 - save complete system information using the given filename
	 ---------------------------------------------------------*/
	void save_config(String filename) {
		System.out.println("s_c:  "+filename);
		/*
		 fSystem.out.printf(fp,"C %f %f %f\n",latt[c_lat][i].x,latt[c_lat][i].y,latt[c_lat][i].z);
		 */
	}
	
	/* ---------------------------------------------------------
	 save_config_xyz:
	 - in xyz format for use with xmol
	 ---------------------------------------------------------*/
	void save_config_xyz(String filename) {
		int i,ii,ij,ik;
		
		/* open file and outpout number of atoms to come */
		try {
			PrintStream out = new PrintStream(new FileOutputStream(filename));
			out.print(" "+n+"\n\n");
			
			/* loop over all particles and output as carbons */
			for ( i = 0; i < n; i++ ) {
				out.println("C "+(latt[c_lat][i].x+disp[i].x)*box.x+" "+(latt[c_lat][i].y+disp[i].y)*box.y+" "+(latt[c_lat][i].z+disp[i].z)*box.z);
			}
			
			/* Loop over the corners of the pbc box, outputting as hydrogens */
			for ( ii = -1; ii <= 1; ii = ii + 2 ) {
				for ( ij = -1; ij <= 1; ij = ij + 2 ) {
					for ( ik = -1; ik <= 1; ik = ik + 2 ) {
						out.println("H "+box.x*(ii*0.5)+" "+box.y*(ij*0.5)+" "+box.z*(ik*0.5));
					}
				}
			}
			
			/* close the file */
			out.close();
		} catch ( Exception e ) {
			System.out.println("save_config_xyz("+filename+") Failed on exception: "+e);
		}
		
	}
	
	/* ---------------------------------------------------------
	 ij_inter:
	 - calculates the pair potential as a function of dr^2.
	 ---------------------------------------------------------*/
	double ij_inter(double dr2) {
		double invr;
		
		/*    if ( dr2 < pot_cut_off2 && dr2 > 0.0 ) {*/
		if ( dr2 > 0.0 ) {
			invr = diameter2/dr2;
			return( lj_eta*(Math.pow(invr,6.0) - Math.pow(invr,3.0)) );
		} else {
			return(0.0);
		}
		
	}
	
	/* ---------------------------------------------------------
	 ij_inter_pow:
	 - calculates given power-law as a function of dr^2.
	 ---------------------------------------------------------*/
	double ij_inter_pow(double dr2, double rpow) {
		
		if ( dr2 > 0.0 ) {
			return( Math.pow(diameter2/dr2,rpow) );
		} else {
			return(0.0);
		}
		
	}
	
	/* ---------------------------------------------------------
	 calc_density:
	 ---------------------------------------------------------*/
	double calc_density() {
		return(n*Math.pow(diameter,3.0)/(Math.sqrt(2.0)*box.x*box.y*box.z));
	}
	
	/* ---------------------------------------------------------
	 calc_e_from_scratch:
	 ---------------------------------------------------------*/
	double calc_e_from_scratch(int il) {
		int i,jc,j;
		double newe, newe12, newe6, ijsepsqrd, ije12, ije6;
		
		newe6 = 0.0; newe12 = 0.0;
		for(i=0; i<n; i++) {
			jc = 0; j = nns[il][i][jc];
			while ( j!=NO_MORE_NNS ) {
				if ( i < j ) {
					//newe += ij_inter( ij_sep2(i,j,il,CUR_POS) );
					//newe -= ij_inter( ij_sep2(i,j,il,ZERO_POS) );
					ijsepsqrd = ij_sep2(i,j,il,CUR_POS);
					ije6 = ij_inter_pow(ijsepsqrd,3.0); ije12 = ij_inter_pow(ijsepsqrd,6.0);
					newe6 += ije6; newe12 += ije12;
					ijsepsqrd = ij_sep2(i,j,il,ZERO_POS);
					ije6 = ij_inter_pow(ijsepsqrd,3.0); ije12 = ij_inter_pow(ijsepsqrd,6.0);
					newe6 -= ije6; newe12 -= ije12;
				}
				jc++; j = nns[il][i][jc];
			}
		}
		
		// store for the lattices:
		if( il == c_lat) {
			e6 = newe6; e12 = newe12;
		} else {
			em6 = newe6; em12 = newe12;
		}
		
		// calc. return energy:
		newe = lj_eta*(newe12 - newe6);
		
		
		return(newe);
	}
	
	/* ---------------------------------------------------------
	 calc_e_order_n_squared:
	 ---------------------------------------------------------*/
	double calc_e_order_n_squared(int il) {
		int i,j;
		double dr2,newe,store_cut,pot;
		
		store_cut = pot_cut_off;
		pot_cut_off = 1.0e20;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		newe = 0.0;
		for(i=0; i<n; i++) {
			for(j=0; j<n; j++) {
				if ( i > j ) {
					dr2 = ij_sep2(i,j,il,CUR_POS);
					pot = ij_inter(dr2);
					newe += pot;
				}
			}
		}
		pot_cut_off = store_cut;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		
		return(newe);
	}
	
	/* ---------------------------------------------------------
	 calc_e_order_n:
	 ---------------------------------------------------------*/
	double calc_e_order_n(int il) {
		int j;
		double dr2,newe,store_cut,pot;
		
		store_cut = pot_cut_off;
		pot_cut_off = 1.0e20;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		newe = 0.0;
		for(j=0; j<n; j++) {
			if ( j != 0 ) {
				dr2 = ij_sep2(0,j,il,CUR_POS);
				pot = ij_inter(dr2);
				newe += pot;
			}
		}
		pot_cut_off = store_cut;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		
		return(newe);
	}
	
	/* ---------------------------------------------------------
	 init_variables:
	 - initializes various variables and stuff.
	 ---------------------------------------------------------*/
	void init_variables() {
		
		/* initial value of m for 'cold' start */
		m = 0;
		
		/* output the two structures and then initialise c_lat */
		c_lat = 0; m_lat = 1 - c_lat;
		save_config_xyz("hcp.xyz");
		c_lat = 1; m_lat = 1 - c_lat;
		save_config_xyz("fcc.xyz");
		c_lat = INIT_PHASE; m_lat = 1 - c_lat;
		
		/* calc initial energies... */
		e = calc_e_from_scratch(c_lat);
		em = calc_e_from_scratch(m_lat);
		System.out.println("TEST: "+(em-e)/(n*beta));
		
		/* initialise time spent in each phase */
		time_in_phase[0] = 0; time_in_phase[1] = 0;
		
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
		//rng = new RanMT(iseed);
		//rng = new Ranmar(iseed);
		rng = new Ran250_521((int)iseed);
		
	}
	
	/* ---------------------------------------------------------
	 calc_local_energy:
	 ---------------------------------------------------------*/
	double calc_local_energy(int i, int il,int trial) {
		int jc,j;
		double dr2,newe;
		
		newe = 0.0;
		jc = 0; j = nns[il][i][jc];
		while ( j!=NO_MORE_NNS ) {
			dr2 = ij_sep2(i,j,il,trial);
			newe += ij_inter(dr2);
			jc++; j = nns[il][i][jc];
		}
		
		return(newe);
	}
	
	/* ---------------------------------------------------------
	 calc_local_energy_pow:
	 ---------------------------------------------------------*/
	double calc_local_energy_pow(int i, int il,int trial, double rpow) {
		int jc,j;
		double dr2,newe;
		
		newe = 0.0;
		jc = 0; j = nns[il][i][jc];
		while ( j!=NO_MORE_NNS ) {
			dr2 = ij_sep2(i,j,il,trial);
			newe += ij_inter_pow(dr2, rpow);
			jc++; j = nns[il][i][jc];
		}
		
		return(lj_eta*newe);
	}
	
	/* ---------------------------------------------------------
	 ---------------------------------------------------------*/
	void check_the_system(int noise) {
		double em_check,e_check,e_bits_check;
		
		/* Output the acceptance rate */
		if( noise == LOUDLY && trial_moves > 0 )
			System.out.println(" C acc_rate = "+acc_moves/trial_moves);
		if( NPT_SIM == 1 && noise == LOUDLY && vol_trial > 0)
			System.out.println(" C dV_acc_rate = "+vol_acc/vol_trial);
		
		
		/* Output the average energy and energy step */
		if (noise == LOUDLY) {
			if( num_e > 0) System.out.println(" C <e> = "+ave_e/num_e);
			if( num_De > 0 ) System.out.println(" C <De> = "+ave_De/num_De);
			if( num_es > 0 ) System.out.println(" C <eS> = "+ave_es/num_es);
		}
		
		
		/* Check the energy running totals*/
		e_bits_check = lj_eta*(e12 - e6);
		e_check = calc_e_from_scratch(c_lat);
		em_check = calc_e_from_scratch(m_lat);
		if (noise == LOUDLY) {
			System.out.println(" C e = "+e+", e_check = "+e_check);
			System.out.println(" C m = "+em+", m_check = "+em_check);
		}
		
		if ( e_check > 0.0 && Math.abs((e_bits_check-e_check)/e_check) > E_TOL ) {
			System.out.println("fatal error: (e_bits_check="+e_bits_check+") != (e_check="+e_check+") by "+Math.abs(e_bits_check-e_check));
			System.exit(1);
		}
		
		if ( e_check > 0.0 && Math.abs((e-e_check)/e_check) > E_TOL ) {
			System.out.println("fatal error: (e="+e+") != (e_check="+e_check+") by "+Math.abs(e-e_check));
			System.exit(1);
		} else {
			e = e_check;
		}
		
		if ( em_check > 0.0 && Math.abs((em-em_check)/em_check) > E_TOL ) {
			System.out.println("fatal error: (em="+em+") != (em_check="+em_check+") by "+Math.abs(em-em_check));
			System.exit(1);
		} else {
			em = em_check;
		}
		
	}
	
	
	
	/* ---------------------------------------------------------
	 calc_gs_by_struct()
	 ---------------------------------------------------------*/
	double calc_gs_by_struct(double densi, int latti ) {
		
		if ( N_SWITCH == 0 ) {
			return( n*beta*calc_Egs_of_struct(densi, latti) );
		} else {
			// fudge! to get the e_shift:
			double e_shift = 0.0;
			if( latti == 1 ) e_shift = EN_SHIFT;
			return( n*beta*calc_Egs_of_struct(densi, INIT_NPHASE ) + e_shift );
		}
	}
	
	
	/* ---------------------------------------------------------
	 mc_volume_move():
	 - try to perform a volume move, fixed aspect only
	 ---------------------------------------------------------*/
	void mc_volume_move() {
		double old_m, new_m, old_e, new_e, new_em;
		double nuVol, nuLinear, oldVol, cond, oldLinear;
		double old_ge, new_ge;
		int ss_decision;
		
		/* Only execute this routine one nth of the time */
		if ( rng.raw() >= one_over_n ) return;
		vol_trial+=1.0;
		
		/* store the old aspect ratio */
		obox[0] = box.x;
		obox[1] = box.y;
		obox[2] = box.z;
		oldVol = box.x*box.y*box.z;
		oldLinear = Math.pow(oldVol/(init_box.x*init_box.y*init_box.z),1.0/3.0);
		/* calc/store old orderparameter/energies */
		old_e = e;
		old_m = (2*c_lat-1)*(em-e);
		old_ge = calc_gs_by_struct(calc_density(), c_lat);
		
		/* create a trial volume move uniform in V */
		nuVol = oldVol + dVol*( (rng.raw()) - 0.5 );
		nuLinear = Math.pow(nuVol/(init_box.x*init_box.y*init_box.z),1.0/3.0);
		box.x = nuLinear*init_box.x;
		box.y = nuLinear*init_box.y;
		box.z = nuLinear*init_box.z;
		/* Old system, uniform in box dimension L:
		 dbox = 1.0 + dVol*( (rng->get()) - 0.5 );
		 box.x *= dbox;
		 box.y *= dbox;
		 box.z *= dbox;
		 nuVol = box.x*box.y*box.z;
		 */
		
		/* calc box side length ^ 2 */
		box2[0] = box.x*box.x;
		box2[1] = box.y*box.y;
		box2[2] = box.z*box.z;
		
		/* calculate the order-parameter/energies */
		//System.out.printf("ch: %f %f %f, %f %f %f\n",old_e,e12-e6,lj_eta*(e12-e6),old_m,em12,em6);
		//System.out.printf("v: %f %f (%f)\n",oldVol,nuVol,nuLinear/oldLinear);
		new_e = lj_eta*(e12*Math.pow(nuLinear/oldLinear,-12.0)-e6*Math.pow(nuLinear/oldLinear,-6.0));
		new_em = lj_eta*(em12*Math.pow(nuLinear/oldLinear,-12.0)-em6*Math.pow(nuLinear/oldLinear,-6.0));
		//new_e = calc_e_from_scratch(c_lat);
		//new_em = calc_e_from_scratch(m_lat);
		new_m = (2*c_lat-1)*(new_em-new_e);
		new_ge = calc_gs_by_struct(calc_density(), c_lat);
		
		/* calc the acceptance probability */
		cond = -((new_ge-old_ge)+(new_e-old_e)+BetaPres*(nuVol-oldVol)-n*Math.log(nuVol/oldVol));
		
		/* accept? - update energies &c */
		ss_decision = ss.control( old_m, new_m, cond );
		if( ss_decision == 1 ) {
			/* accept */
			e = new_e;
			e6 = e6*Math.pow(nuLinear/oldLinear,-6.0);
			e12 = e12*Math.pow(nuVol/oldVol,-4.0);
			em = new_em;
			em6 = em6*Math.pow(nuVol/oldVol,-2.0);
			em12 = em12*Math.pow(nuVol/oldVol,-4.0);
			/* calc RW step in each dirn */
			drobox[0] = dr/box.x;
			drobox[1] = dr/box.y;
			drobox[2] = dr/box.z;
			/* update acc counter */
			vol_acc+=1.0;
		} else {
			/* reject */
			box.x = obox[0];
			box.y = obox[1];
			box.z = obox[2];
			/* calc box side length ^ 2 */
			box2[0] = box.x*box.x;
			box2[1] = box.y*box.y;
			box2[2] = box.z*box.z;
		}
		
		/*
		 System.out.printf("%i %g: %g: (%g) %g %g: %g %g: %g %g\n",
		 ss_decision,
		 cond,
		 new_e-old_e,
		 0.0,
		 new_e,
		 old_e,
		 new_m,
		 old_m,
		 nuVol,
		 oldVol);
		 System.out.printf("V(%f): e %f %f (%f) em %f %f (%f)\n",vol_trial,e12,e6,e12-e6,em12,em6,em12-em6);
		 check_the_system(LOUDLY);
		 */
		
	}
	
	/* ---------------------------------------------------------
	 mc_sphere_move:
	 ---------------------------------------------------------*/
	void mc_sphere_move(int i) {
		double cond;
		double lmold,lmnew,lenew,leold,dE,dEM,old_m,new_m;
		int im, ss_decision;
		double lmold6,lmnew6,lenew6,leold6,lmold12,lmnew12,lenew12,leold12;
		
		/* - select sphere at random */
		//do {
		im = (int)( ((double)(rng.raw())) * ((double)n) );
		//} while( im < 0 );
		if( im < 0 ) {
			System.out.println("ERK! Chosen particle no. "+im+" < 0!");
		}
		
		/* - generate trial move - currently RW only */
		tdisp.x = disp[im].x + (rng.raw()-0.5)*drobox[0];
		tdisp.y = disp[im].y + (rng.raw()-0.5)*drobox[1];
		tdisp.z = disp[im].z + (rng.raw()-0.5)*drobox[2];
		trial_moves++;
		
		/* - calc overlaps and dm */
		/*
		 lmold = calc_local_energy(im,m_lat,CUR_POS);
		 lmnew = calc_local_energy(im,m_lat,TRY_POS);
		 leold = calc_local_energy(im,c_lat,CUR_POS);
		 lenew = calc_local_energy(im,c_lat,TRY_POS);
		 */
		lmold12 = calc_local_energy_pow(im,m_lat,CUR_POS,6.0);
		lmnew12 = calc_local_energy_pow(im,m_lat,TRY_POS,6.0);
		leold12 = calc_local_energy_pow(im,c_lat,CUR_POS,6.0);
		lenew12 = calc_local_energy_pow(im,c_lat,TRY_POS,6.0);
		lmold6 = calc_local_energy_pow(im,m_lat,CUR_POS,3.0);
		lmnew6 = calc_local_energy_pow(im,m_lat,TRY_POS,3.0);
		leold6 = calc_local_energy_pow(im,c_lat,CUR_POS,3.0);
		lenew6 = calc_local_energy_pow(im,c_lat,TRY_POS,3.0);
		
		lmold = lmold12 - lmold6;
		lmnew = lmnew12 - lmnew6;
		leold = leold12 - leold6;
		lenew = lenew12 - lenew6;
		
		dE = (lenew-leold);
		cond = -dE;
		
		dEM = (lmnew-lmold);
		old_m = (2*c_lat-1)*(em-e);
		new_m = (2*c_lat-1)*(em+dEM-(e+dE));
		
		/* - accept? */
		ss_decision = ss.control( old_m, new_m, cond );
		if( ss_decision == 1 ) {
			e = e + dE;
			em = em + dEM;
			e12 = e12 + (lenew12 - leold12)/lj_eta;
			e6 = e6 + (lenew6 - leold6)/lj_eta;
			em12 = em12 + (lmnew12 - lmold12)/lj_eta;
			em6 = em6 + (lmnew6 - lmold6)/lj_eta;
			disp[im].x = tdisp.x;
			disp[im].y = tdisp.y;
			disp[im].z = tdisp.z;
			acc_moves++;
			num_De++;
			ave_De += Math.abs(dE);
		} else if( ss_decision == -1 ) {
			System.out.println("mc_move_sphere: Move rejected!");
		}
		num_e++;
		ave_e += Math.abs(e);
	}
	
	/* ---------------------------------------------------------
	 mc_basis_flip:
	 - flip to other basis, if we can:
	 ---------------------------------------------------------*/
	void mc_basis_flip() {
		double cond, old_m;
		int ss_decision;
		
		if( FLIP_FREQ == 1 ) {
			/* Only execute this routine one nth of the time */
			if ( rng.raw() >= one_over_n ) return;
		} else if( FLIP_FREQ == 2 ) {
			/* Only execute this routine bflipfreq of the time */
			if ( rng.raw() >= bflipfreq ) return;
		}
		
		cond = -(em-e);
		old_m = (2*c_lat-1)*(em-e);
		
		if( NPT_SIM == 1 ) {
			/* include the gs energy diff if in the NPT ensemble */
			cond += -(calc_gs_by_struct(calc_density(), m_lat)
					- calc_gs_by_struct(calc_density(), c_lat));
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
			//swap_doubles(e6, em6);
			temp = e6; e6 = em6; em6 = temp;
			//swap_doubles(e12, em12);
			temp = e12; e12 = em12; em12 = temp;
		} else if( ss_decision == -1 ) {
			/* reject badly! */
			System.out.println("mc_basis_flip: Move rejected, not enuff diagonal space!");
			System.out.println(c_lat+"-"+m_lat+": "+em+" "+e+" ("+(em-e)+")");
		}
		
	}
	
	/* ---------------------------------------------------------
	 ---------------------------------------------------------*/
	void output_sim_info_tail() {
		double dfe;
		
		/* acceptance rate */
		if( trial_moves > 0 )
			System.out.println(" C dr = "+dr+" :: acc_rate = "+acc_moves/trial_moves);
		
		/* average energy and energy step */
		if( num_e != 0.0 )  System.out.println(" C <e> = "+ave_e/num_e);
		if( num_De != 0.0 ) System.out.println(" C <De> = "+ave_De/num_De);
		if( num_es != 0.0 ) System.out.println(" C <eS> = "+ave_es/num_es);
		
		/* Simply free-energy difference calc: */
		if ( time_in_phase[0] > 0 && time_in_phase[1] > 0 ) {
			System.out.println(" Time spent in hcp = "+time_in_phase[0]);
			System.out.println(" Time spent in fcc = "+time_in_phase[1]);
			dfe = Math.log((double)time_in_phase[1]/(double)time_in_phase[0])/(double)n;
			System.out.println(" T Dfe = "+dfe);
		}
	}
	
	/* ---------------------------------------------------------
	 output_nn_dist:
	 ---------------------------------------------------------*/
	void output_nn_dist() {
		int i,max_steps = 250;
		int[][] t_nn = new int[2][250];
		
		for (i=0; i<max_steps; i++) {
			nn_cut_off = 0.999 + 2.0*(double)i/(double)max_steps;
			make_nn_lists();
			t_nn[0][i] = max_inn[0]; t_nn[1][i] = max_inn[1];
			/*	System.out.printf(" %f %i %i\n",nn_cut_off,max_inn[0],max_inn[1]);*/
			if (i>0) System.out.println(nn_cut_off+" "+(t_nn[0][i]-t_nn[0][i-1])+" "+(t_nn[1][i]-t_nn[1][i-1]));
		}
		
		// This is a test routine, and the code should exit after reporting.
		System.exit(1);
	}
	
	/* ---------------------------------------------------------
	 calc_e_v_dens():
	 - code to calc energy as a function of density
	 ---------------------------------------------------------*/
	void calc_e_v_dens() {
		int i,il,max_steps = 10;
		double[] lel = new double[2], onel = new double[2];
		double store_dens;
		
		store_dens = density;
		
		try {
			PrintStream out = new PrintStream(new FileOutputStream("Ediff.v.dens.out"));
			out.println("# density Elhcp Elfcc Ehcp Efcc dEl dE");
			
			for (i=0; i<max_steps; i++) {
				density= 2.0*(double)(i+1)/(double)max_steps;
				diameter = Math.pow(density,1.0/3.0);
				diameter2 = diameter*diameter;
				out.println(density+" ");
				for (il=0; il<2; il++) {
					lel[il] = calc_local_energy(0,il,CUR_POS);
					onel[il] = calc_e_order_n(il);
				}
				out.println(lel[0]+" "+lel[1]+" "+onel[0]+" "+onel[1]+" "+(lel[1]-lel[0])+" "+(onel[1]-onel[0]));
			}
			out.close();
		} catch( Exception e ) {
			System.err.println("calc_e_v_dens failed while writing file with exception: "+e);
		}
		
		density= store_dens;
		diameter = Math.pow(density,1.0/3.0);
		diameter2 = diameter*diameter;
	}
	
	
	/* -------------------------------------------------------
	 calc_com_vect( VectorD3 comvect ):
	 calculate the com vector and put it in the given array:
	 ------------------------------------------------------- */
	void calc_com_vect( VectorD3 comvect) {
		int i;
		
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
		
	}
	
	/* -------------------------------------------------------
	 calc_virial_pressure()
	 calculate and return the virial pressure:
	 ------------------------------------------------------- */
	double calc_virial_pressure() {
		int i,j,iTot;
		double dr2,fij,rij,TotFdotR;
		double vp;
		
		/* Calculate sum(f.r): */
		TotFdotR = 0.0; iTot = 0;
		for(i=0; i<n; i++) {
			for(j=0; j<n; j++) {
				if ( i != j ) {
					dr2 = ij_sep2(i,j,c_lat,CUR_POS);
					rij = Math.sqrt(dr2);
					fij = -12.0*4.0*Math.pow(diameter,12.0)/Math.pow(rij,13.0)
					+6.0*4.0*Math.pow(diameter,6.0)/Math.pow(rij,7.0);
					TotFdotR += fij*rij;
					iTot++;
				}
			}
		}
		//TotFdotR /= (double)iTot;
		TotFdotR /= -2.0;
		
		/* Calc total virial: */
		vp = Math.sqrt(2.0)*calc_density()/beta +
		Math.pow(diameter,3.0)*TotFdotR/(3.0*box.x*box.y*box.z);
		
		/* check line */
		//System.out.printf(" VPi %f %f %f\n",calc_density()/beta, TotFdotR/(3.0*box.x*box.y*box.z),vp);
		
		return(vp);
	}
	
	/* -------------------------------------------------------
	 tweak_for_acc_rates
	 Changes the step-sizes to make the acc_rates right.
	 ------------------------------------------------------- */
	void tweak_for_acc_rates() {
		double vAR, rAR;
		
		/* If position acc_rate is out of bounds, tweak it */
		rAR = ((double)acc_moves)/((double)trial_moves);
		if( rAR < DrARlo || rAR > DrARhi ) {
			if( rAR > 0.0 && rAR < 1.0 ) {
				dr = dr*Math.log(DrARid)/Math.log(rAR);
			} else {
				if( rAR < DrARlo ) dr = dr*(1.0 + DrARtweak);
				if( rAR > DrARhi ) dr = dr*(1.0 - DrARtweak);
			}
		}
		/* calc RW step in each dirn */
		drobox[0] = dr/box.x;
		drobox[1] = dr/box.y;
		drobox[2] = dr/box.z;
		
		if( NPT_SIM == 1 ) {
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
	
	
	/* -------------------------------------------------------
	 - Write out NN lists
	 ------------------------------------------------------- */
	void write_out_nn_list(int ill, String fname) {
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
	
	/* -------------------------------------------------------
	 Converts to N^2 calc:
	 - Make it O(N^2) - THIS MUST BE CHECKED!!!!!
	 ------------------------------------------------------- */
	void make_it_order_n_squared( int latti ) {
		int ns_i, ns_j, ns_n;
		
		System.out.println(" H Making lattice "+latti+" use a O(N^2) neighbour list...\n");
		
		for( ns_i = 0; ns_i < n; ns_i++ ) {
			ns_n = 0;
			for( ns_j = 0; ns_j < n; ns_j++ ) {
				if( ns_i != ns_j ) {
					nns[latti][ns_i][ns_n] = ns_j;
					ns_n++;
				}
			}
			nns[latti][ns_i][ns_n] = NO_MORE_NNS;
		}
	}
	
	/* -------------------------------------------------------
	 Make the LS into a number of neighbours switch:
	 latt[il] and nns[il][i][inn] to be changed, using
	 NO_MORE_NNS as a end flag and checking that
	 MAX_NNS is large enough:
	 ------------------------------------------------------- */
	void alter_for_neighbour_switch() {
		int this_l, tother_l, ns_i, ns_j, ns_n;
		
		// Let the user know...
		System.out.println(" H THIS IS IN N_SWITCH (nn comparison) mode, for structure "+INIT_NPHASE);
		System.out.println(" H en_shift = "+EN_SHIFT);
		
		// Check that there is room at the inn:
		if( n > MAX_NNS ) {
			System.err.println("ARGH! I need "+n+" neighbours in MAX_NNS!\n");
			System.exit(1);
		}
		
		// Make both structures the same:
		this_l = INIT_NPHASE;
		tother_l = 1 - this_l;
		for( ns_i = 0; ns_i < n; ns_i++ ) {
			latt[tother_l][ns_i].x = latt[this_l][ns_i].x;
			latt[tother_l][ns_i].y = latt[this_l][ns_i].y;
			latt[tother_l][ns_i].z = latt[this_l][ns_i].z;
			// copy nn lists as well!
			ns_n = 0; ns_j = nns[this_l][ns_i][ns_n];
			while ( ns_j != NO_MORE_NNS ) {
				nns[tother_l][ns_i][ns_n] = nns[this_l][ns_i][ns_n];
				ns_n++;  ns_j = nns[this_l][ns_i][ns_n];
			}
			nns[tother_l][ns_i][ns_n] = NO_MORE_NNS;
		}
		
		// Alter the neighbour list for structure 1;
		make_it_order_n_squared(1);
	}
	
	/* -------------------------------------------------------
	 Virtual O(N^2) calc:
	 ------------------------------------------------------- */
	
//	Variable overwrites, memory leaks/overwrites.
	//
//	(int) based pbcs differ from +/- 1.0 (float or double).
	//
//	Are there any floats in the calculation?
//	EN_SHIFT adding in a float 0.0?
	//
//	Note, long-double may not actually be doing anything!
	//
	
	void virtual_n_squared() {
		int ns_i, ns_j, di, edge;
		int[] num_inters = new int[2];
		// NOTE: These were long-double, which is not accessible in java.
		double eonn,ijsepsqrd, ije6, ije12, newe6, newe12;
		double denn, nnprb, ijlatsepsqrd, dl[] = new double[3], dp[] = new double[3], gstot;
		double truncEtot, fullEtot, dennOLD;
		// ANJ End of long-doubles.
		double max_circ_trunc2;
		double ijunscaledsepsqrd, unscaled_boxsqrd = 0.0;
		
		// Using a spherical truncation as big as the cell:
		// It implements a _fixed_length_ cut-off, and so the unscaled cell is used later:
		max_circ_trunc2 = Math.pow((double)nz*0.5*(Math.sqrt(2.0/3.0)),2.0);
		
		// count the number of interactions included:
		num_inters[0] = 0;
		num_inters[1] = 0;
		
		newe6 = 0.0; newe12 = 0.0; gstot = 0.0;
		truncEtot = 0.0; fullEtot = 0.0;
		for( ns_i = 0; ns_i < n; ns_i++ ) {
			for( ns_j = ns_i; ns_j < n; ns_j++ ) {
				if( ns_i != ns_j) {
					
					// Calc relative lattice position:
					dl[0] = latt[c_lat][ns_j].x - latt[c_lat][ns_i].x;
					dl[1] = latt[c_lat][ns_j].y - latt[c_lat][ns_i].y;
					dl[2] = latt[c_lat][ns_j].z - latt[c_lat][ns_i].z;
					
					// Calc relative particle position:
					dp[0]=(latt[c_lat][ns_j].x+disp[ns_j].x)-(latt[c_lat][ns_i].x+disp[ns_i].x);
					dp[1]=(latt[c_lat][ns_j].y+disp[ns_j].y)-(latt[c_lat][ns_i].y+disp[ns_i].y);
					dp[2]=(latt[c_lat][ns_j].z+disp[ns_j].z)-(latt[c_lat][ns_i].z+disp[ns_i].z);
					
					// Shift particle to nearest image of the LATTICE SITE!
					// ...while also building up the square of the latt and particle displmt.
					ijsepsqrd = 0.0; ijlatsepsqrd = 0.0;
					ijunscaledsepsqrd = 0.0;
					edge = 0;
					for( di = 0; di < 3; di++ ) {
						if( dl[di] < -0.5-NN_INC_TOL ) {
							dl[di] += 1.0;
							dp[di] += 1.0;
						} else if( dl[di] > 0.5-NN_INC_TOL ) {
							dl[di] -= 1.0;
							dp[di] -= 1.0;
						}
						ijsepsqrd += dp[di]*dp[di]*box2[di];
						ijlatsepsqrd += dl[di]*dl[di]*box2[di];
						if( di == 0 ) unscaled_boxsqrd = init_box.x*init_box.x;
						if( di == 1 ) unscaled_boxsqrd = init_box.y*init_box.y;
						if( di == 2 ) unscaled_boxsqrd = init_box.z*init_box.z;
						ijunscaledsepsqrd += dl[di]*dl[di]*unscaled_boxsqrd;
						if( Math.abs(Math.abs(dl[di])-0.5) < NN_INC_TOL ) {
							edge = 1;
						}
					}
					if( edge == 0 ) {
						// particle-particle interaction:
						ije6 = ij_inter_pow(ijsepsqrd,3.0); ije12 = ij_inter_pow(ijsepsqrd,6.0);
						newe6 += ije6; newe12 += ije12;
						//System.out.printf("ES %i %i %g  %g %g %g\n",ns_i,ns_j,sqrt(ijsepsqrd),ije6,ije12,lj_eta*(ije12-ije6));
						if( ijunscaledsepsqrd < max_circ_trunc2 ) {
							fullEtot += lj_eta*(ije12 - ije6);
							num_inters[0]++;
						}
						if( ijunscaledsepsqrd < max_dr2 ) {
							truncEtot += lj_eta*(ije12 - ije6);
							num_inters[1]++;
						}
						
						// site-site interaction:
						ijsepsqrd = ijlatsepsqrd;
						ije6 = ij_inter_pow(ijsepsqrd,3.0); ije12 = ij_inter_pow(ijsepsqrd,6.0);
						//System.out.printf("GS %i %i %g  %g %g %g\n",ns_i,ns_j,sqrt(ijsepsqrd),ije6,ije12,lj_eta*(ije12-ije6));
						newe6 -= ije6; newe12 -= ije12;
						if( ijunscaledsepsqrd < max_circ_trunc2 ) fullEtot -= lj_eta*(ije12 - ije6);
						if( ijunscaledsepsqrd < max_dr2 ) truncEtot -= lj_eta*(ije12 - ije6);
						
						gstot += 4.0*(ije12-ije6);
						//System.out.printf("TOT %i,%i %e %e %e\n",ns_i,ns_j,newe6,newe12,newe12-newe6);
					}
				}
			}
		}
		
		eonn = lj_eta*(newe12 - newe6) + EN_SHIFT;
//		dennOLD = (eonn - calc_e_from_scratch(c_lat) );
		dennOLD = (eonn - truncEtot);
		denn = fullEtot - truncEtot;
		nnprb = Math.exp(-denn);
		System.out.println(" NN2 "+dennOLD+" "+denn+" "+nnprb+" "+truncEtot+" "+fullEtot);
		if( init_nncalc == 0 ) {
			System.out.println(" NNTEST "+Math.sqrt(max_circ_trunc2)+" "+ 
					2.0*num_inters[0]/216.0+" "+ 2.0*num_inters[1]/216.0);
			init_nncalc = 1;
		}
		return;
	}
	
	/*-------------------------------------------------------------------------------------
	 - Main program
	 -------------------------------------------------------------------------------------*/	
	
	/**
	 * Main routine - used to invoke this simulation.
	 * @param args The command line arguments
	 */
	public static void main(String[] args) {
		// Instanciate this class:
		LJLSOriginal ljls = new LJLSOriginal();
		
		// Run the simulation:
		ljls.run_simulation();
		
	}	
	public void run_simulation() {
		int sweep,i;
		VectorD3 comv = new VectorD3();
		
		/* Initialising */
		read_params_from_user();
		make_lattices();
		/* Add visualization?
		 #if SGL_PLOTS == 1
		 sgl_init_viewer(800, 300);
		 sgl_init_2d_plot(0.0, 0.5, 10000.0, 2.0);
		 sgl_bg_col(0.0, 0.0, 0.0);
		 sgl_clear();
		 sgl_update();
		 sgl_fg_col(1.0, 1.0, 1.0);
		 #endif
		 */
		comv.x = 0.0; comv.y = 0.0; comv.z = 0.0;
		/*output_nn_dist();*/
		make_nn_lists();
		if( N_SWITCH == 1 ) {
			// Change the hcp-fcc LS into a single-phase NN-switch:
			alter_for_neighbour_switch();
		}
		if( MAKEIT_NSQR == 1 ) {
			make_it_order_n_squared(0);
			make_it_order_n_squared(1);
		}
		write_out_nn_list(0,"nns.0.dat");
		write_out_nn_list(1,"nns.1.dat");
		/*calc_e_v_dens();*/
		init_variables();
		/*read_weights_file("weights");*/
		if ( !cold_start ) load_config("init.conf");
		write_output_header();
		/* sampling init */
		ss.initSS(ss_lo,ss_hi,ss_n,ss_nd,ss_lock,samp_type,rng);
		ss.setExpShift(EXPO_SHIFT);
		ss.loadTPM("TPM.in"); // load the old TPM in from a standard file:
		// generate a wf from the tpm:
		if ( samp_type == StrongSampling.MCMC_SAMP ) {
			ss.getMCMCWeightsFromTPM();
			if ( read_weights == 1 ) ss.loadMCMCWeights("weights.in");
			ss.clipMCMCWeights(wf_clip_lo,wf_clip_hi);
			ss.outputMCMCWeights("mcmc.weights.dat");
		}
		
		System.out.println("Cheese");
		
		// reset the TPM to zero, if required:
		if( samp_new_tpm == 1 ) {
			ss.setTPMToZero();
		}
		ss.printSamplingInfo();
		
		System.out.println(" Einit "+e+" "+c_lat+" "+em+" "+((double)(2*c_lat-1)*(em-e)));
		
		for ( sweep = 1; sweep <= total_sweeps; sweep++ ) {        /* Main loop... */
			/* Count the time spent in each: */
			time_in_phase[c_lat]++;
			
			/* Periodically output the system status */
			if ( sweep%output_period == 0 && sweep > equib_sweeps ) {
				num_es++; ave_es += em-e;
				System.out.println(" E "+sweep+" "+e+" "+c_lat+" "+((double)(2*c_lat-1)*(em-e)));
				if( NPT_SIM == 1 ) {
					/*
					 System.out.printf(" V %i %f %f %f %f\n",sweep,box.x*box.y*box.z,box.x,box.y,box.z);
					 */
					System.out.println(" V "+sweep+" "+(box.x*box.y*box.z)+" "+calc_density());
				}
				/* NOTE: Visualization goes here:
				 if SGL_PLOTS == 1
				 sgl_plot_line((double) sweep, 1.0, (double)(sweep+1), calc_density() );
				 //sgl_plot_fcircle((double)(sweep+1), calc_density(), 0.1);
				  sgl_update();
				  endif
				  */
				if( CALC_COM == 1 ) {
					calc_com_vect(comv);
					System.out.println(" CoM "+comv.x+" "+comv.y+" "+comv.z);
				}
				if( CALC_VIRP == 1 ) {
					System.out.println(" VP "+calc_virial_pressure() );
				}
				if( VIRT_NSWC == 1 ) {
					virtual_n_squared();
				}
			}
			/* Periodically check the system */
			if ( sweep%check_period == 0) {
				save_config("midsim.conf");
				save_config_xyz("midsim.xyz");
				ss.outputTPM("TPM.dat");
				ss.outputPDF_TPM("pdf.tpm.dat");
				ss.outputPDF_VS("pdf.vs.dat");
				check_the_system(LOUDLY);
				if( AUTO_STEPS == 1 ) {
					tweak_for_acc_rates();
				}
				/* re-zero the acceptance rate counters */
				acc_moves = 0; trial_moves = 0;
				if( NPT_SIM == 1 ) {
					vol_acc = 0; vol_trial = 0;
				}
				if( N_SWITCH == 1 ) {
					if( time_in_phase[1] > 0 && time_in_phase[0] > 0 )
						System.out.println(" NN "+time_in_phase[0]+" "+time_in_phase[1]+" "+((double)time_in_phase[1]/(double)time_in_phase[0])+" "+(Math.log((double)time_in_phase[1]/(double)time_in_phase[0])/(double)n));
				}
			}
			
			/* Sphere loop... */
			for ( i = 0; i < n; i++ ) {
				//System.out.printf("e %f %f (%f) em %f %f (%f)\n",e12,e6,e12-e6,em12,em6,em12-em6);
				//check_the_system(LOUDLY);
				
				/* Volume move */
				if( NPT_SIM == 1 ) mc_volume_move();
				
				/* Sphere move */
				mc_sphere_move(i);
				
				/* Basis-flip move */				
				if( LAT_SWITCH == 1 ) mc_basis_flip();
				
			}
			
		}
		
		/* Finish off */
		em = calc_e_from_scratch(m_lat);
		check_the_system(LOUDLY);
		save_config("final.conf");
		save_config_xyz("final.xyz");
		output_sim_info_tail();
		
	}
	
	//----
	
	/*--------------------------------------------------------
	 dEgs_fcchcp:
	 - return the static-lattice energy difference 
	 for a given density [fcc-(ideal)hcp].
	 --------------------------------------------------------*/
	double dEgs_fcchcp(double dens) {
		double A12, A6, Uhcp, Ufcc;
		
		/* unit of distance is the 1st NN seperation */
		A12 = 12.132293768711;  A6  = 14.454897093822;
		Uhcp = 2.0*( Math.pow(dens,4.0)*A12 - Math.pow(dens,2.0)*A6 );
		
		A12 = 12.131880196191;  A6 = 14.453920885450;
		Ufcc = 2.0*( Math.pow(dens,4.0)*A12 - Math.pow(dens,2.0)*A6 );
		
		return(Ufcc-Uhcp);
	}
	
	/*--------------------------------------------------------
	 calc_Egs_of_struct:
	 - return the static-lattice energy 
	 for a given density and structure.
	 0 = hcp, 1 = fcc
	 --------------------------------------------------------*/
	double calc_Egs_of_struct(double dens,int i_lat) {
		double A12, A6, U;
		
		/* unit of distance is the 1st NN seperation */
		if( i_lat == 0 ) { 
			A12 = 12.132293768711;  A6  = 14.454897093822; /* hcp */
		} else { 
			A12 = 12.131880196191;  A6 = 14.453920885450; /* fcc */
		}
		U = 2.0*( Math.pow(dens,4.0)*A12 - Math.pow(dens,2.0)*A6 );
		
		return(U);
	}
	
	
}

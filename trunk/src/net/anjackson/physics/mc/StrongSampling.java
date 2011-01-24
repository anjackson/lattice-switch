/**-------------------------------------------------------------
 * jBinLats - StrongSampling.java
 * net.anjackson.physics.mc.StrongSampling
 * 
 * Created on 13-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.physics.mc;

import net.anjackson.maths.*;
import java.io.*;
import edu.cornell.lassp.houle.RngPack.*;

/**
 * General MC sampling routines.
 * 
 * <p>
 * Supporting:
 * <ul>
 * <li>Metropolis importance sampling,</li>
 * <li>Multi-canonical sampling,</li>
 * <li>Strong-sampling,</li>
 * <li>and, hopefully soon, flat-histogram sampling.</li>
 * </ul>
 * </p>
 * 
 * <p>
 * Issues:
 * <ul>
 * <li>TASK Make weight-function left-right clipping automatic (if two clear maxima in WF, clip far left & far right to those maxima).</li>
 * <li>TASK Store only the n-diagonal elements.</li>
 * <li>TASK Make it possible to evolve the weight-function instead of/as well as the TPM?</li>
 * </ul>
 * </p>
 * 
 * @author ajackso1
 * @version $Id: StrongSampling.java 980 2006-08-21 13:22:35Z anj $
 *
 */
public class StrongSampling implements Serializable {
	/* define names/numbers for the sampling types */
	public static int IMP_SAMP     = 0;
	public static int MCMC_SAMP    = 1;
	public static int STRONG_SAMP  = 2;
	static int FLAT_SAMP           = 3;
	/* weight calculations - shooting or SVD matrix approach */
	public static int SHOOT_WCAL   = 0;
	public static int SVD_WCAL     = 1;
	
	/* barrier type */
	static int BARRIER_TYPE = 1; 
	/* 0 = single moving barrier + locking (BAD)
       1 = double barrier system inc locking (OK)
       2 = single barrier moving slowly outwards */
	
	
	/* general parameters and TP matrix */
	int     ss_type, ss_nsamp, ss_init_samp;
	double  ss_tpm_user_lo, ss_tpm_user_hi;
	double  ss_tpm_lo, ss_tpm_hi, ss_tpm_dn;
	int     ss_tpm_n, ss_tpm_doff, ss_tpm_nd;
	double[][]  ss_tpm;
	/** Arbitrary shift factor applied to weights, alters normalization to keep numbers low and manageable */
	double  ss_expo_shift = 100.0;
	
	/* visited histogram stuff */
	double[] ss_h;
	
	/* strong sampling stuff */
	int    ss_cur_b, ss_bar, ss_bar_time;
	
	/* mcmc sampling stuff */
	double[] ss_mcmc_w;
	
	/* Allow block-analysis of pdf histogram */
	double[] blk_h;
	
	/* single-histogram extrapolation stuff  */
	double ss_she_pdf[], ss_she_b0, ss_she_b1;
	
	/* the random number generator */
	RandomSeedable rngn;
	
	/* analysis stuff - accumulators */
	/** Count the number of moves rejected due to the TPM being too narrow on the diagonal */
	int rejected_dtpm = 0;
	/** Count the number of moves rejected due to the overall histogram (M-range) being too narrow */
	int rejected_hsize = 0;
	
	/*-------------------*/
	/* SAMPLING ROUTINES */
	/*-------------------*/
	
	/*----------------------------------------------------
	 ss_init_ss( lo, hi, n, nd, nsamp, stype);
	 ----------------------------------------------------*/
	public void initSSrng(double lo, double hi, int n, int nd, int nsamp, int stype ) {
		RandomSeedable rng = new RanMT(RandomSeedable.ClockSeed());
		initSS(lo, hi, n, nd, nsamp, stype, rng);
		return;
	}
	/*----------------------------------------------------
	 ss_init_ss( lo, hi, n, nd, nsamp, stype, RanGen);
	 ----------------------------------------------------*/
	public void initSS(double lo, double hi, int n, int nd, int nsamp, int stype, RandomSeedable userrngn) {
		int i;
		
		/* remember the system definitions */
		ss_tpm_user_lo = lo;
		ss_tpm_user_hi = hi;
		ss_tpm_lo      = lo - (hi-lo)/(2.0*(double)n);
		ss_tpm_hi      = hi + (hi-lo)/(2.0*(double)n);
		ss_tpm_n       = n;
		ss_tpm_nd      = 2*nd+1;
		ss_tpm_doff    = nd;
		ss_tpm_dn      = (ss_tpm_hi-ss_tpm_lo)/((double)n);
		ss_nsamp       = nsamp;
		ss_type        = stype;
		
		/* remember the rngn */
		rngn = userrngn;
		
		/* allocate nD arrays */
		ss_tpm    = new double[ss_tpm_n][ss_tpm_n];
		ss_mcmc_w = new double[ss_tpm_n];
		ss_h      = new double[ss_tpm_n];
		blk_h     = new double[ss_tpm_n];
		
		/* zero the tpm */
		setTPMToZero();
		
		/* level the mcmc weights and zero the visited-histogram*/
		for( i=0; i<n; i++ ) {
			ss_mcmc_w[i] = 1.0;
			ss_h[i] = 0.0;
		}
		/* zero the block-analysis helper */
		resetBlockHistogram();
		
		/* tell control routines to initialize */
		ss_init_samp = 1;
		
	    /* Initialise accumulators: */
		rejected_dtpm = 0;
		rejected_hsize = 0;
		
	}
	
	/*----------------------------------------------------
	 void ss_set_exp_shift(double):
	 ----------------------------------------------------*/
	public void setExpShift(double user_shift) {
		ss_expo_shift = user_shift;
	}
	
	/*----------------------------------------------------
	 void ss_finish(void):
	 ----------------------------------------------------*/
	public void finish() {
		 /* Display the number of rejections for both reasons */
		if( rejected_dtpm > 0 )
			System.out.println(" SS:WARNING!!! Rejected "+rejected_dtpm+" moves due to them stepping outside the maximum diagonal ("+ss_tpm_doff+") of the TPM.");
		if( rejected_hsize > 0 )
			System.out.println(" SS:WARNING!!! Rejected "+rejected_hsize+" moves due to them stepping outside the TPM histogram limits ("+ss_tpm_lo+","+ss_tpm_hi+")!");
		// Ensure GC kicks in:
		ss_tpm = null;
		ss_mcmc_w = null;
		ss_h = null;
		blk_h = null;
	}
	
	/*----------------------------------------------------
	 int ss_double_to_discrete(double):
	 double ss_discrete_to_double(in):
	 ----------------------------------------------------*/
	public int doubleToDiscrete(double x) {
		return( (int)((x-ss_tpm_lo)/ss_tpm_dn) );
	}
	public double discreteToDouble(int i) {
		return(ss_tpm_lo + (0.5+i)*ss_tpm_dn);
	}
	
	/*----------------------------------------------------
	 ss_printsampling_info();
	 ----------------------------------------------------*/
	public void printSamplingInfo() {
		
		/* general sampling information */
		if ( ss_type == STRONG_SAMP ) System.out.println(" S sampling type = STRONG_SAMP");
		if ( ss_type == MCMC_SAMP ) System.out.println(" S sampling type = MCMC_SAMP");
		if ( ss_type == IMP_SAMP ) System.out.println(" S sampling type = IMP_SAMP");
		System.out.println(" S "+ss_tpm_user_lo+" "+ss_tpm_user_hi+" "+ss_tpm_n+" "+ss_tpm_doff+" "+ss_nsamp);
		
		/* output weight function */
		/*
		if ( ss_type == MCMC_SAMP ) {
			for ( int i = 0; i < ss_tpm_n; i++ ) {
				System.out.println(" W "+" "+i+" "+ss_mcmc_w[i]);
			}
		}
		*/
		
	}
	
	/*------------------------------------------------------------------
	 int ss_set_tpm_to_zero():
	 
	 ------------------------------------------------------------------*/
	public void setTPMToZero() {
		int i,j;
		for( i=0; i<ss_tpm_n; i++ ) {
			for( j=0; j<ss_tpm_nd; j++ ) {
				ss_tpm[i][j] = 0.0;
			}
		}
	}
	
	
	/*------------------------------------------------------------------
	 int ss_load_tpm(fname):
	 - now copes if the def of the TPM has changed slightly (in range).
	 ------------------------------------------------------------------*/
	public void loadTPM(String fname) {
		String id;
		double lo,hi,tp;
		int n,nd,i,j;
		
		try {
		/* Open the input file */
		BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(fname)));
		
		/* read in definition line */
		String line = in.readLine();
		String[] parts = line.split("\\s+");
		id = parts[0];
		lo = Double.parseDouble(parts[1]);
		hi = Double.parseDouble(parts[2]);
		n  = Integer.parseInt(parts[3]);
		nd = Integer.parseInt(parts[4]);
		// These seem to be ignored, and are assumed to be the same as the
		// simulation settings.  Printing out for completeness:
		System.out.println(" H tpm:"+fname+" - "+id+" M=["+lo+","+hi+"] n="+n+" nd="+nd);
		//fscanf(fp,"%s %lf %lf %i %i",id,&lo,&hi,&n,&nd);
		
		/* zero the tpm */
		setTPMToZero();
		
		/* read in the tpm elements */
		while ( in.ready() ) {
			// Read a line from the file:
			line = in.readLine().trim();
			// Parse it into i, j and tp:
			parts = line.split("\\s+");
			i  = Integer.parseInt(parts[0]);
			j  = Integer.parseInt(parts[1]);
			tp = Double.parseDouble(parts[2]);
			// Load it in:
			i = i - n/2 + ss_tpm_n/2;  // cope with range-change.
			j = j - n/2 + ss_tpm_n/2;
			j = j - i + ss_tpm_doff;   // identify correct diagonal element.
			if( i >= 0 && i < ss_tpm_n && j >= 0 && j < ss_tpm_nd ) {
				//System.out.println("D Adding "+tp+" to TPM at "+i+","+j);
				ss_tpm[i][j] = tp;
			}
		}
		
		/* close the file */
		in.close();
		} catch( Exception e ) {
			System.err.println("ss_load_tpm Failed to read '"+fname+"' - caught exception "+e);
		}
		
	}
	
	
	
	/*------------------------------------------------------------------
	 int ss_output_mcmc_weights(fname):
	 
	 ------------------------------------------------------------------*/
	public void outputMCMCWeights(String fname) {
		int i;
		
		if ( ss_type != MCMC_SAMP ) return;
		
		/* Print weight function (to file?): */
		try {
			PrintStream out = new PrintStream( new FileOutputStream( fname ) );
			for( i=0; i<ss_tpm_n; i++) out.println(" "+i+" "+ss_mcmc_w[i]);
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_output_mcmc_weights Failed due to exception "+e);
		}
		
	}
	
	/* ------------------------------------------------------------------*/
	/**
	 * This clips the weight-function so that it is flat at the 
	 * outer edges (M < lo and M > hi ).  Used to avoid sampling 
	 * unlikely areas other than those between the peaks.
	 */
	public void clipMCMCWeights(double lo, double hi) {
		int i,ilo,ihi;
		double wlo;
		
		if ( ss_type != MCMC_SAMP ) return;
		
		/* transform real window onto index window */
		ilo = doubleToDiscrete(lo);
		ihi = doubleToDiscrete(hi);
		
		/* clip left and right to within given window */
		for(i=0; i<ilo; i++) ss_mcmc_w[i] = ss_mcmc_w[ilo];
		for(i=ihi+1; i<ss_tpm_n; i++) ss_mcmc_w[i] = ss_mcmc_w[ihi];
		
		/* go through and shift the minimum to zero */
		wlo = ss_mcmc_w[ilo];
		for(i=ilo+1; i<=ihi; i++) if ( ss_mcmc_w[i] < wlo ) wlo = ss_mcmc_w[i];
		for(i=0; i<ss_tpm_n; i++) ss_mcmc_w[i] -= wlo;
	}
	
	/* ------------------------------------------------------------------*/
	/**
	 * Automatically infers the best clipping window, choosing the 
	 * highest points in the high and low regions of M.
	 */
	public void clipMCMCWeights() {
		/* Find the zero index: */
		int iz = doubleToDiscrete(0.0);
		
		/* Find the highest point for M<0: */
		int lo = iz;
		for( int i=0; i < iz; i++ ) {
			if( ss_mcmc_w[i] > ss_mcmc_w[lo] ) lo = i; 
		}
		
		/* Find the highest point for M>0: */
		int hi = iz;
		for( int i=hi; i < ss_tpm_n; i++ ) {
			if( ss_mcmc_w[i] > ss_mcmc_w[hi] ) hi = i; 
		}
		
		/* Clip the weight function: */
		clipMCMCWeights(discreteToDouble(lo), discreteToDouble(hi));
	}
	
	/*------------------------------------------------------------------
	 int ss_load_mcmc_weights(fname):
	 
	 ------------------------------------------------------------------*/
	public void loadMCMCWeights(String fname) {
		int i;
		double w;
		boolean old_format = false;
		
		if ( ss_type != MCMC_SAMP ) return;
		
		/* level the mcmc weights */
		for( i=0; i<ss_tpm_n; i++ ) {
			ss_mcmc_w[i] = 1.0;
		}
		
		/* apply weights from file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(fname)));
			String line, parts[];
		/* read in the weights function */
		i = 0;
		while ( in.ready() ) {
			// Read and parse a line:	
			line = in.readLine().trim();
			parts = line.split("\\s+");
			// The old format was just a list of numbers:
			if( old_format ) {
				w = Double.parseDouble(parts[1]);
				ss_mcmc_w[i] = w;
				i++;
			} else {
				i = Integer.parseInt(parts[0]);
				w = Double.parseDouble(parts[1]);
				ss_mcmc_w[i] = w;
			}
		}
		
		/* close the file */
		in.close();
		} catch( Exception e ) {
			throw( new RuntimeException("ss_load_mcmc_weights Failed.",e));
		}
		
		/* Print weight function (to file?): */
		try {
			PrintStream out = new PrintStream( new FileOutputStream("mcmc.weights.dat"));
			for( i=0; i<ss_tpm_n; i++) out.println(" "+i+" "+ss_mcmc_w[i]);
			out.close();
		} catch( Exception e ) {
			throw( new RuntimeException("ss_load_mcmc_weights Failed.",e));
		}
		
	}
	
	
	/*------------------------------------------------------------------
	 int ss_use_these_mcmc_weights(w[n]):
	 
	 ------------------------------------------------------------------*/
	public void setMCMCWeights(double[] use_w) {
		int i;
		
		if ( ss_type != MCMC_SAMP ) return;
		
		/* copy weight function into the array */
		for( i=0; i<ss_tpm_n; i++) ss_mcmc_w[i] = use_w[i];
		
	}
	
	
	
	/*------------------------------------------------------------------
	 int ss_mcmc_weights_from_tpm():
	 
	 ------------------------------------------------------------------*/
	public void getMCMCWeightsFromTPM() {
		int i;
		
		if ( ss_type != MCMC_SAMP ) return;
		
		/* determine weight function from the current tpm */
		setMCMCWeightsFromTPM(ss_mcmc_w, StrongSampling.SHOOT_WCAL);
		//setMCMCWeightsFromTPM(ss_mcmc_w, StrongSampling.SVD_WCAL);
		
		/* Print weight function (to file?): */
		try {
			PrintStream out = new PrintStream( new FileOutputStream("mcmc.weights.dat"));
			for( i=0; i<ss_tpm_n; i++) out.println(" "+i+" "+ss_mcmc_w[i]);
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_mcmc_weights_from_tpm Failed due to exception "+e);
		}

		// Set the exp-shift so that the maximum of the WF lies at 1.0.
		double expo_shift = ss_mcmc_w[0];
		for( i=0; i< ss_mcmc_w.length; i++ ) {
			if( ss_mcmc_w[i] > expo_shift ) expo_shift = ss_mcmc_w[i];
		}
		setExpShift(expo_shift);
	}
	
	
	/*------------------------------------------------------------------
	 int ss_control_mcmc(oldb, newb, expp);
	 
	 Decides whether or not to accept a trial move via mcmc sampling.
	 
	 Arguments and Returns: much as ss_control().
	 
	 ------------------------------------------------------------------*/
	public int controlMCMC(int oldb,int newb,double expp) {
		double mcmc_prb;
		int accept,nb;
		
		/* modify the acceptance probability using the weight function */
		/*    if ( oldb > 0 && oldb < ss_tpm_n-1 && newb > 0 && newb < ss_tpm_n-1 ) {*/
		nb = newb+oldb;
		if ( nb >= 0 && nb < ss_tpm_n-1 ) {
			mcmc_prb = Math.exp(expp + ss_mcmc_w[oldb] - ss_mcmc_w[nb]);
		} else {
			mcmc_prb = 0.0;
		}
		
		/* decide whether to accept or not */
		accept = 0; if ( rngn.raw() < mcmc_prb ) accept = 1;
		
		return(accept);
	}
	
	
	/*----------------------------------------------------
	 int ss_generate_new_barrier( b):
	 picks a new l-r barrier position (bar = -/+1)
	 ----------------------------------------------------*/
	public int generateNewBarrier(int b) {
		int bar = 0,again;
		
		again = 1;
		while ( again == 1 ) {
			/* choose a side at random */
			bar = 2*((int) (2.0*rngn.raw()) ) - 1;
			again = 0;
			/* avoid overstepping the bounds of the ss range */
			if ( b <= 1 && bar == +1 ) again = 1;
			if ( b >= (int)ss_tpm_n - 2 && bar == -1 ) again = 1;
		}
		return(bar);
	}
	
	
	/*------------------------------------------------------------------
	 int ss_control_ss(oldb, newb, expp);
	 
	 Decides whether or not to accept a trial move via strong-sampling.
	 
	 Arguments and Returns: much as ss_control().
	 
	 ------------------------------------------------------------------*/
	public int controlSS(int oldb,int newb,double prb) {
		int accept;
		
		/* initialise, the barrier &c */
		if( ss_init_samp == 1 ) {
			ss_cur_b = oldb;
			if(BARRIER_TYPE == 0 || BARRIER_TYPE == 1) {
				ss_bar = generateNewBarrier(oldb);
			} else if (BARRIER_TYPE == 2) {
				ss_bar = 0;
			}
			ss_init_samp = 0;
			ss_bar_time = 0;
		}
		
		/* decide if the move would be accepted */
		accept = 0; if ( rngn.raw() < prb ) accept = 1;
		
		/* sort out the barrier and locking */
		ss_bar_time++;
		
		/* single moving (RW) barrier with locking */
		if( BARRIER_TYPE == 0 ) {
			/* if no change in macrostate, accept as normal */
			/* if there is a change, then... */
			if( newb != 0 ) {
				/* if the move would cross the barrier, block it */
				if( newb/Math.abs(newb) == ss_bar ) {
					accept = 0;
				}
				/* if the move is good... */
				else {
					/* if the lock-time has not expired */
					if ( ss_bar_time < ss_nsamp ) {
						accept = 0;
						/* if lock has expired, then change the barrier */
					} else if ( accept == 1 ){
						ss_bar = generateNewBarrier(oldb+newb);
						ss_bar_time = 0;
					}
				}
			}
			
			/* double barrier with locking */
		} else if ( BARRIER_TYPE == 1 ) {
			/* if no change in macrostate, accept as normal */
			/* if there is a change, then... */
			if( newb != 0 ) {
				if ( ss_bar == 0 ) { 	/* if we are locked (ss_bar==0)... */
					/* the keep locked until it expires */
					if ( ss_bar_time < ss_nsamp ) {
						accept = 0;
					} else {
						/* upon expiration, pick new target macrostate +-1 */
						accept = 0;
						ss_bar = generateNewBarrier(oldb);
					}
				} else {	                             /* else if we are not locked... */
					/* if the current move will be accepted and yields the target macrostate */
					if ( newb == -ss_bar && accept == 1 ) {
						/* restart the lock */
						ss_bar = 0;
						ss_bar_time = 0;
					} else {
						/* otherwise, do not accept the move */
						accept = 0;
					}
				}
			}
			
			/* single barrier moving slowly outwards */
		} else if ( BARRIER_TYPE == 2 ) {
			/* if no change in macrostate, accept as normal */
			/* if there is a change, then... */
			if( newb != 0 ) {
				/* if the move would cross the barrier, block it */
				if( Math.abs(oldb+newb - doubleToDiscrete(0.0)) > ss_bar ) {
					accept = 0;
				}
				/* if lock has expired, then change the barrier */
				if ( ss_bar_time >= ss_nsamp ) {
					ss_bar++;
					ss_bar_time = 0;
					System.out.println(" SSMB ss_bar = "+ss_bar+" (old="+oldb+", new="+oldb+newb+", doff="+doubleToDiscrete(0.0)+")");
				}
			}
			
		}
		
		return(accept);
	}
	
	
	/*----------------------------------------------------------------------
	 int ss_control_im(oldb, newb, prb);
	 
	 Decides whether or not to accept a trial move via importance sampling.
	 
	 Arguments and Returns: much as ss_control().
	 
	 ----------------------------------------------------------------------*/
	public int controlIMP(int oldb,int newb,double prb) {
		int accept;
		
		/* decide whether to accept or not */
		accept = 0; if ( rngn.raw() < prb ) accept = 1;
		
		return(accept);
	}
	
	
	/*------------------------------------------------------------------
	 int ss_control(oldS, newS, expp):
	 
	 Decides whether or not to accept a candidate move via strong-
	 sampling, multicanonical extended sampling or importance sampling.
	 
	 Arguments:
	 oldS = current position in macrostate space,
	 newS = position in macrostate space if move were accepted,
	 prb = normal (Boltzmann) factor [x where p=min(1,exp(x))].
	 
	 Returns:
	 int = 0 if move should be rejected, and
	 = 1 if move should be accepted,
	 = -1 if the move exceeds the diagonal maximum move,
	 = -2 if the move exceeds the histogram bounds.
	 
	 ------------------------------------------------------------------*/
	public int control(double oldS,double newS,double expp) {
		int oldb, newb, accept = 0;
		double prb;
		
		/* map macrostates onto arrays */
		oldb = doubleToDiscrete(oldS);
		newb = doubleToDiscrete(newS) - oldb;
		
		/* Reject big moves in MCMC or Strong sampling */
		if ( Math.abs(newb) > ss_tpm_doff && ss_type != IMP_SAMP ) {
			rejected_dtpm++;
			return(-1);
		}
		
		/* reject big moves outside desired limits */
		if ( ( newb+oldb < 0 || newb+oldb >= ss_tpm_n) && ss_type != IMP_SAMP ) {
			rejected_hsize++;
			return(-2);
		}
		
		/* calculate move acceptance probability */
		prb = Math.exp(expp); if ( prb > 1.0 ) prb = 1.0;
		/* record move if it is within range */
		if( oldb >= 0 && oldb < ss_tpm_n ) {
			/* record the current status */
			ss_h[oldb]++;
			blk_h[oldb]++;
			
			/* measure the TP matrix if it is within range */
			if( Math.abs(newb) < ss_tpm_doff ) {
				if( newb == 0 ) {
					ss_tpm[oldb][newb+ss_tpm_doff] += 1.0;
				} else {
					ss_tpm[oldb][newb+ss_tpm_doff] += prb;
					ss_tpm[oldb][ss_tpm_doff]      += 1.0 - prb;
				}
			}
		}
		
		/* if strong-sampling, call it */
		if ( ss_type == STRONG_SAMP ) accept = controlSS(oldb,newb,prb);
		
		/* if multicanonical extended sampling, call it */
		if ( ss_type == MCMC_SAMP ) accept = controlMCMC(oldb,newb,expp);
		
		/* if importance sampling, just accept with given probability */
		if ( ss_type == IMP_SAMP ) accept = controlIMP(oldb,newb,prb);
		
		/* pass the decision back */
		return(accept);
	}

	/**
	 * Look up the weight function exp(W) for a particular value of the
	 * order parameter.
	 * 
	 * @param S  The value of the order parameter (a double).
	 * @return The weighting at that value, i.e. exp(W).
	 */
	public double getWeightAt( double S ) {
		// If we are in the histogram, return the weight of interest:
		int nb = doubleToDiscrete(S);
		if( nb >= 0 && nb < ss_tpm_n ) return Math.exp(ss_mcmc_w[nb]-this.ss_expo_shift);
		// Otherwise, evenly weight the results:
		return 1.0;
	}
	
	
	/*----------------------------------------------------
	 ss_weights_from_tpm_via_matrix(*w):
	 
	 ----------------------------------------------------*/
	public void setMCMCWeightsFromTPMViaMatrix(double[] wf) {
		int mw_neq,ip,jp,i,j,nx,neq;
		double a[][],b[],x[],tpm[][];
		double w[],v[][];
		double wmax,wmin,xmin,norm;
		
		/* allocate */
		nx = ss_tpm_n;
		neq = ss_tpm_n*ss_tpm_nd;
		tpm = new double[nx][nx];
		a = new double[neq][nx];
		b = new double[neq];
		x = new double[nx];
		w = new double[nx];
		v = new double[nx][nx];
		for ( j = 1; j <= neq; j++ ) {
			for ( i = 1; i <= nx; i++ ) {
				a[j][i] = 0.0;
			}
		}
		
		/* create normalised tpm matrix */
		for ( i = 0; i < ss_tpm_n; i++ ) { 
			norm = 0.0; for ( j = 0; j < ss_tpm_nd; j++ ) norm += ss_tpm[i][j];
			for ( j = 0; j < ss_tpm_nd; j++ ) {
				if ( norm > 0.0 ) { tpm[i][j] = ss_tpm[i][j]/norm; } else { tpm[i][j] = 0.0; }
			}
		}
		
		/* solve as a matrix problem tc.mw = tr: init: */
		mw_neq = 1;
		
		/* loop over TPM, creating transition connections tc and transition rates tr */
		for ( i = 0; i < ss_tpm_n; i++ ) { 
			for ( j = ss_tpm_doff-1; j < ss_tpm_doff; j++ ) {
				ip = i+j-ss_tpm_doff;
				jp = -j+2*ss_tpm_doff;
				if ( ip > 0 && ip < ss_tpm_n ) {
					if ( tpm[i][j] > 0.0 && tpm[ip][jp] > 0.0 ) {
						a[mw_neq][i+1] = +1;
						a[mw_neq][ip+1] = -1;
						b[mw_neq] = Math.log(tpm[ip][jp]/tpm[i][j]);
						mw_neq++;
					}
				}
			}
		}
		mw_neq--;
		
		/* solve for mw via standard routines */
		SingularValueDecomposition.svdcmp(a,mw_neq,nx,w,v);
		wmax = w[1];
		for ( i = 2; i <= nx; i++ ) if ( w[i] > wmax ) wmax = w[i];
		wmin = wmax*1e-10;
		for ( i = 1; i <= nx; i++ ) if ( w[i] < wmin ) w[i] = 0.0;;
		SingularValueDecomposition.svbksb(a,w,v,mw_neq,nx,b,x);
		
		/* place result into return array, repositioning... */
		xmin = x[1]; for ( i = 2; i <= nx; i++ ) if ( x[i] < xmin ) xmin = x[i];
		for ( i = 1; i <= nx; i++ ) wf[i-1] = (x[i] - xmin);
		
	}
	
	/*----------------------------------------------------
	 ss_weights_from_tpm(*w,flag):
	 
	 ----------------------------------------------------*/
	public void setMCMCWeightsFromTPM(double[] w, int flag) {
		int i,j;
		double s1,s2,norm1,norm2;
		double w1[],w2[],wlow1,wlow2;
		
		if ( flag == SVD_WCAL ) {
			/* solve the matrix way */
			System.out.println(" S solving for weights via svd...");
			setMCMCWeightsFromTPMViaMatrix(w);
			/* save to a file */
			try {
				PrintStream out = new PrintStream( new FileOutputStream("w.matrix.dat"));
				for( i=0; i<ss_tpm_n; i++) out.println(" "+i+" "+w[i]);
				out.close();
			} catch( Exception e ) {
				System.err.println("ss_weights_from_tpm Failed due to exception "+e);
			}
			System.out.println(" done.\n");
			
		} else {
			/* initialise */
			w1 = new double[ss_tpm_n];
			w2 = new double[ss_tpm_n];
			for ( i = 0; i<ss_tpm_n; i++ ) { w1[i] = 0.0; w2[i] = 0.0;}
			
			/* shooting method, L to R */
			wlow1 = 0.0;
			for ( i = 0; i <ss_tpm_n-1; i++ ) {
				norm1 = 0.0; for( j=0; j<ss_tpm_nd; j++ ) norm1 += ss_tpm[i][j];
				norm2 = 0.0; for( j=0; j<ss_tpm_nd; j++ ) norm2 += ss_tpm[i+1][j];
				s1 = 0.0; if ( norm1 > 0.0 ) s1 = ss_tpm[i][+1+ss_tpm_doff]/norm1;
				s2 = 0.0; if ( norm2 > 0.0 ) s2 = ss_tpm[i+1][-1+ss_tpm_doff]/norm2;
				if ( s1 != 0.0 && s2 != 0.0 ) {
					w1[i+1] = w1[i] + Math.log(s1)-Math.log(s2);
				} else {
					w1[i+1] = w1[i];
				}
				if (w1[i+1] < wlow1 ) wlow1 = w1[i+1];
			}
			
			/* shooting method, R to L */
			wlow2 = 0.0;
			for ( i = ss_tpm_n-1; i > 0; i-- ) {
				norm1 = 0.0; for( j=0; j<ss_tpm_nd; j++ ) norm1 += ss_tpm[i][j];
				norm2 = 0.0; for( j=0; j<ss_tpm_nd; j++ ) norm2 += ss_tpm[i-1][j];
				s1 = 0.0; if ( norm1 > 0.0 ) s1 = ss_tpm[i][-1+ss_tpm_doff]/norm1;
				s2 = 0.0; if ( norm2 > 0.0 ) s2 = ss_tpm[i-1][+1+ss_tpm_doff]/norm2;
				if ( s1 != 0.0 && s2 != 0.0 ) {
					w2[i-1] = w2[i] + Math.log(s1)-Math.log(s2);
				} else {
					w2[i-1] = w2[i];
				}
				if ( w2[i-1] < wlow2 ) wlow2 = w2[i-1];
			}
			
			/* put the two shoots together */
			try {
				PrintStream out = new PrintStream( new FileOutputStream("w.shoots.dat"));
				for ( i = 0; i < ss_tpm_n; i++ ) { 
					w1[i] = w1[i] - wlow1;	w2[i] = w2[i] - wlow2;
					if ( i == 0 ) {
						w[i] = w2[i];
					} else if ( i == ss_tpm_n-1 ) {
						w[i] = w1[i];
					} else {
						w[i] = (w1[i] + w2[i])/2.0;
					}
					out.println(i+" "+w1[i]+" "+w2[i]+" "+w[i]);
				}
				out.close();
			} catch( Exception e ) {
				System.err.println("ss_weights_from_tpm.shoot Failed due to exception "+e);
			}
			
		}
		
	}
	
	/*----------------------------------------------------
	 ss_pdf_from_tpm(*pdf):
	 
	 ----------------------------------------------------*/
	public double[] getPDF_TPM() {
		int i;
		double norm;
		
		/* allocate */
		double[] w = new double[ss_tpm_n];
		double[] pdf = new double[ss_tpm_n];
		
		/* first estimate weight function from TPM */
		setMCMCWeightsFromTPM( w, SHOOT_WCAL);
		
		/* then turn the weight function into a pdf */
		norm = 0.0;
		for( i=0; i<ss_tpm_n; i++ ) {
			pdf[i] = Math.exp(w[i]-ss_expo_shift);
			norm += pdf[i]*ss_tpm_dn;
		}
		if ( norm > 0.0 ) { for( i=0; i<ss_tpm_n; i++ ) pdf[i] = pdf[i]/norm; }
		
		return pdf;
	}
	
	
	/*----------------------------------------------------
	 ss_pdf_from_vs(pdf[n]):
	 
	 ----------------------------------------------------*/
	public double[] getPDF_VS() {
		int i;
		double norm;
		double pdf[] = new double[ss_tpm_n];
		
		/* copy vs data, add in weights if necessary*/
		norm = 0.0;
		for( i=0; i<ss_tpm_n; i++ ) {
			if ( ss_type == MCMC_SAMP ) {
				pdf[i] = ss_h[i]*Math.exp(ss_mcmc_w[i]-ss_expo_shift);
			} else {
				pdf[i] = ss_h[i];
			}
			norm += pdf[i]*ss_tpm_dn;
		}
		if ( norm > 0.0 ) { for( i=0; i<ss_tpm_n; i++ ) pdf[i] /= norm; }
		
		return pdf;
	}
	
	
	/*----------------------------------------------------
	 ss_output_vs_pdf(fname):
	 
	 ----------------------------------------------------*/
	public void outputPDF_VS(String fname) {
		int i;
		double[] pdf;
		
		/* find vs pdf */
		pdf = getPDF_VS();
		
		/* write out to file */
		try {
			PrintStream out = new PrintStream( new FileOutputStream(fname));
			out.println("# "+ss_tpm_lo+" "+ss_tpm_hi+" "+((int)ss_tpm_n));
			for( i=0; i<ss_tpm_n; i++ )
				out.println(discreteToDouble(i)+" "+pdf[i]+" "+ss_h[i]+" "+ss_mcmc_w[i]);
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_output_vs_pdf Failed due to exception "+e);
		}
		
	}
	
	
	/*----------------------------------------------------
	 ss_output_tpm_pdf(fname):
	 
	 ----------------------------------------------------*/
	public void outputPDF_TPM(String fname) {
		int i;
		double[] pdf;
		
		/* build pdf from tpm data */
		pdf = getPDF_TPM();
		
		/* write out to file */
		try {
			PrintStream out = new PrintStream( new FileOutputStream(fname));
			out.println("# "+ss_tpm_lo+" "+ss_tpm_hi+" "+((int)ss_tpm_n));
			for( i=0; i<ss_tpm_n; i++ ) out.println(discreteToDouble(i)+" "+pdf[i]);
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_output_tpm_pdf Failed due to exception "+e);
		}
		
	}
	
	
	/*----------------------------------------------------
	 ss_output_tpm(fname):
	 
	 ----------------------------------------------------*/
	public void outputTPM(String fname) {
		int i,j;
		
		try {
			PrintStream out = new PrintStream( new FileOutputStream(fname));
			
			out.println("#TPM "+ss_tpm_lo+" "+ss_tpm_hi+" "+((int)ss_tpm_n)+" "+ss_tpm_nd);
			
			for( i=0; i<ss_tpm_n; i++ ) {
				for( j=0; j<ss_tpm_nd; j++ ) {
					if( ss_tpm[i][j] > 0.0 ) 
						out.println(i+" "+(i+j-ss_tpm_doff)+" "+ss_tpm[i][j]);
				}
			}
			
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_output_tpm Failed due to exception "+e);
		}
		
		
	}
	
	/*----------------------------------------------------
	 ss_output_normalized_tpm(fname):
	 
	 ----------------------------------------------------*/
	public void outputNormTPM(String fname) {
		int i,j;
		double norm,output;
		
		try {
			PrintStream out = new PrintStream( new FileOutputStream(fname));
			
			out.println("#TPM "+ss_tpm_lo+" "+ss_tpm_hi+" "+((int)ss_tpm_n)+" "+ss_tpm_nd);
			
			for( i=0; i<ss_tpm_n; i++ ) {
				norm = 0.0; for( j=0; j<ss_tpm_nd; j++ ) norm += ss_tpm[i][j];
				
				for( j=0; j<ss_tpm_nd; j++ ) {
					output = 0.0; if ( norm > 0.0 ) output = ss_tpm[i][j]/norm;
					if ( output > 0.0 ) out.println(i+" "+(i+j-ss_tpm_doff)+" "+output);
				}
			}
			
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_output_normalized_tpm Failed due to exception "+e);
		}
		
	}
	
	
	/*-----------------------------------------*/
	/* SINGLE HISTOGRAM EXTRAPOLATION ROUTINES */
	/*-----------------------------------------*/
	
	public void   initShe(double b0, double b1) {
		int i;
		
		/* claim and zero the extrapolate-pdf */
		ss_she_pdf = new double[ss_tpm_n];
		ss_she_b0 = b0; ss_she_b1 = b1;
		for (i=0; i<ss_tpm_n; i++) ss_she_pdf[i] = 0.0;
		
	}
	
	public void   controlShe(double old, double ener) {
		int mbin;
		
		/* find the right spot in the pdf */
		mbin = doubleToDiscrete(old);
		
		/* increment the extrap-pdf */
		if ( ss_type == MCMC_SAMP ) {
			ss_she_pdf[mbin] += Math.exp((ss_she_b1 - ss_she_b0)*ener+ss_mcmc_w[mbin]-ss_expo_shift);
		}
		
	}
	
	public void   outputShePDF_VS(String fname) {
		int i;
		double norm,output;
		
		try {
			PrintStream out = new PrintStream( new FileOutputStream(fname));
			
			norm = 0.0; for (i=0; i<ss_tpm_n; i++) norm += ss_she_pdf[i]*ss_tpm_dn;
			
			for (i=0; i<ss_tpm_n; i++) {
				output = 0.0;
				if ( norm > 0.0 ) output = ss_she_pdf[i]/norm;
				out.println(discreteToDouble(i)+" "+output);
			}
			
			out.close();
		} catch( Exception e ) {
			System.err.println("ss_output_vs_she_pdf Failed due to exception "+e);
		}
		
	}

	
	/*------------------- Block-analysis helper routines -----------------*/

	/**
	 * This is used to reset the block histogram to zero.
	 */
	public void resetBlockHistogram() {
		for( int i=0; i<ss_tpm_n; i++ )
			blk_h[i] = 0.0;
	}
	
	/**
	 * This can be used to read the current state of the block histogram.
	 * TASK Remove the debugging code from this routine when all is well.
	 * @return An array of doubles containing the block histogram.
	 */
	public double[] getBlockPDF() {
		int i;
		double norm, pdf[] = new double[ss_tpm_n];
		
		// Also output the current un-normalised plot to a file:
		String filename = "BlockPDF-current-raw.dat";
		try {
			PrintStream out = new PrintStream( new FileOutputStream(filename));
			
			/* copy vs data, add in weights if necessary*/
			norm = 0.0;
			for( i=0; i<ss_tpm_n; i++ ) {
				if ( ss_type == MCMC_SAMP ) {
					pdf[i] = blk_h[i]*Math.exp(ss_mcmc_w[i]-ss_expo_shift);// IDEA Remove this expo-shift?
				} else {
					pdf[i] = blk_h[i];
				}
				norm += pdf[i]*ss_tpm_dn;
				// Also output the un-normalised and raw data, for testing.
				out.println(i+" "+pdf[i]+" "+blk_h[i]+" "+Math.exp(ss_mcmc_w[i])+" "+
						Math.exp(ss_mcmc_w[i]-ss_expo_shift));
			}
			if ( norm != 0.0 ) { for( i=0; i<ss_tpm_n; i++ ) pdf[i] /= norm; }
			
			out.close();
		} catch( Exception e ) {
			System.err.println("outputBlockPDF("+filename+") Failed due to exception "+e);
		}
		return pdf;
	}
	
	/**
	 * This measures the normalized weight of the un-weighted distribution
	 * either side of zero.  Used for comparing phases, calculating Df/Dg.
	 * @return double array containing the weights [(-ve peak), (+ve peak)]
	 */
	public double[] getBlockPDFTotals(double[] pdf) {
		// Find the bin corresponding to zero
		int zero_bin = doubleToDiscrete(0.0);
		// Count the weight of each side of zero:
		double zw[] = new double[2];
		zw[0] = 0.0;
		zw[1] = 0.0;
		// Add up:
		for( int i=0; i<zero_bin; i++) zw[0] += pdf[i];
		zw[0] += 0.5*pdf[zero_bin];
		zw[1] += 0.5*pdf[zero_bin];
		for( int i=zero_bin+1; i<ss_tpm_n; i++) zw[1] += pdf[i];
		// Correct the normalization so the totals go to 1.0
		zw[0] *= ss_tpm_dn;
		zw[1] *= ss_tpm_dn;
		
		// Pass the weights back.
		return zw;
	}

	/**
	 * Output the block PDF to the specified file:
	 * @param filename to write the block PDF into.
	 */
	public void outputBlockPDF(String filename) {
		// Get the PDF from the block histogram
		double pdf[] = getBlockPDF();
		try {
			PrintStream out = new PrintStream( new FileOutputStream(filename));						
			for (int i=0; i<pdf.length; i++)
				out.println(discreteToDouble(i)+" "+pdf[i]);
			
			out.close();
		} catch( Exception e ) {
			System.err.println("outputBlockPDF("+filename+") Failed due to exception "+e);
		}
		
	}

}

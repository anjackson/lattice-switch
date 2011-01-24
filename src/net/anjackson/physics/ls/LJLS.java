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

import java.io.FileOutputStream;
import java.io.PrintStream;

import net.anjackson.utils.ParamFileReader;

/**
 * Lennard-Jones Lattice-Switch Monte Carlo Code
 * <p>
 * Notes:
 * <ul>
 * <li>TODO The LJ code really should be tested!</li>
 x* </ul>
 * 
 * <p>
 * First version: 19th August 1999 : Andrew N Jackson
 * </p>
 * 
 * @author ajackso1
 * @version $Id: LJLS.java 963 2006-08-14 21:02:10Z anj $
 *
 */
public class LJLS extends LatticeSwitchSimulation {
	
	/*--------------------------------
	 Global definitions and variables
	 --------------------------------*/
	
	/** init(+always) lattice-id for NN calc */
	static final int INIT_NPHASE =0;                        
	/** shifting the excitation energies */
	static double EN_SHIFT =0.0;                      
		/* 14.0 for 216, 80 for 1728 */
	/** Use the virtual NN-switch move */
	protected static final int VIRT_NSWC = 0;
	/** Option to output the E v. Density data */
	protected static final boolean OUTPUT_E_V_DENSITY = false;

	
	/** lj leading (energy) coefficient */
	protected static final double lj_eta = 4.0;                 
	/** Power-decomposition variables for this potential: two components of the current energy */
	protected double  e12, e6;                                  
	/** Power-decomposition variables for this potential: two components of the conj. energy */
	protected double  em12, em6;                                
	

	/**
	 * read parameters, both general and model-specific:
	 */
	protected void read_params_from_user( ParamFileReader pars ) {
		// Read general parameters:
		super.read_params_from_user(pars);
		
		// Read model parameters (optional):
		if( !NPT_SIM ) {
			try {
				if( pars.isDefined("pot.cut-off") ) pot_cut_off = pars.getDouble("pot.cut-off");
			} catch( ParamFileReader.UndefinedParameterException e ) {
				// If any of them is not defined, cry and grumble.
				throw( new LSSimException(LSSimException.PARAM,
					"Some required parameters were missing: pot.cut-off"));
			}
			pot_cut_off2 = pot_cut_off*pot_cut_off;
		}		
	}
		
	/** ---------------------------------------------------------
	 write_output_header:
	 - write the definition of the simulation parameters to the
	 standard output stream
	 ---------------------------------------------------------*/
	protected void write_output_header() {		
		System.out.println(" H Lennard-Jones lattice-switch monte carlo code: hcp v fcc:");
		System.out.println(" H ");
		// Write out the standard parts of the LSMC header:
		super.write_output_header();
		// Write out LJLS-specific section:
		System.out.println(" H For the Lennard-Jones Lattice-Switch: ");
		System.out.print(" H lj_eta = 4*beta = +"+lj_eta+"\n");
		if( VIRT_NSWC == 1 ) {
			System.out.println(" H VIRTUAL n-switch, structure "+INIT_PHASE+", en_shift = "+EN_SHIFT);
		}
		
		// Optionally output E v. density:
		if( OUTPUT_E_V_DENSITY ) calc_e_v_dens();

	}
		
	/** ---------------------------------------------------------
	 ij_inter:
	 - calculates the pair potential as a function of dr^2.
	 ---------------------------------------------------------*/
	protected double ij_inter(double dr2, int i, int j, int il) {
		double invr;
		
		/*    if ( dr2 < pot_cut_off2 && dr2 > 0.0 ) {*/
		if ( dr2 > 0.0 ) {
			invr = diameter2/dr2;
			return( lj_eta*(Math.pow(invr,6.0) - Math.pow(invr,3.0)) );
		} else {
			return(0.0);
		}
		
	}
	
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
	protected double ij_inter_poly(double dr2, int i, int j, int il) {
		throw new LSSimException(LSSimException.PARAM,"The LJLS class does not support polydisperse calculations at present!");
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
	
	/** ---------------------------------------------------------
	 calc_e_from_scratch:
	 ---------------------------------------------------------*/
	protected double calc_e_from_scratch(int il) {
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
	
	/** ---------------------------------------------------------
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
					pot = ij_inter_act(dr2,i,j,il);
					newe += pot;
				}
			}
		}
		pot_cut_off = store_cut;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		
		return(newe);
	}
	
	/** ---------------------------------------------------------
	 calc_e_order_n:
	 This is invoked by calc_e_v_dens, and is not in use at present.
	 ---------------------------------------------------------*/
	private double calc_e_order_n(int il) {
		int j;
		double dr2,newe,store_cut,pot;
		
		store_cut = pot_cut_off;
		pot_cut_off = 1.0e20;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		newe = 0.0;
		for(j=0; j<n; j++) {
			if ( j != 0 ) {
				dr2 = ij_sep2(0,j,il,CUR_POS);
				pot = ij_inter_act(dr2,0,j,il);
				newe += pot;
			}
		}
		pot_cut_off = store_cut;
		pot_cut_off2 = pot_cut_off*pot_cut_off;
		
		return(newe);
	}
	
	/** ---------------------------------------------------------
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
	
	
	protected void report_simulation(int sweep) {
		// Perform general reporting:
		super.report_simulation(sweep);
		// Perform system-specific reporting:
		if( VIRT_NSWC == 1 ) {
			virtual_n_squared();
		}
	}
	
	/** ---------------------------------------------------------
	 ---------------------------------------------------------*/
	protected void check_the_system(int noise) {
		// Perform general checks:
		super.check_the_system(noise);
		// LJ-specific checks:
		double e_bits_check = lj_eta*(e12 - e6);
		if ( e_check > 0.0 && Math.abs((e_bits_check-e_check)/e_check) > E_TOL ) {
			throw( new LSSimException(LSSimException.CHECK,
					"fatal error: (e_bits_check="+e_bits_check+") != (e_check="+e_check+") by "+
					Math.abs(e_bits_check-e_check)));
		}
		
	}
	
	
	
	/** ---------------------------------------------------------
	 calc_gs_by_struct()
	 ---------------------------------------------------------*/
	protected double calc_gs_by_struct(double densi, int latti ) {
		
		if ( N_SWITCH == 0 ) {
			return( n*beta*LJLS.calc_Egs_of_struct(densi, latti) );
		} else {
			// fudge! to get the e_shift:
			double e_shift = 0.0;
			if( latti == 1 ) e_shift = EN_SHIFT;
			return( n*beta*LJLS.calc_Egs_of_struct(densi, INIT_NPHASE ) + e_shift );
		}
	}
	
	
	protected double calc_e_volume_move( int t_lat, double oldLinear, double nuLinear, double oldVol, double nuVol ) {
		if( npt_fixed_aspect ) {
			if( t_lat == c_lat ) {
				return lj_eta*(e12*Math.pow(nuLinear/oldLinear,-12.0)-e6*Math.pow(nuLinear/oldLinear,-6.0));
			} else {
				return lj_eta*(em12*Math.pow(nuLinear/oldLinear,-12.0)-em6*Math.pow(nuLinear/oldLinear,-6.0));
			}
		} else {
			return calc_e_from_scratch(t_lat);
		}
	}
	
	protected void accepted_volume_move( double oldLinear, double nuLinear, double oldVol, double nuVol ) {
		e6 = e6*Math.pow(nuLinear/oldLinear,-6.0);
		e12 = e12*Math.pow(nuVol/oldVol,-4.0);
		em6 = em6*Math.pow(nuVol/oldVol,-2.0);
		em12 = em12*Math.pow(nuVol/oldVol,-4.0);		
	}
	
	
	/* Particle move */
	
	double[] leold12 = new double[2];
	double[] lenew12 = new double[2];
	double[] leold6 = new double[2];
	double[] lenew6 = new double[2];
	
	public double calc_dE_sphere_move(int ip, int t_lat ) {
		/* - calc overlaps and dm */
		/*
		 lmold = calc_local_energy(im,m_lat,CUR_POS);
		 lmnew = calc_local_energy(im,m_lat,TRY_POS);
		 leold = calc_local_energy(im,c_lat,CUR_POS);
		 lenew = calc_local_energy(im,c_lat,TRY_POS);
		 */
		leold12[t_lat] = calc_local_energy_pow(ip,t_lat,CUR_POS,6.0);
		lenew12[t_lat] = calc_local_energy_pow(ip,t_lat,TRY_POS,6.0);
		//leold12 = calc_local_energy_pow(im,c_lat,CUR_POS,6.0);
		//lenew12 = calc_local_energy_pow(im,c_lat,TRY_POS,6.0);
		leold6[t_lat] = calc_local_energy_pow(ip,t_lat,CUR_POS,3.0);
		lenew6[t_lat] = calc_local_energy_pow(ip,t_lat,TRY_POS,3.0);
		//leold6 = calc_local_energy_pow(im,c_lat,CUR_POS,3.0);
		//lenew6 = calc_local_energy_pow(im,c_lat,TRY_POS,3.0);
		
		
		return ((lenew12[t_lat] - lenew6[t_lat])-(leold12[t_lat] - leold6[t_lat]));
	}
	
	public void accepted_sphere_move() {
		e12 = e12 + (lenew12[c_lat] - leold12[c_lat])/lj_eta;
		e6 = e6 + (lenew6[c_lat] - leold6[c_lat])/lj_eta;
		em12 = em12 + (lenew12[m_lat] - leold12[m_lat])/lj_eta;
		em6 = em6 + (lenew6[m_lat] - leold6[m_lat])/lj_eta;		
	}
	
	protected void accepted_basis_flip() {
		//swap_doubles(e6, em6);
		temp = e6; e6 = em6; em6 = temp;
		//swap_doubles(e12, em12);
		temp = e12; e12 = em12; em12 = temp;		
	}
	

	/** ---------------------------------------------------------
	 calc_e_v_dens():
	 - code to calc (ground-state?) energy as a function of density
	 Not referenced at present.
	 ---------------------------------------------------------*/
	private void calc_e_v_dens() {
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
	 calc_virial_pressure()
	 calculate and return the virial pressure:
	 ------------------------------------------------------- */
	protected double calc_virial_pressure() {
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
		vp = Math.sqrt(2.0)*calc_sys_density()/beta +
		Math.pow(diameter,3.0)*TotFdotR/(3.0*box[c_lat].x*box[c_lat].y*box[c_lat].z);
		
		/* check line */
		//System.out.printf(" VPi %f %f %f\n",calc_density()/beta, TotFdotR/(3.0*box.x*box.y*box.z),vp);
		
		return(vp);
	}
	
	/* -------------------------------------------------------
	 Converts to N^2 calc:
	 ------------------------------------------------------- */
	protected void make_it_order_n_squared( int latti ) {
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
	
	/** -------------------------------------------------------
	 Make the LS into a number of neighbours switch:
	 latt[il] and nns[il][i][inn] to be changed, using
	 NO_MORE_NNS as a end flag and checking that
	 MAX_NNS is large enough:
	 ------------------------------------------------------- */
	protected void alter_for_neighbour_switch() {
		int this_l, tother_l, ns_i, ns_j, ns_n;
		
		// Let the user know...
		System.out.println(" H THIS IS IN N_SWITCH (nn comparison) mode, for structure "+INIT_NPHASE);
		System.out.println(" H en_shift = "+EN_SHIFT);
		
		// Check that there is room at the inn:
		if( n > MAX_NNS ) {
			throw( new LSSimException(LSSimException.PARAM,
			"neighbour switch needs "+n+" neighbours in MAX_NNS!\n"));
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
	
	/** -------------------------------------------------------
	 Virtual O(N^2) calc:
//	Variable overwrites, memory leaks/overwrites.
	//
//	(int) based pbcs differ from +/- 1.0 (float or double).
	//
//	Are there any floats in the calculation?
//	EN_SHIFT adding in a float 0.0?
	//
	 ------------------------------------------------------- */	
	protected void virtual_n_squared() {
		int ns_i, ns_j, di, edge;
		int[] num_inters = new int[2];
		// NOTE: These were long-double, which is not accessible in java.
		// These may not have been doing anything on a given platform anyway.
		double eonn,ijsepsqrd, ije6, ije12, newe6, newe12;
		double denn, nnprb, ijlatsepsqrd, dl[] = new double[3], dp[] = new double[3], gstot;
		double truncEtot, fullEtot, dennOLD;
		// ANJ End of long-doubles.
		double max_circ_trunc2;
		double ijunscaledsepsqrd, unscaled_boxsqrd = 0.0;
		
		// Using a spherical truncation as big as the cell:
		// It implements a _fixed_length_ cut-off, and so the unscaled cell is used later:
		max_circ_trunc2 = Math.pow((double)lsize[0][2]*0.5*(Math.sqrt(2.0/3.0)),2.0);
		
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
						ijsepsqrd += dp[di]*dp[di]*box2[c_lat][di];
						ijlatsepsqrd += dl[di]*dl[di]*box2[c_lat][di];
						if( di == 0 ) unscaled_boxsqrd = init_boxes[c_lat].x*init_boxes[c_lat].x;
						if( di == 1 ) unscaled_boxsqrd = init_boxes[c_lat].y*init_boxes[c_lat].y;
						if( di == 2 ) unscaled_boxsqrd = init_boxes[c_lat].z*init_boxes[c_lat].z;
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
		if( !init_nncalc ) {
			System.out.println(" NNTEST "+Math.sqrt(max_circ_trunc2)+" "+ 
					2.0*num_inters[0]/216.0+" "+ 2.0*num_inters[1]/216.0);
			init_nncalc = true;
		}
		return;
	}
	
	/*--------------------------------
	 Static utility methods
	 --------------------------------*/
	
	/**--------------------------------------------------------
	 dEgs_fcchcp:
	 - return the static-lattice energy difference 
	 for a given density [fcc-(ideal)hcp].
	 --------------------------------------------------------*/
	static double dEgs_fcchcp(double dens) {
		double A12, A6, Uhcp, Ufcc;
		
		/* unit of distance is the 1st NN seperation */
		A12 = 12.132293768711;  A6  = 14.454897093822;
		Uhcp = 2.0*( Math.pow(dens,4.0)*A12 - Math.pow(dens,2.0)*A6 );
		
		A12 = 12.131880196191;  A6 = 14.453920885450;
		Ufcc = 2.0*( Math.pow(dens,4.0)*A12 - Math.pow(dens,2.0)*A6 );
		
		return(Ufcc-Uhcp);
	}
	
	/**--------------------------------------------------------
	 calc_Egs_of_struct:
	 - return the static-lattice energy 
	 for a given density and structure.
	 0 = hcp, 1 = fcc
	 --------------------------------------------------------*/
	static double calc_Egs_of_struct(double dens,int i_lat) {
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

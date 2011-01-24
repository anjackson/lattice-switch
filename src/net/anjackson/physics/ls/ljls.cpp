/*-------------------------------------------------------------------------------------
  Lennard-Jones Lattice-Switch Monte Carlo Code

  v1.1 : 19th August 1999 : Andrew N Jackson
  -------------------------------------------------------------------------------------
  Notes:
   - change to sure arrays instead of structs?
   - more dynamic allocation of multidimensional arrays?
   - change weight file format to XY.
   - change to NAG-style Nd arrays?
   - NPT to be implemented!
   - Load-save to be implemented!
  -------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "strongsamp.h"
#include "histog.h"
#include "nrutil.h"
#include "anjutil.h"

/*--------------------------------
  Global definitions and variables
  --------------------------------*/

/* "Hard-coded" system parameters */
#define NPT_SIM     1                              /* 0=NVT ensemble, 1=NPT ensemble */
#define LAT_SWITCH  1                              /* Lattice-Switch off/on 0/1 */
#define INIT_PHASE  1                              /* Initial phase hcp/fcc 0/1 */
                                                   /* Init NNs in NNmode, 0=some, 1=all */

#define VIRT_NSWC   0                              /* Use the virtual NN-switch move*/
#define N_SWITCH    0                              /* Turn LS into a neighbour-switch */
#define INIT_NPHASE 0                              /* init(+always) lattice-id for NN calc */
#define EN_SHIFT    0.0                           /* shifting the excitation energies */
                                                   /* 14.0 for 216, 80 for 1728 */

#define MOVE_TYPE   1                              /* 1 = RW, 2 = RW + cut-off */
#define FLIP_FREQ   0                              /* frequency of basis-flip attempts:
						      0 = once every move,
						      1 = <once every sweep>,
						      2 = as defined be bflipfreq */
#define MAX_NNS     50                            /* Max number of nearest neighbours */
#define E_TOL       1e-6                           /* Tolerance of energy checks */
#define NN_INC_TOL  1e-6                           /* NN inclusion tolerance */
#define EXPO_SHIFT  0.0                          /* exponential shift for no overflows */
#define SGL_PLOTS   0                              /* use SGL graphical output? */
#define CALC_COM    0                              /* Calc+ouput the CoM vector? */
#define CALC_VIRP   0                              /* Calculate the Virial pressure */
#define AUTO_STEPS  1                              /* Automatically tweak the stepsizes */
#define MAKEIT_NSQR 0                              /* Use O(N^2) calculation */

#if SGL_PLOTS == 1
#include "sgl.h"
#endif

/* useful static flag values */
#define NO_MORE_NNS -1                             /* nnlists: no more nns flag */
#define LOUDLY      0                              /* verbosity: how much output: */
#define QUIETLY     1
#define CUR_POS     0                              /* ij_sep2: calc which seperation: */
#define TRY_POS     1
#define ZERO_POS    2

/* User defined system parameters */
int     nx,ny,nz;                                 /* dimensions in nos of spheres */
double  Pres,BetaPres;                            /* pressure, Beta x pressure */
double  dVol;
double  density;                                  /* density */
double  temp;                                     /* inverse temperature */
double  dr;                                       /* random-walk length parameter */
double  pot_cut_off,pot_cut_off2;                 /* interaction cut-off distance */
double  nn_cut_off, max_dr2;                      /* NN inclusion distance */
int     cold_start;                               /* 1 = cold start, 0 = load old conf */
int     total_sweeps;                             /* total MCS to run */
int     equib_sweeps;                             /* equilibration period in MCS */
int     output_period;                            /* output period in MCS */
int     check_period;                             /* system checking period in MCS */
int     iseed;                                    /* random number generator seed */

/* System storage */
typedef struct { double x,y,z; } vect3d;          /* struct for vectors */
/* Pointers to be malloc'ed */
vect3d* latt[2];                                  /* pointer to 2*lattice array */
vect3d* disp;                                     /* pointer to displacements array */
vect3d tdisp;                                     /* trial displacements */
double* weights;                                  /* pointer to weight function */
int     ***nns;                                  /* nearest neighbour table */

/* General simulation variables */
int     n;                                        /* total number of spheres */
double  beta;                                     /* inverse temperature */
double  lj_eta;                                   /* lj leading (energy) coefficient */
double  e;                                        /* current energy of system */
double  dEgs;                                     /* ground state energy diff. */
double  em,m;                                     /* conjugate energy and order param. */
int     c_lat, m_lat;                             /* identity of the (c)urrent and
                                                     (m) order-param lattices [=1/2]*/
double  one_over_n;                               /* value of 1/n (speed-up calc) */
int     wlo, whi;                                 /* range of the MCMC weight func in m
						     [wlo_m,whi_m] */
int     wzero;
vect3d  init_box;                                 /* inital size of system box a,b,c */
vect3d  box; double box2[3];                      /* current system box a,b,c & ^2 */
double  obox[3];                                  /* old box-size, for trial moves */
double  drobox[3];                                /* RW step in each dirn */
double  diameter,diameter2;                       /* sphere diameter & diameter^2 */
int     max_inn[2];                               /* the observed max no. of NNs */
int     time_in_phase[2];                         /* stores time spent in each phase */
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
RanGen* rng;

/* automatic step-size tweaking parameters */
double DrARlo = 0.3, DrARid = 0.35, DrARhi = 0.4, DrARtweak = 0.05;
double DvARlo = 0.45, DvARid = 0.5, DvARhi = 0.55, DvARtweak = 0.05;

/* ---------------------------------------------------------
   ---------------------------------------------------------*/
inline void swap_doubles( double *aswap, double *bswap ) {
    double cswap;
    cswap = *aswap;
    *aswap = *bswap;
    *bswap = cswap;
}

/*----------------------
  Simulation subroutines
  ----------------------*/

/* ---------------------------------------------------------
   read_params_from_user:
   - read the simulation defining parameters from the user
   - NB file format is quite strict
   ---------------------------------------------------------*/
void read_params_from_user(void) {
    char dummy[256],yn[1];

    scanf("%[^\n]s",dummy);
    scanf("%i%i%i\n\n",&nx,&ny,&nz);
    n = nx*ny*nz;

    scanf("%[^\n]s",dummy);
    scanf("%lf\n\n",&density);
    //if NPT_SIM == 0
    diameter = pow(density,1.0/3.0);
    //else
    //diameter = 1.0;
    //endif
    diameter2 = diameter*diameter;

    scanf("%[^\n]s",dummy);
    scanf("%lf\n\n",&temp);
    beta = 1.0/temp;
    lj_eta = 4.0*beta;

    scanf("%[^\n]s",dummy);
    scanf("%lf\n\n",&dr);

#if NPT_SIM == 1
    scanf("%[^\n]s",dummy);
    scanf("%lf %lf\n\n",&Pres,&dVol);
    BetaPres = beta*Pres/density;
#else
    scanf("%[^\n]s",dummy);
    scanf("%lf\n\n",&pot_cut_off);
    pot_cut_off2 = pot_cut_off*pot_cut_off;
#endif

    scanf("%[^\n]s",dummy);
    scanf("%lf\n\n",&nn_cut_off);

    scanf("%[^\n]s",dummy);
    scanf("%1s\n\n",yn);
    if ( strncmp(yn,"y",1) == 0 ) cold_start = 0;
    if ( strncmp(yn,"n",1) == 0 ) cold_start = 1;

    scanf("%[^\n]s",dummy);
    scanf("%i\n\n",&total_sweeps);

    scanf("%[^\n]s",dummy);
    scanf("%i\n\n",&equib_sweeps);

    scanf("%[^\n]s",dummy);
    scanf("%i\n\n",&output_period);

    scanf("%[^\n]s",dummy);
    scanf("%i\n\n",&check_period);

    scanf("%[^\n]s",dummy);
    scanf("%i\n\n",&iseed);

    scanf("%[^\n]s",dummy);
    scanf("%s %lf %lf %i %i %i %i %i %lf %lf\n\n",
	  dummy,&ss_lo,&ss_hi,&ss_n,&ss_nd,&ss_lock,&samp_new_tpm,
	  &read_weights,&wf_clip_lo,&wf_clip_hi);
    if ( strcmp(dummy,"strong")==0 ) samp_type = STRONG_SAMP;
    if ( strcmp(dummy,"mcmc")==0 ) samp_type = MCMC_SAMP;
    if ( strcmp(dummy,"imp")==0 ) samp_type = IMP_SAMP;

}


/* ---------------------------------------------------------------------
   read_weights_file:
   - read-in or construct the multi-canonical weight function
   - if 'filename' exists, then read it, else use a flat weight function
   - also, only build to cope with N=216, N=1728 or N=5832
   - should change to XY format.
   ---------------------------------------------------------------------*/
void read_weights_file(char *filename) {
    int wi;
    FILE *fp;

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
	fprintf(stderr,"ERROR: read_weights_file cannot cope with n = %i\n",n);
	exit(EXIT_FAILURE);
    }
    /* Calculate where M=0 lies in the [0,whi] array*/
    wzero = whi/2;

    /* Allocate memory for the weights array */
    weights = (double *) calloc(whi-wlo+1,sizeof(double));

    /* If no weights file, define a flat weight function */
    if ( (fp = fopen(filename,"r")) == (FILE *)NULL ) {
	for ( wi = wlo; wi <= whi; wi++ ) {
	    weights[wi] = 1.0;
	}
    }
    /* But if there is a weights file, load it */
    else {
	for ( wi = wlo; wi <= whi; wi++ ) {
	    fscanf(fp,"%lf",&weights[wi]);
	}

    }
    /* close the weights file */
    fclose(fp);

}

/* ---------------------------------------------------------
   write_output_header:
   - write the definition of the simulation parameters to the
     standard output stream
   ---------------------------------------------------------*/

void write_output_header(void) {

    printf(" H Lennard-Jones lattice-switch monte carlo code: hcp v fcc:\n");
    printf(" H \n");
    printf(" H --Code definitions--\n");
    printf(" H ensemble = ");
    if ( NPT_SIM == 0 ) printf("NVT\n");
    if ( NPT_SIM == 1 ) printf("NPT\n");
    printf(" H sphere move mechanism = ");
    if ( MOVE_TYPE == 1 ) printf("random walk\n");
    if ( MOVE_TYPE == 2 ) printf("random walk + cut-off\n");
    if ( MOVE_TYPE == 3 ) printf("top-hat\n");
    printf(" H \n");
    printf(" H --User definitions--\n");
    printf(" H system of (%i,%i,%i) = %i spheres\n",nx,ny,nz,n);
    printf(" H density = %f\n",density);
    printf(" H  ...diameter = %f\n",diameter);
    printf(" H beta = %le\n",beta);
    printf(" H lj_eta = 4*beta = %le\n",lj_eta);
    printf(" H dr = %f\n",dr);
    if( NPT_SIM == 1 ) {
	printf(" H NPT: Pres = %g, dVol = %g\n",Pres,dVol);
    }
    printf(" H interaction cut off = %f\n",pot_cut_off);
    printf(" H nn cut off = %f\n",nn_cut_off);
    if ( cold_start == 0) printf(" H loading initial configuration from init.conf\n");
    if ( cold_start == 1) printf(" H cold start\n");
    printf(" H total_sweeps = %i\n",total_sweeps);
    printf(" H equib_sweeps = %i\n",equib_sweeps);
    printf(" H output_period = %i\n",output_period);
    printf(" H check_period = %i\n",check_period);
    printf(" H iseed = %i\n",iseed);
    printf(" H \n");
#if VIRT_NSWC == 1
    printf(" H VIRTUAL n-switch, structure %i, en_shift = %g\n",INIT_PHASE,EN_SHIFT);
#endif
    printf(" H --Notes--\n");
    printf(" H lat %i: up to %i NNs given a cut-off of %lf.\n",0,max_inn[0],nn_cut_off);
    printf(" H lat %i: up to %i NNs given a cut-off of %lf.\n",1,max_inn[1],nn_cut_off);
    printf(" H dEgs(fcc-hcp): = %le = %leN\n",dEgs,dEgs/((double)n));
    printf(" H There now follows an easily machine-readable list of parameter values\n");
    if ( NPT_SIM == 0 ) printf(" P NVT\n");
    if ( NPT_SIM == 1 ) printf(" P NPT\n");
    printf(" P %i %i %i %i\n",nx,ny,nz,n);
    if ( NPT_SIM == 0 ) printf(" P %.12le %.12le\n",density,beta);
    if ( NPT_SIM == 1 ) printf(" P %.12le %.12le %.12le\n",density,Pres,beta);
    printf(" P %.12le %.12le %.12le\n",0.0 ,0.0 , 0.0/((double)n));
    printf(" P %le %le %le\n",dr,pot_cut_off,nn_cut_off);
    printf(" P %i %i %i %i %i\n",cold_start,total_sweeps,
	   equib_sweeps,output_period,check_period);
    printf(" P %i\n",iseed);
    printf(" P %le %le %i %i %i %i %lf %lf\n",ss_lo,ss_hi,ss_n,ss_nd,ss_lock,
	   samp_type,wf_clip_lo,wf_clip_hi);
}

/*----------------------------------------------------------------
  make_lattices:
  - builds the fcc and hcp lattices, incorporating the LS mapping.
  - requires 6n planes in the z direction for normal pbcs.
  -----------------------------------------------------------------*/
void make_lattices(void) {
    int ilatt,ix,iy,iz,in;
    double xs,ys,zs,t;
    double xo,yo,zo,xoo,yoo;
    vect3d shifts[2];

    /* Allocate memory for the lattices and displacements */
    latt[0] = (vect3d *) calloc(n+1,sizeof(vect3d));
    latt[1] = (vect3d *) calloc(n+1,sizeof(vect3d));
    disp    = (vect3d *) calloc(n+1,sizeof(vect3d));

    /* Define offset in the x-dirn between A and B planes */
    t  = sqrt(3.0)/3.0;
    /* Define distances between neighbouring atoms in the same stacking plane */
    xs = sqrt(3.0)/2.0;
    ys = 1.0;
    zs = sqrt(2.0/3.0);

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
		latt[0][in].x=xoo;
		latt[0][in].y=yoo+ys*((double)(iy-1));
		latt[0][in].z=zo;
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
#if NPT_SIM == 1
    /* if this is an NPT simulation, resize the box to get the desired density */
    box.x = pow(density,-1.0/3.0)*init_box.x;
    box.y = pow(density,-1.0/3.0)*init_box.y;
    box.z = pow(density,-1.0/3.0)*init_box.z;
#else
    /* use standard box-size:*/
    box.x = init_box.x;
    box.y = init_box.y;
    box.z = init_box.z;
#endif
    /* store box side length ^ 2 */
    box2[0] = box.x*box.x;
    box2[1] = box.y*box.y;
    box2[2] = box.z*box.z;
    /* calc RW step in each dirn */
    drobox[0] = dr/box.x;
    drobox[1] = dr/box.y;
    drobox[2] = dr/box.z;

    /* check that the NN cut-off is actually inside the box */
    for( ix = 0; ix < 3; ix++ ) {
	if (sqrt(box2[ix]) < nn_cut_off) {
	    fprintf(stderr,"ERROR: This cut-off (%g) does not fit in box (%i,%g)\n"
		    ,nn_cut_off,ix,sqrt(box2[ix]));
	    exit(EXIT_FAILURE);
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
inline double ij_sep2(int i, int j,int il,int trial) {
    double d[3],idisp[3];

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
void make_nn_lists(void) {
    int il,i,j,inn,max_lnn;
    double dr2;

    /* allocate for the 3D NN array*/
    nns = imatrix3d(0,1,0,n,0,MAX_NNS);

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
		    if ( il == 0 && i == 0 ) printf(" NN %i %i %i %lf %lf\n",il,i,j,dr2,sqrt(dr2));
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
	printf("ljls: too many NNs; %i > MAX_NNS(%i); exiting...\n",max_lnn,MAX_NNS);
	exit(EXIT_FAILURE);
    }

}

/* ---------------------------------------------------------
   ---------------------------------------------------------*/
void load_config(char *filename) {
    printf("l_c:  %s\n",filename);
}

/* ---------------------------------------------------------
   - save complete system information using the given filename
   ---------------------------------------------------------*/
void save_config(char *filename) {
    printf("s_c:  %s\n",filename);
 /*
 fprintf(fp,"C %f %f %f\n",latt[c_lat][i].x,latt[c_lat][i].y,latt[c_lat][i].z);
 */
}

/* ---------------------------------------------------------
   save_config_xyz:
   - in xyz format for use with xmol
   ---------------------------------------------------------*/
void save_config_xyz(char *filename) {
    int i,ii,ij,ik;
    FILE *fp;

    /* open file and outpout number of atoms to come */
    fp = fopen(filename,"w");
    fprintf(fp," %i\n\n",n+8);

    /* loop over all particles and output as carbons */
    for ( i = 0; i < n; i++ ) {
	fprintf(fp,"C %f %f %f\n",(latt[c_lat][i].x+disp[i].x)*box.x,(latt[c_lat][i].y+disp[i].y)*box.y,(latt[c_lat][i].z+disp[i].z)*box.z);
    }

    /* Loop over the corners of the pbc box, outputting as hydrogens */
    for ( ii = -1; ii <= 1; ii = ii + 2 ) {
	for ( ij = -1; ij <= 1; ij = ij + 2 ) {
	    for ( ik = -1; ik <= 1; ik = ik + 2 ) {
		fprintf(fp,"H %f %f %f\n",box.x*(ii*0.5),box.y*(ij*0.5),box.z*(ik*0.5));
	    }
	}
    }

    /* close the file */
    fclose(fp);

}

/* ---------------------------------------------------------
   ij_inter:
   - calculates the pair potential as a function of dr^2.
   ---------------------------------------------------------*/
inline double ij_inter(double dr2) {
    double invr;

    /*    if ( dr2 < pot_cut_off2 && dr2 > 0.0 ) {*/
    if ( dr2 > 0.0 ) {
	invr = diameter2/dr2;
	return( lj_eta*(pow(invr,6.0) - pow(invr,3.0)) );
    } else {
	return(0.0);
    }

}

/* ---------------------------------------------------------
   ij_inter_pow:
   - calculates given power-law as a function of dr^2.
   ---------------------------------------------------------*/
inline double ij_inter_pow(double dr2, double rpow) {

    if ( dr2 > 0.0 ) {
	return( pow(diameter2/dr2,rpow) );
    } else {
	return(0.0);
    }

}

/* ---------------------------------------------------------
   calc_density:
   ---------------------------------------------------------*/
double calc_density(void) {
    return(n*pow(diameter,3.0)/(sqrt(2.0)*box.x*box.y*box.z));
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
void init_variables(void) {

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
    printf("TEST: %g\n",(em-e)/(n*beta));

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
    do {
	rng = new RanGen();
	if( iseed > 0 ) {
	    rng->setSeed(iseed);
	} else {
	    iseed = rng->getSeed();
	}
    } while( rng->get() < 0.0 );

}

/* ---------------------------------------------------------
   calc_local_energy:
   ---------------------------------------------------------*/
inline double calc_local_energy(int i, int il,int trial) {
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
inline double calc_local_energy_pow(int i, int il,int trial, double rpow) {
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
      printf(" C acc_rate = %.12f.\n",acc_moves/trial_moves);
    if( NPT_SIM == 1 && noise == LOUDLY && vol_trial > 0)
      printf(" C dV_acc_rate = %.12f.\n",vol_acc/vol_trial);


    /* Output the average energy and energy step */
    if (noise == LOUDLY) {
	if( num_e > 0) printf(" C <e> = %f.\n",ave_e/num_e);
	if( num_De > 0 ) printf(" C <De> = %f.\n",ave_De/num_De);
	if( num_es > 0 ) printf(" C <eS> = %f.\n",ave_es/num_es);
    }


    /* Check the energy running totals*/
    e_bits_check = lj_eta*(e12 - e6);
    e_check = calc_e_from_scratch(c_lat);
    em_check = calc_e_from_scratch(m_lat);
    if (noise == LOUDLY) {
	printf(" C e = %f, e_check = %f.\n",e,e_check);
	printf(" C m = %f, m_check = %f.\n",em,em_check);
    }

    if ( e_check > 0.0 && fabs((e_bits_check-e_check)/e_check) > E_TOL ) {
	printf("fatal error: (e_bits_check=%f) != (e_check=%f) by %e.\n",e_bits_check,e_check,fabs(e_bits_check-e_check));
	exit(EXIT_FAILURE);
    }

    if ( e_check > 0.0 && fabs((e-e_check)/e_check) > E_TOL ) {
	printf("fatal error: (e=%f) != (e_check=%f) by %e.\n",e,e_check,fabs(e-e_check));
	exit(EXIT_FAILURE);
    } else {
	e = e_check;
    }

    if ( em_check > 0.0 && fabs((em-em_check)/em_check) > E_TOL ) {
	printf("fatal error: (em=%f) != (em_check=%f) by %e.\n",em,em_check,fabs(em-em_check));
	exit(EXIT_FAILURE);
    } else {
	em = em_check;
    }

}



/* ---------------------------------------------------------
   calc_gs_by_struct()
   ---------------------------------------------------------*/
inline double calc_gs_by_struct(double densi, int latti ) {

#if N_SWITCH == 0
   return( n*beta*calc_Egs_of_struct(densi, latti) );
#else
    // fudge! to get the e_shift:
    if( latti == 1 ) {
	e_shift = EN_SHIFT;
    } else {
	e_shift = 0.0;
    }
    return( n*beta*calc_Egs_of_struct(densi, INIT_NPHASE ) + e_shift );
#endif
}


/* ---------------------------------------------------------
   mc_volume_move():
   - try to perform a volume move, fixed aspect only
   ---------------------------------------------------------*/
inline void mc_volume_move(void) {
    double old_m, new_m, old_e, new_e, new_em;
    double nuVol, nuLinear, oldVol, cond, oldLinear;
    double old_ge, new_ge;
    int ss_decision;

    /* Only execute this routine one nth of the time */
    if ( rng->get() >= one_over_n ) return;
    vol_trial+=1.0;

    /* store the old aspect ratio */
    obox[0] = box.x;
    obox[1] = box.y;
    obox[2] = box.z;
    oldVol = box.x*box.y*box.z;
    oldLinear = pow(oldVol/(init_box.x*init_box.y*init_box.z),1.0/3.0);
    /* calc/store old orderparameter/energies */
    old_e = e;
    old_m = (2*c_lat-1)*(em-e);
    old_ge = calc_gs_by_struct(calc_density(), c_lat);

    /* create a trial volume move uniform in V */
    nuVol = oldVol + dVol*( (rng->get()) - 0.5 );
    nuLinear = pow(nuVol/(init_box.x*init_box.y*init_box.z),1.0/3.0);
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
    //printf("ch: %f %f %f, %f %f %f\n",old_e,e12-e6,lj_eta*(e12-e6),old_m,em12,em6);
    //printf("v: %f %f (%f)\n",oldVol,nuVol,nuLinear/oldLinear);
    new_e = lj_eta*(e12*pow(nuLinear/oldLinear,-12.0)-e6*pow(nuLinear/oldLinear,-6.0));
    new_em = lj_eta*(em12*pow(nuLinear/oldLinear,-12.0)-em6*pow(nuLinear/oldLinear,-6.0));
    //new_e = calc_e_from_scratch(c_lat);
    //new_em = calc_e_from_scratch(m_lat);
    new_m = (2*c_lat-1)*(new_em-new_e);
    new_ge = calc_gs_by_struct(calc_density(), c_lat);

    /* calc the acceptance probability */
    cond = -((new_ge-old_ge)+(new_e-old_e)+BetaPres*(nuVol-oldVol)-n*log(nuVol/oldVol));

    /* accept? - update energies &c */
    ss_decision = ss_control( old_m, new_m, cond );
    if( ss_decision == 1 ) {
	/* accept */
	e = new_e;
	e6 = e6*pow(nuLinear/oldLinear,-6.0);
	e12 = e12*pow(nuVol/oldVol,-4.0);
	em = new_em;
	em6 = em6*pow(nuVol/oldVol,-2.0);
	em12 = em12*pow(nuVol/oldVol,-4.0);
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
    printf("%i %g: %g: (%g) %g %g: %g %g: %g %g\n",
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
    printf("V(%f): e %f %f (%f) em %f %f (%f)\n",vol_trial,e12,e6,e12-e6,em12,em6,em12-em6);
    check_the_system(LOUDLY);
    */

}

/* ---------------------------------------------------------
   mc_sphere_move:
   ---------------------------------------------------------*/
inline void mc_sphere_move(int i) {
    double cond;
    double lmold,lmnew,lenew,leold,dE,dEM,old_m,new_m;
    int im, ss_decision;
    double lmold6,lmnew6,lenew6,leold6,lmold12,lmnew12,lenew12,leold12;

    /* - select sphere at random */
    //do {
    im = (int)( ((double)(rng->get())) * ((double)n) );
	//} while( im < 0 );
    if( im < 0 ) {
	printf("ERK! Chosen particle no. %i < 0!!!!\n",im);
	exit(EXIT_FAILURE);
    }

    /* - generate trial move - currently RW only */
    tdisp.x = disp[im].x + (rng->get()-0.5)*drobox[0];
    tdisp.y = disp[im].y + (rng->get()-0.5)*drobox[1];
    tdisp.z = disp[im].z + (rng->get()-0.5)*drobox[2];
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
    ss_decision = ss_control( old_m, new_m, cond );
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
	ave_De += fabs(dE);
    } else if( ss_decision == -1 ) {
	printf("mc_move_sphere: Move rejected!\n");
    }
    num_e++;
    ave_e += fabs(e);
}

/* ---------------------------------------------------------
   mc_basis_flip:
   - flip to other basis, if we can:
   ---------------------------------------------------------*/
inline void mc_basis_flip(void) {
    double cond, old_m;
    int ss_decision;

#if FLIP_FREQ == 1
    /* Only execute this routine one nth of the time */
    if ( rng->get() >= one_over_n ) return;
#elif FLIP_FREQ == 2
    /* Only execute this routine bflipfreq of the time */
    if ( rng->get() >= bflipfreq ) return;
#endif

    cond = -(em-e);
    old_m = (2*c_lat-1)*(em-e);

#if NPT_SIM == 1
    /* include the gs energy diff if in the NPT ensemble */
    cond += -(calc_gs_by_struct(calc_density(), m_lat)
		     - calc_gs_by_struct(calc_density(), c_lat));
#endif

    /* pass this to the SS routines */
    ss_decision = ss_control( old_m, old_m, cond );
    if( ss_decision == 1 ) {
	/* accept */
	c_lat = m_lat;
	m_lat = 1 - c_lat;
	swap_doubles(&e, &em);
	swap_doubles(&e6, &em6);
	swap_doubles(&e12, &em12);
    } else if( ss_decision == -1 ) {
	/* reject badly! */
	printf("mc_basis_flip: Move rejected, not enuff diagonal space!\n");
	printf("%i-%i: %g %g (%g)\n",c_lat,m_lat,em,e,em-e);
    }

}

/* ---------------------------------------------------------
   ---------------------------------------------------------*/
void output_sim_info_tail(void) {
    double dfe;

    /* acceptance rate */
    if( trial_moves > 0 )
      printf(" C dr = %f :: acc_rate = %.12f.\n",dr,acc_moves/trial_moves);

    /* average energy and energy step */
    if( num_e != 0.0 )  printf(" C <e> = %f.\n",ave_e/num_e);
    if( num_De != 0.0 ) printf(" C <De> = %f.\n",ave_De/num_De);
    if( num_es != 0.0 ) printf(" C <eS> = %f.\n",ave_es/num_es);

    /* Simply free-energy difference calc: */
    if ( time_in_phase[0] > 0 && time_in_phase[1] > 0 ) {
	printf(" Time spent in hcp = %i\n",time_in_phase[0]);
	printf(" Time spent in fcc = %i\n",time_in_phase[1]);
	dfe = log((double)time_in_phase[1]/(double)time_in_phase[0])/(double)n;
	printf(" T Dfe = %f\n",dfe);
    }
}

/* ---------------------------------------------------------
   output_nn_dist:
   ---------------------------------------------------------*/
void output_nn_dist(void) {
    int i,max_steps = 250,t_nn[2][250];

    for (i=0; i<max_steps; i++) {
	nn_cut_off = 0.999 + 2.0*(double)i/(double)max_steps;
	make_nn_lists();
	t_nn[0][i] = max_inn[0]; t_nn[1][i] = max_inn[1];
	/*	printf(" %f %i %i\n",nn_cut_off,max_inn[0],max_inn[1]);*/
	if (i>0) printf(" %f %i %i\n",nn_cut_off,t_nn[0][i]-t_nn[0][i-1],t_nn[1][i]-t_nn[1][i-1]);
    }

    exit(0);
}

/* ---------------------------------------------------------
   calc_e_v_dens():
    - code to calc energy as a function of density
   ---------------------------------------------------------*/
void calc_e_v_dens() {
    int i,il,max_steps = 10;
    double lel[2],onel[2],store_dens;
    FILE *fp;

    store_dens = density;

    fp = fopen("Ediff.v.dens.out","w");
    fprintf(fp,"# density Elhcp Elfcc Ehcp Efcc dEl dE\n");

    for (i=0; i<max_steps; i++) {
	density= 2.0*(double)(i+1)/(double)max_steps;
	diameter = pow(density,1.0/3.0);
	diameter2 = diameter*diameter;
	fprintf(fp,"%e ",density);
	for (il=0; il<2; il++) {
	    lel[il] = calc_local_energy(0,il,CUR_POS);
	    onel[il] = calc_e_order_n(il);
	}
	fprintf(fp,"%e %e %e %e %e %e\n",lel[0],lel[1],onel[0],onel[1],lel[1]-lel[0],onel[1]-onel[0]);
    }
    fclose(fp);

    density= store_dens;
    diameter = pow(density,1.0/3.0);
    diameter2 = diameter*diameter;
}


/* -------------------------------------------------------
   calc_com_vect( vect3d comvect ):
   calculate the com vector and put it in the given array:
   ------------------------------------------------------- */
void calc_com_vect( vect3d* comvect) {
    int i;

    comvect->x = 0.0;
    comvect->y = 0.0;
    comvect->z = 0.0;

    for( i=0; i<n; i++ ) {
	/*
	comvect->x += latt[c_lat][i].x + disp[i].x;
	comvect->y += latt[c_lat][i].y + disp[i].y;
	comvect->z += latt[c_lat][i].z + disp[i].z;
	*/
	comvect->x += disp[i].x;
	comvect->y += disp[i].y;
	comvect->z += disp[i].z;
    }

    comvect->x /= (double)n;
    comvect->y /= (double)n;
    comvect->z /= (double)n;

}

/* -------------------------------------------------------
   calc_virial_pressure()
   calculate and return the virial pressure:
   ------------------------------------------------------- */
double calc_virial_pressure( void ) {
    int i,j,iTot;
    double dr2,fij,rij,TotFdotR;
    double vp;

    /* Calculate sum(f.r): */
    TotFdotR = 0.0; iTot = 0;
    for(i=0; i<n; i++) {
	for(j=0; j<n; j++) {
	    if ( i != j ) {
		dr2 = ij_sep2(i,j,c_lat,CUR_POS);
		rij = sqrt(dr2);
		fij = -12.0*4.0*pow(diameter,12.0)/pow(rij,13.0)
		       +6.0*4.0*pow(diameter,6.0)/pow(rij,7.0);
		TotFdotR += fij*rij;
		iTot++;
	    }
	}
    }
    //TotFdotR /= (double)iTot;
    TotFdotR /= -2.0;

    /* Calc total virial: */
    vp = sqrt(2.0)*calc_density()/beta +
      pow(diameter,3.0)*TotFdotR/(3.0*box.x*box.y*box.z);

    /* check line */
    //printf(" VPi %f %f %f\n",calc_density()/beta, TotFdotR/(3.0*box.x*box.y*box.z),vp);

    return(vp);
}

/* -------------------------------------------------------
   tweak_for_acc_rates
   Changes the step-sizes to make the acc_rates right.
   ------------------------------------------------------- */
void tweak_for_acc_rates( void ) {
    double vAR, rAR;

    /* If position acc_rate is out of bounds, tweak it */
    rAR = ((double)acc_moves)/((double)trial_moves);
    if( rAR < DrARlo || rAR > DrARhi ) {
	if( rAR > 0.0 && rAR < 1.0 ) {
	    dr = dr*log(DrARid)/log(rAR);
	} else {
	    if( rAR < DrARlo ) dr = dr*(1.0 + DrARtweak);
	    if( rAR > DrARhi ) dr = dr*(1.0 - DrARtweak);
	}
    }
    /* calc RW step in each dirn */
    drobox[0] = dr/box.x;
    drobox[1] = dr/box.y;
    drobox[2] = dr/box.z;

#if NPT_SIM == 1
    /* If volume acc_rate is out of bounds, tweak it */
    vAR = ((double)vol_acc)/((double)vol_trial);
    if( vAR < DvARlo || vAR > DvARhi ) {
	if( vAR > 0.0 && vAR < 1.0 ) {
	    dVol = dVol*log(DvARid)/log(vAR);
	} else {
	    if( vAR < DvARlo ) dVol = dVol*(1.0 + DvARtweak);
	    if( vAR > DvARhi ) dVol = dVol*(1.0 - DvARtweak);
	}
    }
    /* Notify */
    printf(" AR dr = %g dVol = %g\n", dr, dVol);
#else
    /* Notify */
    printf(" AR dr = %g\n", dr, dVol);
#endif

}


/* -------------------------------------------------------
   Converts to N^2 calc:
    - THIS MUST BE CHECKED!!!!!
   ------------------------------------------------------- */
void write_out_nn_list(int ill, char* fname) {
    int nni, nnj, nnn;
    FILE* nnfp;

    nnfp = fopen(fname,"w");
    for( nni = 0; nni < n; nni++ ) {
	nnn = 0; nnj = nns[ill][nni][nnn];
	while( nnj != NO_MORE_NNS ) {
	    fprintf(nnfp,"%i %i\n",nni,nnj);
	    nnn++; nnj = nns[ill][nni][nnn];
	}
    }
    fclose(nnfp);
}

/* -------------------------------------------------------
   Converts to N^2 calc:
    - THIS MUST BE CHECKED!!!!!
   ------------------------------------------------------- */
void make_it_order_n_squared( int latti ) {
    int ns_i, ns_j, ns_n;

    printf(" H Making lattice %i use a O(N^2) neighbour list...\n",latti);

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
void alter_for_neighbour_switch(void) {
    int this_l, tother_l, ns_i, ns_j, ns_n;

    // Let the user know...
    printf(" H THIS IS IN N_SWITCH (nn comparison) mode, for structure %i.\n",INIT_NPHASE);
    printf(" H en_shift = %g.\n",EN_SHIFT);

    // Check that there is room at the inn:
    if( n > MAX_NNS ) {
	fprintf(stderr,"ARGH! I need %i neighbours in MAX_NNS!\n",n);
	exit(EXIT_FAILURE);
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

//Variable overwrites, memory leaks/overwrites.
//
//(int) based pbcs differ from +/- 1.0 (float or double).
//
//Are there any floats in the calculation?
// EN_SHIFT adding in a float 0.0?
//
//Note, long-double may not actually be doing anything!
//

void virtual_n_squared( void ) {
    int ns_i, ns_j, di, edge;
    int num_inters[2];
    long double eonn,ijsepsqrd, ije6, ije12, newe6, newe12;
    long double denn, nnprb, ijlatsepsqrd, dl[3], dp[3], gstot;
    long double truncEtot, fullEtot, dennOLD;
    double max_circ_trunc2;
    static int init_nncalc = 0;
    double ijunscaledsepsqrd, unscaled_boxsqrd;

    // Using a spherical truncation as big as the cell:
    // It implements a _fixed_length_ cut-off, and so the unscaled cell is used later:
    max_circ_trunc2 = pow((double)nz*0.5*(sqrt(2.0/3.0)),2.0);

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
		    if( fabs(fabs(dl[di])-0.5) < NN_INC_TOL ) {
		      edge = 1;
		    }
		}
		if( edge == 0 ) {
		  // particle-particle interaction:
		  ije6 = ij_inter_pow(ijsepsqrd,3.0); ije12 = ij_inter_pow(ijsepsqrd,6.0);
		  newe6 += ije6; newe12 += ije12;
		//printf("ES %i %i %g  %g %g %g\n",ns_i,ns_j,sqrt(ijsepsqrd),ije6,ije12,lj_eta*(ije12-ije6));
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
		//printf("GS %i %i %g  %g %g %g\n",ns_i,ns_j,sqrt(ijsepsqrd),ije6,ije12,lj_eta*(ije12-ije6));
		  newe6 -= ije6; newe12 -= ije12;
		  if( ijunscaledsepsqrd < max_circ_trunc2 ) fullEtot -= lj_eta*(ije12 - ije6);
		  if( ijunscaledsepsqrd < max_dr2 ) truncEtot -= lj_eta*(ije12 - ije6);

		  gstot += 4.0*(ije12-ije6);
		//printf("TOT %i,%i %e %e %e\n",ns_i,ns_j,newe6,newe12,newe12-newe6);
		}
	    }
	}
    }

    eonn = lj_eta*(newe12 - newe6) + EN_SHIFT;
//dennOLD = (eonn - calc_e_from_scratch(c_lat) );
    dennOLD = (eonn - truncEtot);
    denn = fullEtot - truncEtot;
    nnprb = exp(-denn);
    printf(" NN2 %.20Le %.20Le %.20Le %.20Le %.20Le\n", dennOLD, denn, nnprb, truncEtot, fullEtot);
    if( init_nncalc == 0 ) {
	printf(" NNTEST %.20e %.12e %.12e\n",sqrt(max_circ_trunc2), 
	       2.0*num_inters[0]/216.0, 2.0*num_inters[1]/216.0);
	init_nncalc = 1;
    }
    return;
}

/*-------------------------------------------------------------------------------------
   - Main program
  -------------------------------------------------------------------------------------*/
int main( void ) {
    int sweep,i;
    vect3d comv;

    /* Initialising */
    read_params_from_user();
    make_lattices();
#if SGL_PLOTS == 1
    sgl_init_viewer(800, 300);
    sgl_init_2d_plot(0.0, 0.5, 10000.0, 2.0);
    sgl_bg_col(0.0, 0.0, 0.0);
    sgl_clear();
    sgl_update();
    sgl_fg_col(1.0, 1.0, 1.0);
#endif
    comv.x = 0.0; comv.y = 0.0; comv.z = 0.0;
    /*output_nn_dist();*/
    make_nn_lists();
#if N_SWITCH == 1
    // Change the hcp-fcc LS into a single-phase NN-switch:
    alter_for_neighbour_switch();
#endif
#if MAKEIT_NSQR == 1
    make_it_order_n_squared(0);
    make_it_order_n_squared(1);
#endif
    write_out_nn_list(0,"nns.0.dat");
    write_out_nn_list(1,"nns.1.dat");
    /*calc_e_v_dens();*/
    init_variables();
    /*read_weights_file("weights");*/
    if ( cold_start == 0 ) load_config("init.conf");
    write_output_header();
    /* sampling init */
    ss_init_ss(ss_lo,ss_hi,ss_n,ss_nd,ss_lock,samp_type,rng);
    ss_set_exp_shift(EXPO_SHIFT);
    ss_load_tpm("TPM.in"); // load the old TPM in from a standard file:
    // generate a wf from the tpm:
    if ( samp_type == MCMC_SAMP ) {
	ss_mcmc_weights_from_tpm();
	if ( read_weights == 1 ) ss_load_mcmc_weights("weights.in");
	ss_clip_mcmc_weights(wf_clip_lo,wf_clip_hi);
	ss_output_mcmc_weights("mcmc.weights.dat");
    }

    // reset the TPM to zero, if required:
    if( samp_new_tpm == 1 ) {
	ss_set_tpm_to_zero();
    }
    ss_print_sampling_info();

    printf(" Einit %e %i %e %e\n",e,c_lat,em,(double)(2*c_lat-1)*(em-e));

    for ( sweep = 1; sweep <= total_sweeps; sweep++ ) {        /* Main loop... */
	/* Count the time spent in each: */
	time_in_phase[c_lat]++;

	/* Periodically output the system status */
	if ( sweep%output_period == 0 && sweep > equib_sweeps ) {
	    num_es++; ave_es += em-e;
	    printf(" E %i %e %i %e\n",sweep,e,c_lat,(double)(2*c_lat-1)*(em-e));
#if NPT_SIM == 1
	    /*
	    printf(" V %i %f %f %f %f\n",sweep,box.x*box.y*box.z,box.x,box.y,box.z);
	    */
	    printf(" V %i %f %f\n",sweep,box.x*box.y*box.z,calc_density());
#endif
#if SGL_PLOTS == 1
	    sgl_plot_line((double) sweep, 1.0, (double)(sweep+1), calc_density() );
	    //sgl_plot_fcircle((double)(sweep+1), calc_density(), 0.1);
	    sgl_update();
#endif
#if CALC_COM == 1
	    calc_com_vect(&comv);
	    printf(" CoM %e %e %e\n",comv.x,comv.y,comv.z);
#endif
#if CALC_VIRP == 1
	    printf(" VP %f\n", calc_virial_pressure() );
#endif
#if VIRT_NSWC == 1
	    virtual_n_squared();
#endif
	}
	/* Periodically check the system */
	if ( sweep%check_period == 0) {
	    save_config("midsim.conf");
	    save_config_xyz("midsim.xyz");
	    ss_output_tpm("TPM.dat");
	    ss_output_tpm_pdf("pdf.tpm.dat");
	    ss_output_vs_pdf("pdf.vs.dat");
	    check_the_system(LOUDLY);
#if AUTO_STEPS == 1
	    tweak_for_acc_rates();
#endif
	    /* re-zero the acceptance rate counters */
	    acc_moves = 0; trial_moves = 0;
#if NPT_SIM == 1
	    vol_acc = 0; vol_trial = 0;
#endif
#if N_SWITCH == 1
	    if( time_in_phase[1] > 0 && time_in_phase[0] > 0 )
	      printf(" NN %i %i %g %g\n", time_in_phase[0], time_in_phase[1], (double)time_in_phase[1]/(double)time_in_phase[0], log((double)time_in_phase[1]/(double)time_in_phase[0])/(double)n);
#endif
	}

	for ( i = 0; i < n; i++ ) { 	                       /* Sphere loop... */
	    //printf("e %f %f (%f) em %f %f (%f)\n",e12,e6,e12-e6,em12,em6,em12-em6);
	    //check_the_system(LOUDLY);

#if NPT_SIM == 1
	    mc_volume_move();                                  /* Volume move */
#endif

	    mc_sphere_move(i);	                               /* Sphere move */

#if LAT_SWITCH == 1
	    mc_basis_flip();	                               /* Basis-flip move */
#endif

	}

    }

    /* Finish off */
    em = calc_e_from_scratch(m_lat);
    check_the_system(LOUDLY);
    save_config("final.conf");
    save_config_xyz("final.xyz");
    output_sim_info_tail();

    exit(EXIT_SUCCESS);
}




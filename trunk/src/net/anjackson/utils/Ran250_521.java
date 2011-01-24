/**-------------------------------------------------------------
 * jBinLats - Ran250_521.java
 * net.anjackson.utils.Ran250_521
 * 
 * Created on 19-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.utils;

import java.util.Date;

import edu.cornell.lassp.houle.RngPack.RandomSeedable;

/** RAN250/521, a.k.a RANVEC.
 * Vectorized version of a (hopefully) improved shift-register    
 * pseudorandom number generator.                                 
 * Based on two independent generators with different pairs       
 * of "magic numbers", here chosen as (250,103) and (521,168).    
 * The output of those is XORed together and yields the final     
 * random number. For theoretical reasons, this generator should  
 * be significantly less hampered by correlations than the simple 
 * good old R250.                                                 
 *                                                                
 * Disclaimer of warranty:                                        
 * I will not be held responsible for any problems whatsoever     
 * which this program may cause.                                  
 *                                                                
 * B. Duenweg, July 9, 1996.                                      
 *                                                                
 * Some (not all!) other possible magic umbers are (see also      
 * Kirkpatrick and Stoll, Journ. Comp. Phys. 40, 517 (1981), and  
 * N. Zierler, Information and Control 15, p. 67 (1969)):         
 * (98, 27)                                                       
 * (521, 32)                                                      
 * (521, 48)                                                      
 * (521, 158)                                                     
 * (607, 105)                                                     
 * (607, 147)                                                     
 * (607, 273)                                                     
 * (1279, 216)                                                    
 * (1279, 418)                                                    
 * (2281, 715)                                                    
 * (2281, 915)                                                    
 * (2281, 1029)                                                   
 * (9689, 4187)                                                   
 *                                                                
 * The program can be changed very easily to other magic numbers  
 * by just changing the parameters                                
 * BIGMAGIC1, BIGMAGIC2, SMALLMAGIC1 and SMALLMAGIC2.             
 * Since it is explicitly written for 31 bit integers, it         
 * should produce the SAME sequence on any architecture.          
 *                                                                
 * This is how it works:                                          
 *                                                                
 * - Specify a "seed" value iseed.                                
 *                                                                
 * - Specify how many random numbers should be generated          
 *   on ONE call of the generator. Call this number, say,         
 *   nrand. A typical value is, say, 100000.                      
 *                                                                
 * - Invoke the function                                          
 *   init_vector_random_generator(iseed,nrand).                   
 *   This does the following:                                     
 *   - A simple congruential generator is run for NWARM           
 *     (in our case 10000) times.                                 
 *   - Two integer working arrays of size BIGMAGIC1 + nrand,      
 *     BIGMAGIC2 + nrand, named rand_w_array1 and rand_w_array2   
 *     are created (i.e. memory is allocated).                    
 *   - The first BIGMAGICX array elements are                     
 *     filled with (congruential) random numbers.                 
 *   - The bit columns are treated for linear independence.       
 *   - The shift-register generators are run for nrand times      
 *     (in order to also warm them up) and are then ready to go.  
 *                                                                
 * - Invoke the function                                          
 *   vector_random_generator(nrand, random_numbers).              
 *   Result: Normalized double precision random numbers in [0:1]  
 *   are written on the array random_numbers.                     
 *   The current "state" of the generator is "coded" in the first 
 *   BIGMAGICX elements of rand_w_arrayX.                         
 *                                                                
 * - There are also routines provided to save this status to a    
 *   file, and read it from there. Have care: In case you read    
 *   the status, make sure to first run the function              
 *   init_vector_random_generator in order to assure proper       
 *   memory allocation!                                           
 *
 * N.B. Original C version is available from the author's homepage:
 * http://www.mpip-mainz.mpg.de/~duenweg/Public/preprints.html
 *
 * @author ajackso1
 * @version $Id: Ran250_521.java 605 2006-02-27 17:35:45Z anj $
 *
 */
public class Ran250_521 extends RandomSeedable {

	/* Anj's notes: 
	   - The #pragma ivdep directive enables loop splitting,
	   vectorization, cache_bypass, and other optimizations. It's a CRAY
	   T3E thing.  I've commented them out for now.
	   - The routine could occasionally returned 1.0, because it used to use:
	   #define FACTOR 4.6566128752457969e-10, which is 1/BIGINT.  I've changed
	   it to 1/(BIGINT+1), so now it works properly.
	   - However, it is possible that the routine now returns neither 0 or 1,
	   which would be wrong, but not tragically so.
	   - Also, I'm adding an alternative RNG, SNWS.(14/3/2000).
	   */

	/*                                                                */	
	
	static int DEFAULT_SIZE =10000;
//	 definitionings for RANMAR:
	static int BIGMAGIC1 =250;                       /* magic numbers for the      */
	static int SMALLMAGIC1 =103;                     /*   first generator          */
	static int BIGMAGIC2 =521;                       /* magic numbers for the      */
	static int SMALLMAGIC2 =168;                     /*   second generator         */
	static int NBIT =32;                             /* use only (NBIT - 1) bits   */
	static int BIGINTEGER =2147483647;               /* = largest integer          */
	static float BIGFLOAT =2147483647.f;                /* same in float              */
	static float FACTOR =4.6566128730773925781e-10f;    /* = 1. / (largest integer+1)   */
	static float MULTIPLY =16807.f;                     /* for congruential generator */
	static int NWARM =10000;                         /* number of empty runs, c.g. */
	
	double random_numbers[];
	int    rand_w_array1[], rand_w_array2[];
	int    bign,snws_count;
	
    int total = DEFAULT_SIZE;
    int entry = total;
    long max_seed = BIGINTEGER;
    long seed = -1;
    
    /* Constructors ---------------------------------------------------*/
	
    public Ran250_521() {
    	this.setSeed(Ran250_521.ClockSeed());
    }
    
    public Ran250_521( long _seed ) {
    	this.setSeed(_seed);
    }
    
    public Ran250_521( Date d ) {
    	this.setSeed(d.getTime());
    }

	/*--------------------------------------------------------------------------*/
	void initGenerator() {
	  double rmod;
	  int i, ihlp, imask1, imask2;
	  int icyc, ncyc, nrest, ibas1, ibas2, ibas3;
	  long iseed;
	  int nrand;

	  // additions by anj:
	  iseed = getSeed();
	  nrand = getTotal();
	  random_numbers = new double[ getTotal() ];

	  if( iseed >= BIGINTEGER ) iseed = iseed%BIGINTEGER;
	  System.out.println(" D RAN250/521 Got seed "+iseed);
	  if(iseed <=0 )
	    {
	      System.err.print("Message from random number initialization:\n");
	      System.err.println("Please specify a seed greater than zero and smaller than "+BIGINTEGER);
	      System.exit(0);
	    }
	  if(nrand <=0)
	    {
	      System.err.print("Message from random number initialization:\n");
	      System.err.print("Please specify a positive number of random numbers\n");
	      System.exit(0);
	    }

	  rmod = (double) (iseed);

	/* Warm up the congruential generator */

	  for(i = 0; i < NWARM; ++i)
	    {
	      rmod = MULTIPLY * rmod;
	      rmod = rmod - ( (double) ( (int) (rmod * FACTOR) ) ) * BIGFLOAT;
	      ihlp = (int) (rmod + 0.1);      /* This is done to get rid of */
	      rmod = (double) (ihlp);         /* possible roundoff errors   */
	    }

	/* Allocate memory for the working arrays */

	  rand_w_array1 = new int[ BIGMAGIC1 + nrand ];
	  rand_w_array2 = new int[ BIGMAGIC2 + nrand ];

	/* Put congruential random numbers onto the working arrays */

	  for(i = 0; i < BIGMAGIC1; ++i)
	    {
	      rmod = MULTIPLY * rmod;
	      rmod = rmod - ( (double) ( (int) (rmod * FACTOR) ) ) * BIGFLOAT;
	      ihlp = (int) (rmod + 0.1);
	      rmod = (double) (ihlp);
	      rand_w_array1[i] = ihlp;
	    }

	  for(i = 0; i < BIGMAGIC2; ++i)
	    {
	      rmod = MULTIPLY * rmod;
	      rmod = rmod - ( (double) ( (int) (rmod * FACTOR) ) ) * BIGFLOAT;
	      ihlp = (int) (rmod + 0.1);
	      rmod = (double) (ihlp);
	      rand_w_array2[i] = ihlp;
	    }

	/* Linear independence of the bit columns for both generators. */
	/* Put ones on the main diagonal, and zeroes above.            */
	/* & is the bitwise AND                                        */
	/* | is the bitwise OR                                         */
	/* ^ is the bitwise XOR                                        */

	  imask1 = 1;
	  imask2 = BIGINTEGER;
	  for(i = NBIT - 2; i > 0; --i)
	    {
	      rand_w_array1[i] = ( rand_w_array1[i] | imask1 ) & imask2;
	      rand_w_array2[i] = ( rand_w_array2[i] | imask1 ) & imask2;
	      imask2 = imask2 ^ imask1;
	      imask1 = imask1 * 2;
	    }
	  rand_w_array1[0] = imask1;    /* This last element is treated separately */
	  rand_w_array2[0] = imask1;    /* in order to avoid overflow in imask1    */

	/* Warm up. Same structure as in vector_random_generator.      */
	/* Double loop structure to enable vectorization of inner loop */

	/* First, generator one */

	  ncyc  = nrand / SMALLMAGIC1;
	  nrest = nrand - SMALLMAGIC1 * ncyc;

	  ibas3 = BIGMAGIC1;                 /* position of first new random number */
	  ibas2 = BIGMAGIC1 - SMALLMAGIC1;   /* position of first input for this    */
	  ibas1 = 0;                         /* position of second input for this   */

	  for(icyc = 0; icyc < ncyc; ++icyc)
	    {
		//#pragma ivdep
	      for(i = 0; i < SMALLMAGIC1; ++i)
		{
		  rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] 
	                                   ^ rand_w_array1[ibas2 + i];
		}
	      ibas1 = ibas1 + SMALLMAGIC1;
	      ibas2 = ibas2 + SMALLMAGIC1;
	      ibas3 = ibas3 + SMALLMAGIC1;
	    }

	  if(nrest > 0)
	    {
		//#pragma ivdep
	      for(i = 0; i < nrest; ++i)
		{
		  rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] 
	                                   ^ rand_w_array1[ibas2 + i];
		}
	    }

	/* Put last elements to the beginning */

	  //#pragma ivdep
	  for(i = 0; i < BIGMAGIC1; ++i)
	    {
	      rand_w_array1[i] = rand_w_array1[nrand + i];
	    }

	/* Now the same for the second generator */

	  ncyc  = nrand / SMALLMAGIC2;
	  nrest = nrand - SMALLMAGIC2 * ncyc;

	  ibas3 = BIGMAGIC2;                 /* position of first new random number */
	  ibas2 = BIGMAGIC2 - SMALLMAGIC2;   /* position of first input for this    */
	  ibas1 = 0;                         /* position of second input for this   */

	  for(icyc = 0; icyc < ncyc; ++icyc)
	    {
		//#pragma ivdep
	      for(i = 0; i < SMALLMAGIC2; ++i)
		{
		  rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] 
	                                   ^ rand_w_array2[ibas2 + i];
		}
	      ibas1 = ibas1 + SMALLMAGIC2;
	      ibas2 = ibas2 + SMALLMAGIC2;
	      ibas3 = ibas3 + SMALLMAGIC2;
	    }

	  if(nrest > 0)
	    {
		//#pragma ivdep
	      for(i = 0; i < nrest; ++i)
		{
		  rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] 
	                                   ^ rand_w_array2[ibas2 + i];
		}
	    }

	/* Put last elements to the beginning */

	  //#pragma ivdep
	  for(i = 0; i < BIGMAGIC2; ++i)
	    {
	      rand_w_array2[i] = rand_w_array2[nrand + i];
	    }

	/* Initialization complete */

	  return;
	}

	/*--------------------------------------------------------------------------*/
	void generateNumbers() {
	  int i, icyc, ncyc, nrest, ibas1, ibas2, ibas3, nrand;

	  // Additions by anj:
	  nrand = getTotal();

	/* First, run generator one */

	  ncyc  = nrand / SMALLMAGIC1;
	  nrest = nrand - SMALLMAGIC1 * ncyc;

	  ibas3 = BIGMAGIC1;                 /* position of first new random number */
	  ibas2 = BIGMAGIC1 - SMALLMAGIC1;   /* position of first input for this    */
	  ibas1 = 0;                         /* position of second input for this   */

	  for(icyc = 0; icyc < ncyc; ++icyc)
	    {
		//#pragma ivdep
	      for(i = 0; i < SMALLMAGIC1; ++i)
		{
		  rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] 
	                                   ^ rand_w_array1[ibas2 + i];
		}
	      ibas1 = ibas1 + SMALLMAGIC1;
	      ibas2 = ibas2 + SMALLMAGIC1;
	      ibas3 = ibas3 + SMALLMAGIC1;
	    }

	  if(nrest > 0)
	    {
		//#pragma ivdep
	      for(i = 0; i < nrest; ++i)
		{
		  rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] 
	                                   ^ rand_w_array1[ibas2 + i];
		}
	    }

	/* Put last elements to the beginning */

	  //#pragma ivdep
	  for(i = 0; i < BIGMAGIC1; ++i)
	    {
	      rand_w_array1[i] = rand_w_array1[nrand + i];
	    }

	/* Now the same for the second generator */

	  ncyc  = nrand / SMALLMAGIC2;
	  nrest = nrand - SMALLMAGIC2 * ncyc;

	  ibas3 = BIGMAGIC2;                 /* position of first new random number */
	  ibas2 = BIGMAGIC2 - SMALLMAGIC2;   /* position of first input for this    */
	  ibas1 = 0;                         /* position of second input for this   */

	  for(icyc = 0; icyc < ncyc; ++icyc)
	    {
		//#pragma ivdep
	      for(i = 0; i < SMALLMAGIC2; ++i)
		{
		  rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] 
	                                   ^ rand_w_array2[ibas2 + i];
		}
	      ibas1 = ibas1 + SMALLMAGIC2;
	      ibas2 = ibas2 + SMALLMAGIC2;
	      ibas3 = ibas3 + SMALLMAGIC2;
	    }

	  if(nrest > 0)
	    {
		//#pragma ivdep
	      for(i = 0; i < nrest; ++i)
		{
		  rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] 
	                                   ^ rand_w_array2[ibas2 + i];
		}
	    }

	/* Put last elements to the beginning */

	  //#pragma ivdep
	  for(i = 0; i < BIGMAGIC2; ++i)
	    {
	      rand_w_array2[i] = rand_w_array2[nrand + i];
	    }

	/* Generate normalized random numbers:                        */
	/* Take output from generator one and combine it with         */
	/* that from generator two, via a simple XOR                  */

	  //#pragma ivdep
	  for(i = 0; i < nrand; ++i)
	    {
	      random_numbers[i] = FACTOR *(
		(rand_w_array1[i + BIGMAGIC1] ^ rand_w_array2[i + BIGMAGIC2]));
	    }

	  return;
	}

	  long getSeed() { return( seed ); } // grab the seed.

	  int getTotal() { return( total ); } // enumerate the number of numbers number.

	
	/* (non-Javadoc)
	 * @see edu.cornell.lassp.houle.RngPack.RandomElement#raw()
	 */
	public double raw() {
		if( entry == total ) {
			generateNumbers();
			entry = 0;
		}
		entry++;
		return( random_numbers[entry-1] );
	}
	
	  // Setting the seed:
	  private void setSeed( long nuseed ) {
	      seed = nuseed;
	      this.initGenerator();
	  }

	  // allow the user to set the total:
	  public void setTotal( int nutotal ) {
		  total = nutotal;
	      this.initGenerator();
	  }

}

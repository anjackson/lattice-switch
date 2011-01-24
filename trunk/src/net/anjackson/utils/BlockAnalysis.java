/**-------------------------------------------------------------
 * jBinLats - BlockAnalysis.java
 * net.anjackson.utils.BlockAnalysis
 * 
 * Created on 17-Feb-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.utils;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.Vector;
import edu.cornell.lassp.houle.RngPack.Ranmar;
;

/**
 * This class is given a long list of measured values of some observable.
 * The measurements may or may not be correlated.
 * 
 * <p>
 * When asked to analyse it, it breaks the sample down into blocks, and calculates
 * the mean/standard deviation etc for many different block sizes.
 * </p>
 * <p>
 * The block-size over which the error estimates are roughly constant 
 * will give a good estimate of the true value and error.
 * </p>
 * <p>
 * For more information see section 3.4 of ANJ's thesis.
 * </p>
 * <p>
 * To Do List
 * <ul>
 * <li>IDEA Could also attempt to calculate auto-correlation time.</li>
 * </ul>
 * </p>
 * 
 * @author ajackso1
 * @version $Id:$
 *
 */
public class BlockAnalysis implements Serializable {
	// The vector that will hold the list of measurements:
	/* Java 1.5:
	Vector<Double> times = new Vector<Double>();
	Vector<Double> values = new Vector<Double>();
	*/
	Vector times = new Vector();
	Vector values = new Vector();
	// Has the data been analysed?
	boolean analysed = false;
	// Properties of the data set:
	/* Java 1.5:
	Vector<Integer> bN = null;
	Vector<Double> bmean = null;
	Vector<Double> bstddev = null;
	Vector<Double> bstderr = null;
    */
	Vector bN = null;
	Vector bmean = null;
	Vector bstddev = null;
	Vector bstderr = null;
	double mean, stderr, stddev, varience;
	int bestblocksize;
	
	/**
	 * Add a new observation to the list:
	 * @param time The simulation time at which the measurement was taken.
	 * @param value The double-precision value of the measurement for this observable.
	 */
	public void add( double time, double value ) {
		analysed = false;
		times.add(new Double(time));
		values.add(new Double(value));
	}
	
	/**
	 * Determine whether this block analyser has been given any data.
	 * @return false if no data has been collected, true otherwise.
	 */
	public boolean hasData() {
		if (values.size() > 0) return true;
		return false;
	}

	/**
	 * Perform the data analysis.
	 */
	private void analyse() {
		if( !hasData() ) return;
		// The block-size parameters:
		int min_blocks = 5; // The minimum number of blocks to break the dataset into.
		int min_block_len = 1; // The smallest size block to use.
		String mstr = ""; // String used while manipulating m.
		// Vectors to store the results in:
		/* Java 1.5:
		bN = new Vector<Integer>();
		bmean = new Vector<Double>();
		bstderr = new Vector<Double>();
		bstddev = new Vector<Double>();
		*/
		bN = new Vector();
		bmean = new Vector();
		bstderr = new Vector();
		bstddev = new Vector();
		
		// Perform block analysis
		Double ix; // The ith value:
		double bx, bx2, x, x2, xx2; // Variables to hold the block x, x and x2 sums, and stderr.
		int ib = 0; // The number of blocks for the given block size.
		int m = min_block_len; // The block size - the number of data-points per block.
		do {
			// Debug output:
			// Loop over the dataset in blocks
			x = 0; x2 = 0; xx2 = 0.0;
			ib = 0;
			// While there is another block to analyse:
			while( m*(ib+1) <= values.size() ) {
				bx = 0.0; bx2 = 0.0;
				for( int i = m*ib; i < m*(ib+1); i++) {
					ix = (Double)values.get(i);
					bx += ix.doubleValue();
					bx2 += ix.doubleValue()*ix.doubleValue();
				}
				// Calculate <x>B and StdDev of the sums
				bx = bx/(double)m;
				bx2 = Math.sqrt( bx2/(double)m - bx*bx);
				// Calculate the mean:
				x += bx;
				// For the standard-error of the mean:
				x2 += bx*bx;
				// For the standard deviation of the distribution:
				xx2 += bx2;
				// Move on to next block:
				ib++;
			}
			// Calculate <x> and sigma(x) for this block
			if( ib > 0 ) {
				mean = x/(double)ib;
				x2 = Math.sqrt( (x2/(double)ib) - mean*mean);
				stderr = x2/Math.sqrt((double)ib);
				xx2 = xx2/(double)ib;
				if( new Double(stderr).isNaN() ) stderr = 0.0;
				// Store this result for this block size:
				bN.add(new Integer(m));
				bmean.add(new Double(mean));
				bstddev.add(new Double(xx2));
				bstderr.add(new Double(stderr));
				//System.out.println(" DEBUG "+m+" "+ib+" "+mean+" "+x2+" "+stderr+" "+xx2);
			}
			// Increase the block size, using 'add 1 to the first digit' at present.  
			// IDEA This should be easier in log_10(x) I think.
			mstr = Integer.toString(m);
			m  = (int)( (Integer.parseInt(mstr.substring(0,1))+1)*Math.pow(10.0, mstr.length()-1));
		} while ( ib != 0 && ib >= min_blocks );
		
		// Now go through and pick the 'best' result (should be were stderr is roughly stable)
		// TASK Currently using simple central-difference numerical integration, but this is risky and a suitable fitting function should be used to ensure stability.
		// IDEA Get Skewness and Kurtosis (3rd and 4th moments) as close to zero (Gaussian) as possible?
		// Skewness: http://en.wikipedia.org/wiki/Skewness
		// Kurtosis: http://en.wikipedia.org/wiki/Kurtosis
		int ibest = 0;
		if( bN.size() > 2 ) {
			// For each block size, compute the disagreement in stderr with the neighbouring results (differential wrt block size)
			double diff[] = new double[bN.size()];
			double mindiff = 0.0, fordiff, backdiff;
			// Differentiate via central difference.
			for( int iB = 1; iB < bN.size()-1; iB++ ) {
				// Central difference: (1/2){ [ (f(x+h) - f(x)) / h ] + [ (f(x) - f(x-h)) / h ] }
				fordiff = (((Double)bstderr.get(iB+1)).doubleValue() - ((Double)bstderr.get(iB)).doubleValue())/
				(((Integer)bN.get(iB+1)).intValue() - ((Integer)bN.get(iB)).doubleValue());
				backdiff = (((Double)bstderr.get(iB)).doubleValue() - ((Double)bstderr.get(iB-1)).doubleValue())/
				(((Integer)bN.get(iB)).intValue() - ((Integer)bN.get(iB-1)).intValue());
				diff[iB] = 0.5*(fordiff+backdiff);
				// Choose the one closest to zero.
				if( iB == 1 || mindiff > Math.abs(diff[iB])) {
					mindiff = Math.abs(diff[iB]);
					ibest = iB;
				}
			}
		}
		// Store the best estimate in the class members:
		bestblocksize = ((Integer)bN.get(ibest)).intValue();
		mean = ((Double)bmean.get(ibest)).doubleValue();
		stddev = ((Double)bstddev.get(ibest)).doubleValue();
		varience = stddev*stddev;
		stderr = ((Double)bstderr.get(ibest)).doubleValue();
		//System.out.println(" DEBUG Best block guessed to be "+bestblocksize+", "+mean+", "+stderr);
		
		// Set flag to indicated data has been analysed.
		analysed = true;
		return;
	}
	
	/**
	 * Write the plot of block-size v. measurement mean and error out to a file.
	 * @param filename The file to write the data out to.  X Y1 Y2 Y3 xmgrace-compatible format.
	 */
	public void writeBlockAnalysisGraph( String filename ) {
		if( !hasData() ) return;
		if( !analysed ) analyse();
		try {
			// Open the file:
			PrintStream out = new PrintStream( new FileOutputStream(filename));

			// Output a header:
			out.println("# blocksize <O> StdErr_b(O) StdDev_b(O)");
		
			// Loop over the block sizes and output the data:
			for( int i = 0; i < bN.size(); i++ ) {
				double m = ((Integer)bN.get(i)).intValue();
				double mmean = ((Double)bmean.get(i)).doubleValue();
				double mstddev = ((Double)bstddev.get(i)).doubleValue();
				double mstderr = ((Double)bstderr.get(i)).doubleValue();
				out.println(" "+m+" "+mmean+" "+mstderr+" "+mstddev);
			}
		
			// Close the file.
			out.close();
		} catch( Exception e ) {
			System.err.println("writeBlockAnalysisGraph("+filename+") failed due to exception "+e);
		}
		// Let the use know this was done:
		System.out.println("writeBlockAnalysisGraph("+filename+") succeeded.");
	}

	/**
	 * Get the best estimate for the mean value of the raw data
	 * @return The estimated mean (average)
	 */
	public double getMean() {
		if( !hasData() ) return 0;
		if( !analysed ) analyse();
		return mean;
	}

	/**
	 * Get the best estimate for the standard deviation of the raw data
	 * @return The standard deviation
	 */
	public double getStdDev() {
		if( !hasData() ) return 0;
		if( !analysed ) analyse();
		return stddev;
	}
	
	/**
	 * Get the best estimate for the standard error of the mean
	 * @return The standard error of the mean
	 */
	public double getStdErr() {
		if( !hasData() ) return 0;
		if( !analysed ) analyse();
		return stderr;
	}
	
	/**
	 * Get the best estimate for the varience
	 * @return The varience of the measurements
	 */
	public double getVarience() {
		if( !hasData() ) return 0;
		if( !analysed ) analyse();
		return varience;
	}
	
	/**
	 * Write the raw measurements to a data file.
	 * @param filename The file to write the data out to.  Time and value.
	 */
	public void writeRawDatafile( String filename ) {
		if( !hasData() ) return;
		try {
			// Open the file:
			PrintStream out = new PrintStream( new FileOutputStream(filename));

			// Output a header:
			out.println("# Simulation_Time Measurement_Value ");
		
			// Loop over the block sizes and output the data:
			for( int i = 0; i < values.size(); i++ ) {
				double t = ((Double)times.get(i)).doubleValue();
				double x = ((Double)values.get(i)).doubleValue();
				out.println(" "+t+" "+x);
			}
		
			// Close the file.
			out.close();
		} catch( Exception e ) {
			System.err.println("writeRawDatafile("+filename+") failed due to exception "+e);
		}
		
	}

	/**
	 * Empty vectors of any recorded data.
	 */
	public void clear() { 
		values.clear();
		times.clear();
		analysed = false;
	}
	
	/**
	 * A simple algorithm to run some basic tests on the block-analysis code.
	 * IDEA: Use proper unit-testing framework (JUnit) for this?
	 */
	private void runTests() {
		// Count the number of failed tests:
		int fails = 0;
		// Rounding tolerance for tests:
		double fptol = 1e-10;
		// A random number generator for the tests:
		Ranmar rng = new Ranmar();
		
		// Simple test, is the average of lots of the same number correct?
		for( int t = 0; t < 10; t++ ) {
			// Clear the vectors:
			this.clear();
			// Generate a number and a number of times to try:
			double val = Math.pow(rng.raw(),rng.uniform(-20,20));
			int len = (int)rng.uniform(1000,10000);
			// Add the number to the vectors:
			for( int i = 0; i < len; i++ ) {
				this.add(i,val);
			}
			// Get the mean value:
			double tval = this.getMean();
			// Check it makes sense:
			if( Math.abs((val-tval)/val) > fptol ) {
				System.out.println("TEST FAILED: Mean of "+len+" x "+val+" gave "+tval+" Delta = "+(val-tval));
				fails++;
				this.writeBlockAnalysisGraph("last-failed-test.dat");
			} else {
				System.out.println("TEST PASSED: Mean of "+len+" x "+val+" Delta = "+(val-tval));			
			}
		}
		
		// A slightly harder test:  Generate a lot of points for a known distribution and check the mean and stdev all make sense:
		this.clear();
		// The StdDev and mean of the Gaussian:
		double stddev = 7.0, mean = 3.0;
		// The number of sample points to draw:
		int len = 400000;
		// Statistical error tolerance: The factor of ten gives sufficent room to pass the test, but hopefully not too much.
		double statol = 10/Math.sqrt(len);
		for( int i = 0; i < len; i++ ) {
			this.add(i,rng.gaussian(stddev)+mean);
		}
		if( len < 2000 ) this.writeRawDatafile("stat-test-raw.dat");
		double tmean = this.getMean();
		double tstddev = this.getStdDev();
		if( Math.abs((tmean-mean)) > statol ) {
			System.out.println("TEST FAILED: Mean of ("+mean+", "+stddev+") gave "+tmean+" Delta = "+(mean-tmean)+" > "+statol);
			fails++;
			this.writeBlockAnalysisGraph("last-failed-test.dat");
		} else {
			System.out.println("TEST PASSED: Mean of ("+mean+", "+stddev+") gave "+tmean);
		}
		if( Math.abs((tstddev-stddev)) > statol ) {
			System.out.println("TEST FAILED: Stddev of ("+mean+", "+stddev+") gave "+tstddev+" Delta = "+(stddev-tstddev)+" > "+statol);
			fails++;
			this.writeBlockAnalysisGraph("last-failed-test.dat");
		} else {
			System.out.println("TEST PASSED: Stddev of ("+mean+", "+stddev+") gave "+tstddev);
		}
		
		
		
		// Check the number of fails and exit if there are many:
		if( fails > 0 ) {
			System.out.println("UNIT TEST FAILURE! "+fails+" tests failed.");
			System.exit(1);
		}		
	}

	/**
	 *  Routine to load a t v. x data set into the vector arrays.
	 */
	public void loadFromFile(String filename) {
		// Clear any stored content:
		this.clear();
		// Open the file:
		
		// Load in each value as t and x:
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			String line, parts[];
			while ( in.ready() ) {
				// Read and parse a line:	
				line = in.readLine().trim();
				int index = 0;
				// Skip comments:
				if( !"#".equals(line.substring(0,1))) {
					parts = line.split("\\s+");
					// Two-columns, assume T v X.
					if( parts.length >= 2 ) {
						this.add( Double.parseDouble(parts[0]), Double.parseDouble(parts[1]));
					} else {
						// Assuming take mean of first column:
						this.add(index,Double.parseDouble(parts[0]));
					}
					// Increment the index:
					index++;
				}
			}
			/* close the file */
			in.close();
		} catch( Exception e ) {
			throw( new RuntimeException("BlockAnalysis.loadFromFile Failed! - ",e));
		}
	}
	
	/**
	 * A simple runner, to perform test and file block analyses
	 * @param args The command line arguments.
	 */
	public static void main( String[] args ) {
		// Check CLargs
		if( args.length != 1 ) {
			System.err.println("Usage: BlockAnalysis [<filename>|-test]");
			System.exit(1);
		}

		// Create an instance of this class:
		BlockAnalysis ban = new BlockAnalysis();
		
		// Either run some tests:
		if( "-test".equalsIgnoreCase(args[0])) {
			ban.runTests();
		} else {
			// Or analyse a file.
			ban.loadFromFile(args[0]);
			System.out.println(" Average is estimated to be "+ban.getMean()+" +/- "+ban.getStdErr());
			System.out.println(" i.e. 1StdDev interval is ["+(ban.getMean()-ban.getStdErr())+", "+(ban.getMean()+ban.getStdErr())+"]");
			System.out.println(" Standard deviation is estimated to be "+ban.getStdDev());
			String graphFilename = args[0]+".ba.dat";
			System.out.println(" Writing block analysis graph to "+graphFilename);
			ban.writeBlockAnalysisGraph(graphFilename);
		}
	}
}

package net.anjackson.physics.lattices;

import java.io.Serializable;

import net.anjackson.maths.*;

/**
 * This is the abstract class defining a Lattice.  It contains the basic fields
 * to describe a simple repeating spacial lattice in three-dimensions.  It includes 
 * a number of methods for accessing these properties.
 * 
 * It is based on a generation-on-demand form.  i.e. you start up an instance (or
 * rather a sub-class which has implemented the generate function), and the lattice
 * is generated when you request the points themselves, using the properties as 
 * they are at that time.
 * 
 * As well as implementing the generate function, sub-classes should also declare 
 * the unit-cell dimensions at unit site seperation.  For multiple lattices, this 
 * is the major-lattice seperation.
 * 
 * IDEA Move UC declaration into a method called by constructor, so that it is always set.
 *
 * @author ajackso1
 * @version $Id: Lattice.java 969 2006-08-15 17:50:02Z anj $
 *
 */
public abstract class Lattice implements Serializable {
	// Unit cell dimensions:
	protected double uc[] = new double[3];  
	// The radius/seperation of neighbours:
	// By convention, sub-classes should set this to the unit cell for site-seperation = 1.0.
	protected double uc_sep = 1.0;
	// Number of sites per unit cell:
	protected int nspc = 1;
	// The name of this lattice
	protected String name = null;

	// Store the lattice as an array of vectors:
	protected VectorD3[] lat = null;
	// Unit cell scaling:
	protected double ucx = 1.0, ucy = 1.0, ucz = 1.0;
	// Lattice origin offset:
	protected double loox = 0.0, looy = 0.0, looz = 0.0;
	// System size in terms of the number of cells:
	protected int nx = -1, ny = -1, nz = -1;
	// Number of sites:
	protected int na = -1;

	/**
	 * Has the class been given enough information to be used?
	 * @return false if there is not enough information to call generate().
	 */
	protected boolean isInitialised() {
		if( nx < 0 || ny < 0 || nz < 0 ) return false;
		return true;
	}
	
	/**
	 * Look up the name of the lattice - a simple text string to describe it.
	 * @return The name of this lattice, e.g. 'FCC'
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * Set a new unit cell scaling factor.
	 * 
	 * @param ucx The scaling of the unit cell in the x direction.
	 * @param ucy The scaling of the unit cell in the y direction.
	 * @param ucz The scaling of the unit cell in the z direction.
	 */
	public void setUnitCellScaling( double ucx, double ucy, double ucz ) {
		this.ucx = ucx;
		this.ucy = ucy;
		this.ucz = ucz;
	}

	/**
	 * Get the size of the unit cell:
	 * @return A double[3] containing the size in the x,y,z directions ([0,1,2]).
	 */
	public double[] getBasicUnitCell() {
		return uc;
	}
	
	/**
	 * Get the size of the unit cell - including the scaling factors.
	 * @return A double[3] containing the size in the x,y,z directions ([0,1,2]).
	 */
	public double[] getUnitCell() {
		double ucs[] = new double[3];
		ucs[0] = uc[0]*ucx;
		ucs[1] = uc[1]*ucy;
		ucs[2] = uc[2]*ucz;
		return ucs;
	}
	
	/**
	 * Set a new origin for the lattice.
	 * @param _loox The origin (x-coordinate).
	 * @param _looy The origin (y-coordinate).
	 * @param _looz The origin (z-coordinate).
	 */
	public void setOriginOffset( double _loox, double _looy, double _looz ) {
		this.loox = _loox;
		this.looy = _looy;
		this.looz = _looz;
	}
	
	/**
	 * Tell the generator how big a crystal to generate.
	 * @param _nx Number of unit-cells in the x-direction.
	 * @param _ny Number of unit-cells in the y-direction.
	 * @param _nz Number of unit-cells in the z-direction.
	 */
	public void setSizeInCells( int _nx, int _ny, int _nz ) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		na = nspc * nx*ny*nz;
		if( ! checkSize(nx,ny,nz) ) {
			System.out.println("This lattice ("+name+") cannot be constructed with the given size ("+nx+", "+ny+", "+nz+").");			
			System.exit(1);
		}
	}

	/**
	 * This is the general generator for a 3D lattice of nx by ny by nx
	 * unit cells.
	 */
	abstract protected void generate();

	/**
	 * Get the vector array that describes this lattice.
	 * @return The array of VectorD3 site vectors.
	 */
	public VectorD3[] getLattice() {
    	if( !isInitialised() ) return null;
    	if( lat == null ) generate();
		return lat;
	}

	/**
	 * Return the number of lattices - one for this base class:
	 * @return 1 (one lattice)
	 */
    public int getNumberOfLattices() {
    	return 1;
    }
 
    /**
     * Return the total number of sites:
     * @return Total number of sites in this lattice.
     */
    public int getNumberOfSites() {
    	if( !isInitialised() ) return -1;
    	return na;
    }
    
    public double[] getSystemSize() {
    	if( !isInitialised() ) return null;
    	double[] ldim = new double[3];
    	double ucell[] = getBasicUnitCell();
		ldim[0] = nx*ucell[0];
		ldim[1] = ny*ucell[1];
		ldim[2] = nz*ucell[2];
		return ldim;
    }
    
    /**
     * Can this lattice be constructed (in PBC) at this size:
     * @param nx Cells in X direction.
     * @param ny Cells in Y direction.
     * @param nz Cells in Z direction.
     * @return FALSE if this size would not be allowed.  True if all is well.
     */
    public boolean checkSize(int nx, int ny, int nz) {
    	//System.out.println("Using generic checkSize...");
    	return true;
    }
    
}

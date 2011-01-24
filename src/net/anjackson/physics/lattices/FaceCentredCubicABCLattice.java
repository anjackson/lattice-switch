package net.anjackson.physics.lattices;

/**
 * 
 * This class generates face-centred-cubic lattices. It is a special case of
 * the generic close-packing generator {@link ClosePackedLattice}.  i.e. it is
 * built from close-packed layers in an ABC pattern.
 *
 * @author ajackso1
 * @version $Id: FaceCentredCubicABCLattice.java 981 2006-08-21 14:40:43Z anj $
 *
 */
public class FaceCentredCubicABCLattice extends ClosePackedLattice {
	
	// The constructor, used to define the unit-seperation unit cell.
	public FaceCentredCubicABCLattice() {
		uc[0] = Math.sqrt(3.0)/2.0;
		uc[1] = 1.0;
		uc[2] = Math.sqrt(2.0/3.0);
		uc_sep = 1.0;
		name = "fcc";
	}
	
	/**
	 * Extend the generate mechanism, specify the
	 * FCC lattice to the CP generator.
	 */
	protected void generate() {
		cpGenerator( nx, ny, nz, ClosePackedLattice.FCC);
	}

	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public boolean checkSize(int nx, int ny, int nz) {
    	if( nx%2 != 0 ) return false;
    	if( nz%3 != 0 ) return false;
    	return true;
    }
	
}

package net.anjackson.physics.lattices;

/**
 * 
 * This class generates hexagonal close-packed lattices. It is a special case of
 * the generic close-packing generator {@link ClosePackedLattice}.
 *
 * @author ajackso1
 * @version $Id: HexagonalClosePackedLattice.java 968 2006-08-15 16:54:11Z anj $
 *
 */
public class HexagonalClosePackedLattice extends ClosePackedLattice {
	
	// The constructor, used to define the unit-seperation unit cell.
	public HexagonalClosePackedLattice() {
		uc[0] = Math.sqrt(3.0)/2.0;
		uc[1] = 1.0;
		uc[2] = Math.sqrt(2.0/3.0);
		uc_sep = 1.0;
		name = "hcp";
	}
	
	/**
	 * Extend the generate mechanism, specify the
	 * HCP lattice to the CP generator.
	 */
	protected void generate() {
		cpGenerator( nx, ny, nz, ClosePackedLattice.HCP);
	}
	
	/**
	 * Disallow system sizes that will not fit in the PBC.
	 * @return FALSE if the system size is bad.
	 */
    public boolean checkSize(int nx, int ny, int nz) {
    	if( nx%2 != 0 ) return false;
    	if( nz%2 != 0 ) return false;
    	return true;
    }

}

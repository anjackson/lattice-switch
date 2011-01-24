/**-------------------------------------------------------------
 * jBinLats - VectorD3.java
 * net.anjackson.maths.VectorD3
 * 
 * Created on 23-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.maths;

import java.io.Serializable;

/**
 * A simple double-precision vector.
 *
 * @author ajackso1
 * @version $Id: VectorD3.java 573 2006-02-20 20:03:56Z anj $
 *
 */
public class VectorD3 implements Serializable {
   public double x,y,z;
   
    public VectorD3( ) {
    	x=0; y=0; z=0;
    }
    
    public VectorD3( double _x, double _y, double _z ) {
    	x = _x;
    	y = _y;
    	z = _z;
    }
    
    /**
     * Multiple this vector by a scalar, returning a new object:
     * @param scalar The amount to scale all coordinates by:
     */
    public VectorD3 times( double scalar ) {
    	return new VectorD3(this.x*scalar, this.y*scalar, this.z*scalar);
    }
    
    /**
     * Create a clone of this object
     */
    public VectorD3 copy() {
    	return new VectorD3(this.x, this.y, this.z);    	
    }
}

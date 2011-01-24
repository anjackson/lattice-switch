/**-------------------------------------------------------------
 * jBinLats - VectorI3.java
 * net.anjackson.maths.VectorI3
 * 
 * Created on 23-Jan-2006 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2006 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.maths;

import java.io.Serializable;

/**
 * A simple integer vector.
 *
 * @author ajackso1
 * @version $Id:$
 *
 */
public class VectorI3 implements Serializable {
	public int x,y,z;
	
	public VectorI3( ) {
		x=0; y=0; z=0;
	}
	
	public VectorI3( int _x, int _y, int _z ) {
		x = _x;
		y = _y;
		z = _z;
	}
}

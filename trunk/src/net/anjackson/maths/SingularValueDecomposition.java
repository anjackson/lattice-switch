/**-------------------------------------------------------------
 * jBinLats - SingularValueDecomposition.java
 * net.anjackson.maths.SingularValueDecomposition
 * 
 * Created on 14-Dec-2005 by ajackso1.
 * -------------------------------------------------------------
 * Copyright (c) 2005 Andrew N. Jackson.
 * Licence: 
 * -------------------------------------------------------------
 */
package net.anjackson.maths;

/**
 * This is a straightforward port of the Numerical Recipies Singular Value
 * Decomposition routines.
 * <p>
 * See Numerical Recipies in C - <a href="http://www.nr.com">www.nr.com</a>.
 * Also see <a href="http://www.library.cornell.edu/nr/bookcpdf/c2-6.pdf">NR section 2.6 Singular Value Decomposition</a>.
 * </p>
 * 
 * @author Ported by ajackso1
 * @version $Id: SingularValueDecomposition.java 518 2006-01-25 18:21:41Z anj $
 *
 */
public class SingularValueDecomposition {

	/**
	 * Straight port of the Numerical Recipies Singular Value Decomposition routine.
	 * ---
	 * Given a matrix a[1..m][1..n], this routine computes its singular value 
	 * decomposition, A = U·W·V T. The matrix U replaces a on output. The 
	 * diagonal matrix of singular values W is output as a vector w[1..n]. 
	 * The matrix V (not the transpose V T ) is output as v[1..n][1..n].
	 * ---
	 * @param a The matrix to be de-composed.  Replaced with U on return.
	 * @param m The first dimension of the matrix (number of slowest-running indexes).
	 * @param n The second dimension of the matrix (number of fastest-running indexes).
	 * @param w On return, this contains the singular values.
	 * @param v On return, this contains the V matrix of the decomposition.
	 */
	public static void svdcmp(double a[][], int m, int n, double w[], double v[][]) {
		int flag,i,its,j,jj,k,l= 0,nm = 0;
		double anorm,c,f,g,h,s,scale,x,y,z,rv1[];
	
		rv1 = new double[n];
		g=scale=anorm=0.0;
		// Householder reduction to bidiagonal form.
		for (i=1;i<=n;i++) {
			l=i+1;
			rv1[i]=scale*g;
			g=s=scale=0.0;
			if (i <= m) {
				for (k=i;k<=m;k++) scale += Math.abs(a[k][i]);
				if (scale!=0) {
					for (k=i;k<=m;k++) {
						a[k][i] /= scale;
						s += a[k][i]*a[k][i];
					}
					f=a[i][i];
					g = -Utils.SIGN(Math.sqrt(s),f);
					h=f*g-s;
					a[i][i]=f-g;
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
					for (k=i;k<=m;k++) a[k][i] *= scale;
				}
			}
			w[i]=scale *g;
			g=s=scale=0.0;
			if (i <= m && i != n) {
				for (k=l;k<=n;k++) scale += Math.abs(a[i][k]);
				if (scale!=0) {
					for (k=l;k<=n;k++) {
						a[i][k] /= scale;
						s += a[i][k]*a[i][k];
					}
					f=a[i][l];
					g = -Utils.SIGN(Math.sqrt(s),f);
					h=f*g-s;
					a[i][l]=f-g;
					for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
					for (k=l;k<=n;k++) a[i][k] *= scale;
				}
			}
			anorm=Utils.DMAX(anorm,(Math.abs(w[i])+Math.abs(rv1[i])));
		}
		// Accumulation of right-hand transformations.
		for (i=n;i>=1;i--) {
			if (i < n) {
				if (g!=0) {
					for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g; // Double division to avoid possible underflow.
					for (j=l;j<=n;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
						for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
					}
				}
				for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
			}
			v[i][i]=1.0;
			g=rv1[i];
			l=i;
		}
		// Accumulation of left-hand transformations.
		for (i=Utils.IMIN(m,n);i>=1;i--) {
			l=i+1;
			g=w[i];
			for (j=l;j<=n;j++) a[i][j]=0.0;
			if (g!=0) {
				g=1.0/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (j=i;j<=m;j++) a[j][i] *= g;
			} else for (j=i;j<=m;j++) a[j][i]=0.0;
			++a[i][i];
		}
		// Diagonalization of the bidiagonal form: Loop over
		// singular values, and over allowed iterations.
		for (k=n;k>=1;k--) {
			for (its=1;its<=30;its++) {
				flag=1;
				for (l=k;l>=1;l--) {
					// Test for splitting.
					// Note that rv1[1] is always zero.
					nm=l-1;
					if ((double)(Math.abs(rv1[l])+anorm) == anorm) {
						flag=0;
						break;
					}
					if ((double)(Math.abs(w[nm])+anorm) == anorm) break;
				}
				if (flag!=0) {
					// Cancellation of rv1[l], if l > 1.
					c=0.0;
					s=1.0;
					for (i=l;i<=k;i++) {
						f=s*rv1[i];
						rv1[i]=c*rv1[i];
						if ((double)(Math.abs(f)+anorm) == anorm) break;
						g=w[i];
						h=Utils.dpythag(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s = -f*h;
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
				z=w[k];
				// Convergence.
				if (l == k) {
					if (z < 0.0) {
						// Singular value is made non-negative.
						w[k] = -z;
						for (j=1;j<=n;j++) v[j][k] = -v[j][k];
					}
					break;
				}
				if (its == 30) System.out.println("no convergence in 30 dsvdcmp iterations");
				x=w[l]; // Shift from bottom 2-by-2 minor.
				nm=k-1;
				y=w[nm];
				g=rv1[nm];
				h=rv1[k];
				f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
				g=Utils.dpythag(f,1.0);
				f=((x-z)*(x+z)+h*((y/(f+Utils.SIGN(g,f)))-h))/x;
				c=s=1.0; // Next QR transformation:
				for (j=l;j<=nm;j++) {
					i=j+1;
					g=rv1[i];
					y=w[i];
					h=s*g;
					g=c*g;
					z=Utils.dpythag(f,h);
					rv1[j]=z;
					c=f/z;
					s=h/z;
					f=x*c+g*s;
					g = g*c-x*s;
					h=y*s;
					y *= c;
					for (jj=1;jj<=n;jj++) {
						x=v[jj][j];
						z=v[jj][i];
						v[jj][j]=x*c+z*s;
						v[jj][i]=z*c-x*s;
					}
					z=Utils.dpythag(f,h);
					w[j]=z;
					// Rotation can be arbitrary if z = 0.
					if (z!=0) {
						z=1.0/z;
						c=f*z;
						s=h*z;
					}
					f=c*g+s*y;
					x=c*y-s*g;
					for (jj=1;jj<=m;jj++) {
						y=a[jj][j];
						z=a[jj][i];
						a[jj][j]=y*c+z*s;
						a[jj][i]=z*c-y*s;
					}
				}
				rv1[l]=0.0;
				rv1[k]=f;
				w[k]=x;
			}
		}
	}

	/**
	 * This is a straight port of the SVD back-substitution algorithm from
	 * Numerical Recipies.
	 * ---
	 * Solves A·X = B for a vector X, where A is specified by the 
	 * arrays u[1..m][1..n], w[1..n],v[1..n][1..n] as returned by svdcmp. 
	 * m and n are the dimensions of a, and will be equal for square matrices. 
	 * b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
	 * No input quantities are destroyed, so the routine may be called 
	 * sequentially with different b’s.
	 * ---
	 * @param u The decomposed matrix, U, from svdcmp.
	 * @param w The singular values, W, from svdcmp.
	 * @param v The decomposed matrix, V, from svdcmp.
	 * @param m The slowest-running dimension of the A/U/V matricies.
	 * @param n The fastest-running dimension of the A/U/V matricies.
	 * @param b The right-hand side of the A.x = b we wish to solve.
	 * @param x Om return, the values of x for which A.x = b.
	 */
	public static void svbksb(double u[][], double w[], double v[][], int m, int n, double b[], double x[]) {
		int jj,j,i;
		double s,tmp[];
	
		tmp = new double[n];
		// Calculate U^T.B.
		for (j=1;j<=n;j++) {
			s=0.0;
			// Nonzero result only if wj is nonzero.
			if (w[j]!=0) {
				for (i=1;i<=m;i++) s += u[i][j]*b[i];
				s /= w[j]; // This is the divide by wj.
			}
			tmp[j]=s;
		}
		// Matrix multiply by V to get answer.
		for (j=1;j<=n;j++) {
			s=0.0;
			for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
			x[j]=s;
		}
		
	}

}

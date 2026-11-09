/*
 * Code written by Richard Brent and John Burkardt and released under
 * the LGPL. Original available from:
 *
 * https://people.math.sc.edu/Burkardt/c_src/praxis/praxis.html
 *
 * Modified, 2026, by Michael Carley to add passing of data to
 * objective functions, workspace to replace malloc, and to stop
 * multiple #include s
 * 
 */

#ifndef __PRAXIS_H_INCLUDED__
#define __PRAXIS_H_INCLUDED__

void minfit ( int n, double tol, double a[], double q[] );
double praxis ( double t0, double h0, int n, int prin, double x[], 
		double f ( double x[], int n , void *d), void *d,
		double *work );

#endif /*__PRAXIS_H_INCLUDED__*/

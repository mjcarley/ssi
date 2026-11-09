/*
 * Code translated to C and modified from the original C++ of Ramsey,
 * Potter, and Hansen, Ray Bilinear Patch Intersections, Journal of
 * Graphics Tools, 2004, https://doi.org/10.1080/10867651.2004.10504896
 *
 * Original code from
 * 
 * http://shaunramsey.com/research/bp/
 *
 * Original code released under MIT License, below
 */

/*
 * The ray bilinear patch intersection software are "Open Source"
 * according to the MIT License located at:
 * http://www.opensource.org/licenses/mit-license.php

 * Copyright (c) 2003 Shaun David Ramsey, Kristin Potter, Charles
 * Hansen

 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sel copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

/*
 * modification to allow for special case of z-component of tangent
 * equal to zero, Michael Carley, 2026
 * 
 */

/* #include <stdio.h> */
#include <math.h>

#include "bilinear.h"

#define ABS(_x) ((_x) > 0 ? (_x) : (-(_x)))

int bilinear_patch_srf_eval(BILINEAR_REAL *P00,
			    BILINEAR_REAL *P10,			   
			    BILINEAR_REAL *P11,
			    BILINEAR_REAL *P01,
			    BILINEAR_REAL u, BILINEAR_REAL v,
			    BILINEAR_REAL *x)

{
  x[0] =
    (1.0 - u)*(1.0 - v)*P00[0] +
    (1.0 - u)*(      v)*P01[0] +
    (      u)*(1.0 - v)*P10[0] +
    (      u)*(      v)*P11[0] ;
  x[1] =
    (1.0 - u)*(1.0 - v)*P00[1] +
    (1.0 - u)*(      v)*P01[1] +
    (      u)*(1.0 - v)*P10[1] +
    (      u)*(      v)*P11[1] ;
  x[2] =
    (1.0 - u)*(1.0 - v)*P00[2] +
    (1.0 - u)*(      v)*P01[2] +
    (      u)*(1.0 - v)*P10[2] +
    (      u)*(      v)*P11[2] ;

  return 0 ;
}

/* static void bilinear_tangent_u(BILINEAR_REAL *P00, */
/* 			       BILINEAR_REAL *P10,			    */
/* 			       BILINEAR_REAL *P11, */
/* 			       BILINEAR_REAL *P01, */
/* 			       BILINEAR_REAL u, BILINEAR_REAL v, */
/* 			       BILINEAR_REAL *tanu) */

/* { */
/*   tanu[0] = (1.0 -v)*(P10[0] - P00[0]) + v*(P11[0] - P01[0]) ; */
/*   tanu[1] = (1.0 -v)*(P10[1] - P00[1]) + v*(P11[1] - P01[1]) ; */
/*   tanu[2] = (1.0 -v)*(P10[2] - P00[2]) + v*(P11[2] - P01[2]) ; */
  
/*   return ; */
/* } */

/* static void bilinear_tangent_v(BILINEAR_REAL *P00, */
/* 			       BILINEAR_REAL *P10,			    */
/* 			       BILINEAR_REAL *P11, */
/* 			       BILINEAR_REAL *P01, */
/* 			       BILINEAR_REAL u, BILINEAR_REAL v, */
/* 			       BILINEAR_REAL *tanv) */

/* { */
/*   tanv[0] = (1.0 -u)*(P01[0] - P00[0]) + u*(P11[0] - P10[0]) ; */
/*   tanv[1] = (1.0 -u)*(P01[1] - P00[1]) + u*(P11[1] - P10[1]) ; */
/*   tanv[2] = (1.0 -u)*(P01[2] - P00[2]) + u*(P11[2] - P10[2]) ; */
  
/*   return ; */
/* } */

/* static void bilinear_normal(BILINEAR_REAL *P00, */
/* 			    BILINEAR_REAL *P10,			    */
/* 			    BILINEAR_REAL *P11, */
/* 			    BILINEAR_REAL *P01, */
/* 			    BILINEAR_REAL u, BILINEAR_REAL v, */
/* 			    BILINEAR_REAL *n) */

/* { */
/*   BILINEAR_REAL tanu[3], tanv[3] ; */

/*   bilinear_tangent_u(P00, P01, P11, P10, u, v, tanu) ; */
/*   bilinear_tangent_u(P00, P01, P11, P10, u, v, tanv) ; */

/*   bilinear_vector_cross(n,tanu,tanv) ; */
  
/*   return ; */
/* } */

static BILINEAR_REAL getu(BILINEAR_REAL v, BILINEAR_REAL M1, BILINEAR_REAL M2,
			  BILINEAR_REAL J1,BILINEAR_REAL J2,
			  BILINEAR_REAL K1, BILINEAR_REAL K2,
			  BILINEAR_REAL R1, BILINEAR_REAL R2)
{

  BILINEAR_REAL denom = (v*(M1-M2)+J1-J2);
  BILINEAR_REAL d2 = (v*M1+J1);

  if( ABS(denom) > ABS(d2)) // which denominator is bigger
    {
      return (v*(K2-K1)+R2-R1)/denom;
    }
  return -(v*K1+R1)/d2;
}

static BILINEAR_REAL computet(BILINEAR_REAL *dir, BILINEAR_REAL *orig,
			      BILINEAR_REAL *srfpos)

{
  // if x is bigger than y and z
  if ( ABS(dir[0]) >= ABS(dir[1]) && ABS(dir[0]) >= ABS(dir[2]) )
    return (srfpos[0] - orig[0]) / dir[0];
  // if y is bigger than x and z
  else if ( ABS(dir[1]) >= ABS(dir[2])) // && ABS(dir[1]) >= ABS(dir[0]))
    return (srfpos[1] - orig[1]) / dir[1];
  // otherwise x isn't bigger than both and y isn't bigger than both
  else  //if( ABS(dir[2]) >= ABS(dir[0]) && ABS(dir[2]) >= ABS(dir[1]))
    return (srfpos[2] - orig[2]) / dir[2];    
}


// a x ^2 + b x + c = 0
// in this case, the root must be between min and max
// it returns the # of solutions found
// x = [ -b +/- sqrt(b*b - 4 *a*c) ] / 2a
// or x = 2c / [-b +/- sqrt(b*b-4*a*c)]
int QuadraticRoot(BILINEAR_REAL a, BILINEAR_REAL b, BILINEAR_REAL c, 
		  BILINEAR_REAL min, BILINEAR_REAL max, BILINEAR_REAL *u)
{
  u[0] = u[1] = min-min; // make it lower than min
  if (a == 0.0) // then its close to 0
    {
      if (b != 0.0) // not close to 0
	{
	  u[0] = - c / b;
	  if(u[0] > min && u[0] < max) //its in the interval
	    return 1; //1 soln found
	  else  //its not in the interval
	    return 0;
	}
      else
	return 0;
    }
  BILINEAR_REAL d = b*b - 4*a*c; //discriminant
  if(d <= 0.0) // single or no root
    {
      if(d == 0.0) // close to 0
	{
	  u[0] = -b / a;
	  if(u[0] > min && u[0] < max) // its in the interval
	    return 1;
	  else //its not in the interval
	    return 0;
	}
      
      else // no root d must be below 0
	return 0;
    }
  
  BILINEAR_REAL q = -0.5  * (b + copysign(sqrt(d),b));
  u[0] = c / q;
  u[1] = q / a;

  if(     (u[0] > min && u[0] < max)
	  &&  (u[1] > min && u[1] < max))
    return 2;
  else if(u[0] > min && u[0] < max) //then one wasn't in interval
    return 1;
  else if(u[1] > min && u[1] < max)
    {  // make it easier, make u[0] be the valid one always
      BILINEAR_REAL dummy;
      dummy = u[0];
      u[0] = u[1];
      u[1] = dummy; // just in case somebody wants to check it
      return 1;
    }
  return 0;

}

//+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
//             RayPatchIntersection
// intersect rays of the form p = r + t q where t is the parameter
// to solve for. With the patch pointed to by *this
// for valid intersections:
// place the u,v intersection point in uv[0] and uv[1] respectively.
// place the t value in uv[2]
// return true to this function
// for invalid intersections - simply return false uv values can be 
// anything
//+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+

int bilinear_ray_patch_intersection(BILINEAR_REAL *P00,
				    BILINEAR_REAL *P10,
				    BILINEAR_REAL *P11,
				    BILINEAR_REAL *P01,
				    BILINEAR_REAL *r,
				    BILINEAR_REAL *q,
				    BILINEAR_REAL *uv)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Equation of the patch:
  // P(u, v) = (1-u)(1-v)P00 + (1-u)vP01 + u(1-v)P10 + uvP11
  // Equation of the ray:
  // R(t) = r + tq
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BILINEAR_REAL pos1[3], pos2[3] ; //BILINEAR_REAL *pos = ro + t*rd;
  int num_sol; // number of solutions to the quadratic
  BILINEAR_REAL vsol[2]; // the two roots from quadraticroot
  BILINEAR_REAL t2,u; // the t values of the two roots

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Variables for substitition
  // a = P11 - P10 - P01 + P00
  // b = P10 - P00
  // c = P01 - P00
  // d = P00  (d is shown below in the #ifdef raypatch area)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

  if ( ABS(q[2]) < 1e-9 ) {
    /*allow for possible spurious miss when z component of tangent
      is zero: rotate the coordinates and call the function again*/
    BILINEAR_REAL pt00[3], pt10[3], pt11[3], pt01[3], rt[3], qt[3] ;

    pt00[0] = P00[1] ; pt00[1] = P00[2] ; pt00[2] = P00[0] ; 
    pt10[0] = P10[1] ; pt10[1] = P10[2] ; pt10[2] = P10[0] ; 
    pt11[0] = P11[1] ; pt11[1] = P11[2] ; pt11[2] = P11[0] ; 
    pt01[0] = P01[1] ; pt01[1] = P01[2] ; pt01[2] = P01[0] ; 
    rt[0] = r[1] ; rt[1] = r[2] ; rt[2] = r[0] ; 
    qt[0] = q[1] ; qt[1] = q[2] ; qt[2] = q[0] ; 
    
    return bilinear_ray_patch_intersection(pt00, pt10, pt11, pt01, rt, qt, uv) ;
  }
  
  // Find a w.r.t. x, y, z
  BILINEAR_REAL ax = P11[0] - P10[0] - P01[0] + P00[0];
  BILINEAR_REAL ay = P11[1] - P10[1] - P01[1] + P00[1];
  BILINEAR_REAL az = P11[2] - P10[2] - P01[2] + P00[2];


  // Find b w.r.t. x, y, z
  BILINEAR_REAL bx = P10[0] - P00[0];
  BILINEAR_REAL by = P10[1] - P00[1];
  BILINEAR_REAL bz = P10[2] - P00[2];

  // Find c w.r.t. x, y, z
  BILINEAR_REAL cx = P01[0] - P00[0];
  BILINEAR_REAL cy = P01[1] - P00[1];
  BILINEAR_REAL cz = P01[2] - P00[2];

  BILINEAR_REAL rx = r[0];
  BILINEAR_REAL ry = r[1];
  BILINEAR_REAL rz = r[2];

  // Retrieve the xyz of the q part of ray
  BILINEAR_REAL qx = q[0];
  BILINEAR_REAL qy = q[1];
  BILINEAR_REAL qz = q[2];

  // Find d w.r.t. x, y, z - subtracting r just after  
  BILINEAR_REAL dx = P00[0] - r[0];
  BILINEAR_REAL dy = P00[1] - r[1];
  BILINEAR_REAL dz = P00[2] - r[2];
  

  // Find A1 and A2
  BILINEAR_REAL A1 = ax*qz - az*qx;
  BILINEAR_REAL A2 = ay*qz - az*qy;

  // Find B1 and B2
  BILINEAR_REAL B1 = bx*qz - bz*qx;
  BILINEAR_REAL B2 = by*qz - bz*qy;

  // Find C1 and C2
  BILINEAR_REAL C1 = cx*qz - cz*qx;
  BILINEAR_REAL C2 = cy*qz - cz*qy;

  // Find D1 and D2
  BILINEAR_REAL D1 = dx*qz - dz*qx;
  BILINEAR_REAL D2 = dy*qz - dz*qy;
 
  /* BILINEAR_REAL *dir(qx,qy,qz), orig(rx,ry,rz); */
  BILINEAR_REAL dir[3], orig[3] ;
  BILINEAR_REAL A = A2*C1 - A1*C2;
  BILINEAR_REAL B = A2*D1 -A1*D2 + B2*C1 -B1*C2;
  BILINEAR_REAL C = B2*D1 - B1*D2;

  dir[0] = qx ; dir[1] = qy ; dir[2] = qz ;
  orig[0] = rx ; orig[1] = ry ; orig[2] = rz ;
  uv[0] = -2 ; uv[1] = -2 ; uv[2] = -2;

  num_sol = QuadraticRoot(A,B,C,-ray_epsilon,1+ray_epsilon,vsol);

  switch(num_sol) {
  case 0:
    return 0; // no solutions found
  case 1:
    uv[1] = vsol[0] ;
    uv[0] = getu(uv[1],A2,A1,B2,B1,C2,C1,D2,D1) ;
    /* pos1 = SrfEval(uv[0],uv[1]); */
    bilinear_patch_srf_eval(P00, P10, P11, P01, uv[0], uv[1], pos1) ;
    uv[2] = computet(dir,orig,pos1) ;
    if(uv[0] < 1+ray_epsilon && uv[0] > -ray_epsilon && uv[2] > 0)//vars okay?
      return 1;
    else
      return 0; // no other soln - so ret false
  case 2: // two solutions found
    uv[1] = vsol[0] ;
    uv[0] = getu(uv[1],A2,A1,B2,B1,C2,C1,D2,D1) ;
    /* pos1 = SrfEval(uv[0],uv[1]); */
    bilinear_patch_srf_eval(P00, P10, P11, P01, uv[0], uv[1], pos1) ;
    uv[2] = computet(dir,orig,pos1) ; 
    if(uv[0] < 1+ray_epsilon && uv[0] > -ray_epsilon && uv[2] > 0)
      {
	u = getu(vsol[1],A2,A1,B2,B1,C2,C1,D2,D1);
	if(u < 1+ray_epsilon && u > ray_epsilon)
	  {
	    /* pos2 = SrfEval(u,vsol[1]); */
	    bilinear_patch_srf_eval(P00, P10, P11, P01, u, vsol[1], pos2) ;
	    t2 = computet(dir,orig,pos2);
	    if(t2 < 0 || uv[2] < t2) // t2 is bad or t1 is better
	      return 1; 
	    // other wise both t2 > 0 and t2 < t1
	    uv[1] = vsol[1];
	    uv[0] = u ;
	    uv[2] = t2 ;
	    return 1;
	  }
	return 1; // u2 is bad but u1 vars are still okay
	}
    else // doesn't fit in the root - try other one
      {
	uv[1] = vsol[1] ;
	uv[0] = getu(vsol[1],A2,A1,B2,B1,C2,C1,D2,D1) ;
	/* pos1 = SrfEval(uv[0],uv[1]); */
	bilinear_patch_srf_eval(P00, P10, P11, P01, uv[0], uv[1], pos1) ;
	uv[2] = computet(dir,orig,pos1) ; 
	if(uv[0] < 1+ray_epsilon && uv[0] > -ray_epsilon &&uv[2] > 0)
	  return 1;
	  else
	    return 0;
      }
    break;
  }
  
  /* cout << " ERROR: We don't get here in Ray Patch Intersection" << endl; */
  return 0; 
}


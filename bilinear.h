/*
 * Code translated to C and modified from the original C++ of Ramsey,
 * Potter, and Hansen, Ray Bilinear Patch Intersections, Journal of
 * Graphics Tools, 2004, https://doi.org/10.1080/10867651.2004.10504896
 *
 * Original code from
 * 
 * http://shaunramsey.com/research/bp/
 *
 * Original code released under MIT License
 */


#ifndef __BILINEAR_H_INCLUDED__
#define __BILINEAR_H_INCLUDED__

#define ray_epsilon 1e-12

#define bilinear_vector_scalar(_A,_B)				\
  (((_A)[0])*((_B)[0])+((_A)[1])*((_B)[1])+((_A)[2])*((_B)[2]))

#define bilinear_vector_length(_A)		\
  (sqrt(bilinear_vector_scalar((_A),(_A))))

#define bilinear_vector_normalize(_A)				\
  do {								\
  BILINEAR_REAL _len = bilinear_vector_length((_A)) ;		\
  ((_A)[0]) /= _len ; ((_A)[1]) /= _len ; ((_A)[2]) /= _len ;	\
  } while (0)

#define bilinear_vector_cross(_C,_A,_B)				\
  ((_C)[0] = (_A)[1]*(_B)[2] - (_A)[2]*(_B)[1],			\
   (_C)[1] = (_A)[2]*(_B)[0] - (_A)[0]*(_B)[2],			\
   (_C)[2] = (_A)[0]*(_B)[1] - (_A)[1]*(_B)[0])

int bilinear_patch_srf_eval(BILINEAR_REAL *P00,
			    BILINEAR_REAL *P10,			   
			    BILINEAR_REAL *P11,
			    BILINEAR_REAL *P01,
			    BILINEAR_REAL u, BILINEAR_REAL v,
			    BILINEAR_REAL *x) ;
int bilinear_ray_patch_intersection(BILINEAR_REAL *P00,
				    BILINEAR_REAL *P10,
				    BILINEAR_REAL *P11,
				    BILINEAR_REAL *P01,
				    BILINEAR_REAL *r,
				    BILINEAR_REAL *q,
				    BILINEAR_REAL *uv) ;

#define bilinear_shape_func(_u,_v,_L00,_L10,_L11,_L01)		\
  do {								\
    (_L00) = (1.0 - (_u))*(1.0 - (_v)) ;			\
    (_L01) = (1.0 - (_u))*(      (_v)) ;			\
    (_L10) = (      (_u))*(1.0 - (_v)) ;			\
    (_L11) = (      (_u))*(      (_v)) ;			\
  } while (0)

#endif /*__BILINEAR_H_INCLUDED__*/

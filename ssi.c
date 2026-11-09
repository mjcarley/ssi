#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <float.h>

#include "bilinear.h"
#include "praxis.h"

#include "crawl.h"
#include "ssi.h"

void qsort_r(void *base, size_t nmemb, size_t size,
	     int (*compar)(const void *, const void *, void *),
	     void *arg);

static void swap_double_array(double *a, double *b, int n)

{
  int i ;
  double tmp ;

  for ( i = 0 ; i < n ; i ++ ) {
    tmp = a[i] ; a[i] = b[i] ; b[i] = tmp ;
  }
  
  return ;
}

double ssi_mapping_linear(double u, void *data)

{
  double *ulim = data ;
  double um ;

  fprintf(stderr, "untested code\n") ;
  exit(1) ;
  
  um = ulim[0] + u*(ulim[1] - ulim[0]) ;
  
  return um ;
}

double ssi_mapping_cosine(double u, void *data)

{
  double th ;
  double um ;

  th = M_PI + u*M_PI ;
  um = 0.5*(1+cos(th)) ;
  
  return um ;
}

double ssi_mapping_split_cosine(double u, void *data)

{
  double th ;
  double um ;

  if ( u < 0.5  ) {
    th = -M_PI + 2.0*M_PI*u ;
    um = 0.5*(cos(th) - 1) ;
  } else {
    th = 2.0*(u - 0.5)*M_PI ;
    um = 0.5*(1.0 - cos(th)) ;
  }

  return um ;
}

ssi_bvh_t *ssi_bvh_new(int depth)

{
  ssi_bvh_t *b ;
  int h, nb, nbt, np ;
  
  b = (ssi_bvh_t *)malloc(sizeof(ssi_bvh_t)) ;

  for ( h = nbt = 0 ; h <= depth ; h ++ ) {
    /*number of boxes at this level*/
    nb = (1 << h) ; nb *= nb ;
    b->off[h] = nbt ;
    nbt += nb ;
  }

  /*number of points in sampled surface*/
  np = (1 << depth) + 1 ; np *= np ;
  
  b->x = (double *)malloc((np + 2*nbt)*3*sizeof(double)) ;
  b->xb = &(b->x[3*np]) ;
  
  ssi_bvh_umin(b) = 0 ; ssi_bvh_umax(b) = 1 ;
  ssi_bvh_vmin(b) = 0 ; ssi_bvh_vmax(b) = 1 ;

  b->depth = depth ;

  /* ssi_bvh_mapping_u(b) = SSI_MAPPING_LINEAR ; */
  /* ssi_bvh_mapping_v(b) = SSI_MAPPING_LINEAR ; */
  ssi_bvh_mapping_u(b) = ssi_mapping_linear ;
  ssi_bvh_mapping_data_u(b) = &(ssi_bvh_umin(b)) ;
  ssi_bvh_mapping_v(b) = ssi_mapping_linear ;
  ssi_bvh_mapping_data_v(b) = &(ssi_bvh_vmin(b)) ;
  
  return b ;
}

int ssi_bvh_fill_grid(ssi_bvh_t *bvh, ssi_surface_func_t sfunc, void *data)

{
  int i, j, n, ret ;
  double u, v ;
  ssi_mapping_func_t umfunc, vmfunc ;
  void *umdata, *vmdata ;
  
  umfunc = ssi_bvh_mapping_u(bvh) ;
  umdata = ssi_bvh_mapping_data_u(bvh) ;
  vmfunc = ssi_bvh_mapping_v(bvh) ;
  vmdata = ssi_bvh_mapping_data_v(bvh) ;
  
  n = 1 << bvh->depth ;
  for ( i = 0 ; i <= n ; i ++ ) {
    /* u = ssi_bvh_u_global(bvh, (double)i/n) ; */
    u = umfunc((double)i/n, umdata) ;
    for ( j = 0 ; j <= n ; j ++ ) {
      /* v = ssi_bvh_v_global(bvh, (double)j/n) ; */
      v = vmfunc((double)j/n, vmdata) ;
      if ( (ret = sfunc(u, v, &(bvh->x[(i*(n+1)+j)*3]), NULL, NULL, data)) ) {
	return ret ;
      }
    }
  }

  bvh->sfunc = sfunc ; bvh->sdata = data ;
  
  return 0 ;
}

int ssi_bvh_write_grid(FILE *f, ssi_bvh_t *bvh)

{
  int i, j, n ;
  double u, v, *x ;

  n = 1 << bvh->depth ;
  for ( i = 0 ; i <= n ; i ++ ) {
    u = ssi_bvh_u_global(bvh, (double)i/n) ;
    for ( j = 0 ; j <= n ; j ++ ) {
      v = ssi_bvh_v_global(bvh, (double)j/n) ;
      x = &(bvh->x[(i*(n+1)+j)*3]) ;
      fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e\n",
	      x[0], x[1], x[2], u, v) ;
    }
  }

  return 0 ;
}

int ssi_surface_write_grid(FILE *f, 
			   ssi_surface_func_t sfunc2, void *data,
			   double umin, double umax,
			   double vmin, double vmax,
			   int nu, int nv)

{
  int i, j ;
  double u, v, x[3] ;

  for ( i = 0 ; i <= nu ; i ++ ) {
    u = umin + (umax - umin)*i/nu ;
    for ( j = 0 ; j <= nv ; j ++ ) {
      v = vmin + (vmax - vmin)*j/nv ;
      sfunc2(u, v, x, NULL, NULL, data) ;
      fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e\n",
	      x[0], x[1], x[2], u, v) ;
    }
  }

  return 0 ;
}

int ssi_bvh_patch_bbox(double *x1, double *x2, double *x3, double *x4,
		       double del, double *b1, double *b2)

{
  double xmin, xmax, ymin, ymax, zmin, zmax ;

  xmin = MIN(x1[0],MIN(x2[0],MIN(x3[0],x4[0]))) ;
  xmax = MAX(x1[0],MAX(x2[0],MAX(x3[0],x4[0]))) ;
  ymin = MIN(x1[1],MIN(x2[1],MIN(x3[1],x4[1]))) ;
  ymax = MAX(x1[1],MAX(x2[1],MAX(x3[1],x4[1]))) ;
  zmin = MIN(x1[2],MIN(x2[2],MIN(x3[2],x4[2]))) ;
  zmax = MAX(x1[2],MAX(x2[2],MAX(x3[2],x4[2]))) ;

  b1[0] = xmin - del ; b1[1] = ymin - del ; b1[2] = zmin - del ;
  b2[0] = xmax + del ; b2[1] = ymax + del ; b2[2] = zmax + del ;
  
  return 0 ;
}

int ssi_bvh_make_bboxes(ssi_bvh_t *bvh, int level, double del)

{
  int i, j, i1, j1, idx, idx1, idx2, idx3, idx4, n, off, off1, n1 ;
  double xmin, xmax, ymin, ymax, zmin, zmax ;

  off = bvh->off[level] ;
  n = 1 << level ;
  if ( level == bvh->depth ) {
    /*at full tree depth, we have to box the surface patches*/
    for ( i = 0 ; i < n ; i ++ ) {
      for ( j = 0 ; j < n ; j ++ ) {
	idx  = off + i*n + j ;
	idx1 = (i+0)*(n+1) + j + 0 ;
	idx2 = (i+1)*(n+1) + j + 0 ;
	idx3 = (i+1)*(n+1) + j + 1 ;
	idx4 = (i+0)*(n+1) + j + 1 ;
	ssi_bvh_patch_bbox(&(bvh->x[3*idx1]),
			   &(bvh->x[3*idx2]),
			   &(bvh->x[3*idx3]),
			   &(bvh->x[3*idx4]),
			   del,
			   &(bvh->xb[6*idx+0]),
			   &(bvh->xb[6*idx+3])) ;
      }
    }

    return 0 ;
  }

  /*at higher levels we box the boxes at the next level down*/
  off1 = bvh->off[level+1] ;
  n1   = 1 << (level + 1) ;
  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ ) {
      idx  = off + i*n + j ;
      idx1 = off1 + (2*i+0)*n1 + 2*j + 0 ;
      xmin = bvh->xb[6*idx1+0] ;
      ymin = bvh->xb[6*idx1+1] ;
      zmin = bvh->xb[6*idx1+2] ;
      xmax = bvh->xb[6*idx1+3] ;
      ymax = bvh->xb[6*idx1+4] ;
      zmax = bvh->xb[6*idx1+5] ;

      for ( i1 = 0 ; i1 < 2 ; i1 ++ ) {
	for ( j1 = 0 ; j1 < 2 ; j1 ++ ) {
	  idx1 = off1 + (2*i+i1)*n1 + 2*j + j1 ;
	  xmin = MIN(xmin,bvh->xb[6*idx1+0]) ;
	  ymin = MIN(ymin,bvh->xb[6*idx1+1]) ;
	  zmin = MIN(zmin,bvh->xb[6*idx1+2]) ;
	  xmax = MAX(xmax,bvh->xb[6*idx1+3]) ;
	  ymax = MAX(ymax,bvh->xb[6*idx1+4]) ;
	  zmax = MAX(zmax,bvh->xb[6*idx1+5]) ;
	}
      }

      bvh->xb[6*idx+0] = xmin ;
      bvh->xb[6*idx+1] = ymin ;
      bvh->xb[6*idx+2] = zmin ;
      bvh->xb[6*idx+3] = xmax ;
      bvh->xb[6*idx+4] = ymax ;
      bvh->xb[6*idx+5] = zmax ;
    }
  }
  
  return 0 ;
}

int ssi_bvh_write_bboxes(FILE *f, ssi_bvh_t *bvh, int level)

{
  int i, j, idx, n, off ;

  off = bvh->off[level] ;
  n = 1 << level ;
  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ ) {
      idx  = off + i*n + j ;
      fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
	      bvh->xb[6*idx+0], bvh->xb[6*idx+1], bvh->xb[6*idx+2], 
	      bvh->xb[6*idx+3], bvh->xb[6*idx+4], bvh->xb[6*idx+5]) ;
    }
  }
  
  return 0 ;
}

int ssi_bbox_ray_slab_test(double *p0, double *p1,
			   double *orig, double *invr,
			   double rtmin, double rtmax)

{
  /*
   * from Peter Shirley, Ingo Wald & Adam Marrs, Ray Axis-Aligned
   * Bounding Box Intersection, 2021
   * 
   * https://doi.org/10.1007/978-1-4842-7185-8_2
   * 
   */
  double tlower[3], tupper[3], tmins[4], tmaxes[4], tboxmin, tboxmax ;

  tlower[0] = (p0[0] - orig[0])*invr[0] ;
  tlower[1] = (p0[1] - orig[1])*invr[1] ;
  tlower[2] = (p0[2] - orig[2])*invr[2] ;
  tupper[0] = (p1[0] - orig[0])*invr[0] ;
  tupper[1] = (p1[1] - orig[1])*invr[1] ;
  tupper[2] = (p1[2] - orig[2])*invr[2] ;

  tmins[0] = MIN(tlower[0], tupper[0]) ; tmaxes[0] = MAX(tlower[0], tupper[0]) ;
  tmins[1] = MIN(tlower[1], tupper[1]) ; tmaxes[1] = MAX(tlower[1], tupper[1]) ;
  tmins[2] = MIN(tlower[2], tupper[2]) ; tmaxes[2] = MAX(tlower[2], tupper[2]) ;
  tmins[3] =                     rtmin ; tmaxes[3] =                     rtmax ;

  tboxmin = MAX(tmins[0],MAX(tmins[1],MAX(tmins[2],tmins[3]))) ;
  tboxmax = MIN(tmaxes[0],MIN(tmaxes[1],MIN(tmaxes[2],tmaxes[3]))) ;
  
  return tboxmin < tboxmax ;
}

static int intersection_check(ssi_bvh_t *bvh, int level, int nb,
			      int i, int j,
			      double *orig, double *r, double *invr,
			      double rtmin, double rtmax, 
			      int *hits, double *uvhits, int uvstr,
			      int *nhits)

{
  int idx = i*nb + j ;
  int off = bvh->off[level] ;
  double *p00, *p10, *p11, *p01, *uv ;
  
  if ( ssi_bbox_ray_slab_test(&(bvh->xb[6*(off+idx)+0]),
			      &(bvh->xb[6*(off+idx)+3]),
			      orig, invr, rtmin, rtmax) == 0 ) {
    return 0 ;
  }

  /*slab hit, check for panel hit or child nodes*/
  if ( level == bvh->depth ) {
    /*check panel on leaf node*/
    p00 = &(bvh->x[3*((i+0)*(nb+1) + j + 0)]) ;
    p10 = &(bvh->x[3*((i+1)*(nb+1) + j + 0)]) ;
    p11 = &(bvh->x[3*((i+1)*(nb+1) + j + 1)]) ;
    p01 = &(bvh->x[3*((i+0)*(nb+1) + j + 1)]) ;
    uv = &(uvhits[uvstr*(*nhits)]) ;
    if ( bilinear_ray_patch_intersection(p00, p10, p11, p01, orig, r, uv) ) {
      if ( uv[2] >= rtmin && uv[2] <= rtmax ) {
	hits[(*nhits)] = idx ;
	(*nhits) ++ ;
	double x[3] ;
	bilinear_patch_srf_eval(p00, p10, p11, p01, uv[0], uv[1], x) ;
      }
    }

    return 0 ;
  }

  /*not a leaf node: check the children*/
  intersection_check(bvh, level+1, 2*nb, 2*i+0, 2*j+0, orig, r, invr,
		     rtmin, rtmax, hits, uvhits, uvstr, nhits) ;
  intersection_check(bvh, level+1, 2*nb, 2*i+1, 2*j+0, orig, r, invr,
		     rtmin, rtmax, hits, uvhits, uvstr, nhits) ;
  intersection_check(bvh, level+1, 2*nb, 2*i+1, 2*j+1, orig, r, invr,
		     rtmin, rtmax, hits, uvhits, uvstr, nhits) ;
  intersection_check(bvh, level+1, 2*nb, 2*i+0, 2*j+1, orig, r, invr,
		     rtmin, rtmax, hits, uvhits, uvstr, nhits) ;
  
  return 0 ;
}

int ssi_bvh_ray_intersection(ssi_bvh_t *bvh,
			     double *orig, double *r,
			     double rtmin, double rtmax,
			     double *uv, int uvstr)

{
  double invr[3], L00, L10, L11, L01, d, u, v ;
  int hits[128], i, nb, i0, j0, ni ;
  
  invr[0] = 1.0/r[0] ; invr[1] = 1.0/r[1] ; invr[2] = 1.0/r[2] ; 

  ni = 0 ;
  intersection_check(bvh, 0, 1, 0, 0, orig, r, invr, rtmin, rtmax,
		     hits, uv, uvstr, &ni) ;

  nb = 1 << bvh->depth ; d = 1.0/nb ;
  /*convert local coordinates to global*/
  for ( i = 0 ; i < ni ; i ++ ) {
    i0 = hits[i]/nb ; j0 = hits[i] % nb ;

    u = uv[i*uvstr+0] ; v = uv[i*uvstr+1] ;
    bilinear_shape_func(u, v, L00, L10, L11, L01) ;
    
    uv[i*uvstr+0] = d*((i0+0)*L00 + (i0+1)*L10 + (i0+1)*L11 + (i0+0)*L01) ;
    uv[i*uvstr+1] = d*((j0+0)*L00 + (j0+0)*L10 + (j0+1)*L11 + (j0+1)*L01) ;
  }
  
  return ni ;
}

int ssi_bvh_ray_intersection_slow(ssi_bvh_t *bvh,
				  double *orig, double *r,
				  double rtmin, double rtmax,
				  double *uv, int uvstr, int *ni)

{
  int nb, i, j ;
  double *p00, *p10, *p11, *p01 ;

  nb = 1 << bvh->depth ;

  for ( i = 0 ; i < nb ; i ++ ) {
    for ( j = 0 ; j < nb ; j ++ ) {
      p00 = &(bvh->x[3*((i+0)*(nb+1) + j + 0)]) ;
      p10 = &(bvh->x[3*((i+1)*(nb+1) + j + 0)]) ;
      p11 = &(bvh->x[3*((i+1)*(nb+1) + j + 1)]) ;
      p01 = &(bvh->x[3*((i+0)*(nb+1) + j + 1)]) ;
      bilinear_ray_patch_intersection(p00, p10, p11, p01, orig, r, uv) ;
    }
  }
  
  return 0 ;
}

static double refine_func(double *x, int n, void *d)

{
  void **data = d ;
  ssi_bvh_t *bvh1 = data[0] ;
  ssi_bvh_t *bvh2 = data[1] ;
  double uvfix = *((double *)data[2]) ;
  int fix = *((int *)data[3]) ;
  int i ;
  double uv[4], u1, v1, u2, v2, x1[3], x2[3], R ;

  uv[fix] = uvfix ;
  for ( i = 0 ; i < fix ; i ++ ) {
    uv[i] = x[i] ;
  }
  for ( i = fix+1 ; i < 4 ; i ++ ) {
    uv[i] = x[i-1] ;
  }
  
  u1 = ssi_bvh_u_global(bvh1,uv[0]) ; v1 = ssi_bvh_v_global(bvh1,uv[1]) ;
  u2 = ssi_bvh_u_global(bvh2,uv[2]) ; v2 = ssi_bvh_v_global(bvh2,uv[3]) ;

  bvh1->sfunc(u1, v1, x1, NULL, NULL, bvh1->sdata) ;
  bvh2->sfunc(u2, v2, x2, NULL, NULL, bvh2->sdata) ;  

  R = ssi_vector_distance2(x1, x2) ;
  
  return R ;
}

static int refine_intersection(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
			       double *uvst, double tol, int n,
			       int fix, double *work)

{
  double h0, x[3], r ;
  int i, nc ;
  void *data[8] ;

  if ( tol == 0 ) return 0 ;

  nc = ( fix < 4 ? 3 : 4 ) ;
  
  data[0] = bvh1 ; data[1] = bvh2 ; data[3] = &fix ;
  data[2] = &(uvst[fix]) ;
  for ( i =       0 ; i < fix ; i ++ ) x[i  ] = uvst[i] ;
  for ( i = fix + 1 ; i <   4 ; i ++ ) x[i-1] = uvst[i] ;

  h0 = 1.0/n ;
  
  r = praxis(tol, h0, nc, 0, x, refine_func, data, work) ;

  for ( i =       0 ; i < fix ; i ++ ) uvst[i] = x[i  ] ;
  for ( i = fix + 1 ; i <   4 ; i ++ ) uvst[i] = x[i-1] ;

  if ( r > tol ) return 1 ;

  if ( uvst[0] < 0 || uvst[0] > 1 ) return 1 ;
  if ( uvst[1] < 0 || uvst[1] > 1 ) return 1 ;
  if ( uvst[2] < 0 || uvst[2] > 1 ) return 1 ;
  if ( uvst[3] < 0 || uvst[3] > 1 ) return 1 ;
  
  return 0 ;
}

static int ray_intersection(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
			    int n,
			    int i0, int j0, int i1, int j1,
			    int fix, double *uvst)

/*
 * find intersection of a grid line on bvh2 with a surface patch on bvh1
 */
  
{
  double *x0, *x1, r[3], uv[32] ;
  
  x0 = &(bvh2->x[(i0*(n+1)+j0)*3]) ;
  x1 = &(bvh2->x[(i1*(n+1)+j1)*3]) ;
  ssi_vector_diff(r,x1,x0) ;

  if ( ssi_bvh_ray_intersection(bvh1, x0, r, 0, 1, uv, 3) ) {
    uvst[0] = uv[0] ;
    uvst[1] = uv[1] ;

    uvst[2] = (double)i0/n ;
    uvst[3] = (double)j0/n ;
    /* if ( fix == 3 ) uvst[2] += uv[2]*(i1 - i0) ; */
    /* if ( fix == 2 ) uvst[3] += uv[2]*(j1 - j0) ; */
    /* uvst[3] = d*j0 + uv[2]*(j1 - j0) ; */
    /* if ( uvst[3] <= 1.0 && uvst[3] >= 0.0 ) return 1 ; */
    return 1 ;
  }

  return 0 ;
}

static int find_intersection(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
			     int n,
			     int i0, int j0, int i1, int j1,
			     double *uvst, double tol, int fix,
			     double *work)

{
  if ( ray_intersection(bvh1, bvh2, n, i0, j0, i1, j1, fix, uvst) == 0 ) {
    return 0 ;
  }

  if ( refine_intersection(bvh1, bvh2, uvst, tol, n, fix, work) == 0 )  {
    return 1 ;
  }

  return 0 ;
}

static int find_intersection_swap(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
				  int n,
				  int i0, int j0, int i1, int j1,
				  double *uvst, double tol, int fix,
				  double *work)

{
  if ( ray_intersection(bvh1, bvh2, n, i0, j0, i1, j1, fix, uvst) == 0 ) {
    return 0 ;
  }

  if ( refine_intersection(bvh1, bvh2, uvst, tol, n, fix, work) == 0 )  {
    swap_double_array(&(uvst[0]), &(uvst[2]), 2) ;
    return 1 ;
  }

  return 0 ;
}

static int compare_points(const void *a, const void *b, void *data)

{
  const double *c1 = a, *c2 = b ;

  if ( c1[0] < c2[0] ) return -1 ;
  if ( c1[0] > c2[0] ) return  1 ;
  if ( c1[1] < c2[1] ) return -1 ;
  if ( c1[1] > c2[1] ) return  1 ;

  if ( c1[2] < c2[2] ) return -1 ;
  if ( c1[2] > c2[2] ) return  1 ;
  if ( c1[3] < c2[3] ) return -1 ;
  if ( c1[3] > c2[3] ) return  1 ;

  return 0 ;
}

int ssi_surface_surface_intersection(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
				     double *c, int cstr, int *nc, int ncmax,
				     double tol, double *work)

/*
 * minimum size of work is 27 doubles
 */
  
{
  int i, j, n ;

  if ( cstr < 4 ) {
    fprintf(stderr, "%s: cstr (%d) must be at least 4\n", __FUNCTION__, cstr) ;
    exit(1) ;
  }
  
  n = (1 << bvh2->depth) ;

  for ( i = 1 ; i < n ; i ++ ) {
    for ( j = 1 ; j < n ; j ++ ) {
      (*nc) += find_intersection(bvh1, bvh2, n, i+0, j+0, i+0, j+1,
				 &(c[(*nc)*cstr]), tol, 2, work) ;
      (*nc) += find_intersection(bvh1, bvh2, n, i+0, j+0, i+1, j+0,
				 &(c[(*nc)*cstr]), tol, 3, work) ;
    }
  }

  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j <= n ; j += n ) {
      (*nc) += find_intersection(bvh1, bvh2, n, j+0, i+0, j+0, i+1,
				 &(c[(*nc)*cstr]), tol, 2, work) ;
      (*nc) += find_intersection(bvh1, bvh2, n, i+0, j+0, i+1, j+0,
				 &(c[(*nc)*cstr]), tol, 3, work) ;
    }
  }

  n = (1 << bvh1->depth) ;

  for ( i = 1 ; i < n ; i ++ ) {
    for ( j = 1 ; j < n ; j ++ ) {
      (*nc) += find_intersection_swap(bvh2, bvh1, n, i+0, j+0, i+0, j+1,
				      &(c[(*nc)*cstr]), tol, 2, work) ;
      (*nc) += find_intersection_swap(bvh2, bvh1, n, i+0, j+0, i+1, j+0,
				      &(c[(*nc)*cstr]), tol, 3, work) ;
    }
  }

  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j <= n ; j += n ) {
      (*nc) += find_intersection_swap(bvh2, bvh1, n, j+0, i+0, j+0, i+1,
				      &(c[(*nc)*cstr]), tol, 2, work) ;
      (*nc) += find_intersection_swap(bvh2, bvh1, n, i+0, j+0, i+1, j+0,
				      &(c[(*nc)*cstr]), tol, 3, work) ;
      /* (*nc) += find_intersection_swap(bvh2, bvh1, n, j+0, i+0, j+0, i+1, */
      /* 				      &(c[(*nc)*cstr]), tol, 3, work) ; */
      /* (*nc) += find_intersection_swap(bvh2, bvh1, n, i+0, j+0, i+1, j+0, */
      /* 				      &(c[(*nc)*cstr]), tol, 2, work) ; */
    }
  }
  
  return 0 ;
}

int ssi_curve_sort_points(double *c, int cstr, int nc, int *dup, int *nd)

{
  int i, j ;

  *nd = 0 ;
  qsort_r(c, nc, cstr*sizeof(double), compare_points, NULL) ;

  /* return 0 ; */
  
  for ( i = 0 ; i < nc ; i ++ ) {
    for ( j = i + 1 ; j < nc - (*nd) ; j ++ ) {
      if ( ((c[i*cstr+0] == c[j*cstr+0]) &&
	    (c[i*cstr+1] == c[j*cstr+1]))
	   /* || */
	   /* ((c[i*cstr+2] == c[j*cstr+2]) && */
	   /*  (c[i*cstr+3] == c[j*cstr+3])) */
	   ) {
	/* swap_double_array(&(c[j*cstr]), &(c[(j+1)*cstr]), cstr) ; */
	/* swap_double_array(&(c[(j+1)*cstr]), &(c[(nc-(*nd)-1)*cstr]), cstr) ; */
	swap_double_array(&(c[j*cstr]), &(c[(nc-(*nd)-1)*cstr]), cstr) ;
	dup[2*(*nd)+0] = i ; dup[2*(*nd)+1] = nc-(*nd)-1 ;
	(*nd) ++ ;
      }
    }  
  }

  return 0 ;
}

static int point_is_duplicate(int *dup, int nd, int p, int *d)

{
  int i ;
  
  for ( i = 0 ; i < nd ; i ++ ) {
    if ( (dup[2*i+0] == p) || (dup[2*i+1] == p) ) {
      *d = i ; return 1 ;
    }
  }
  
  return 0 ;
}

static void closest_point(double *c, int cstr, int p, int dup[],
			  int *pc, int *pf)

{
  double di, dj ;
  int i, j ;

  i = dup[0] ; j = dup[1] ;
  di =
    (c[p*cstr+0] - c[i*cstr+0])*(c[p*cstr+0] - c[i*cstr+0]) + 
    (c[p*cstr+1] - c[i*cstr+1])*(c[p*cstr+1] - c[i*cstr+1]) + 
    (c[p*cstr+2] - c[i*cstr+2])*(c[p*cstr+2] - c[i*cstr+2]) + 
    (c[p*cstr+3] - c[i*cstr+3])*(c[p*cstr+3] - c[i*cstr+3]) ;
  dj =
    (c[p*cstr+0] - c[j*cstr+0])*(c[p*cstr+0] - c[j*cstr+0]) + 
    (c[p*cstr+1] - c[j*cstr+1])*(c[p*cstr+1] - c[j*cstr+1]) + 
    (c[p*cstr+2] - c[j*cstr+2])*(c[p*cstr+2] - c[j*cstr+2]) + 
    (c[p*cstr+3] - c[j*cstr+3])*(c[p*cstr+3] - c[j*cstr+3]) ;

  if ( di < dj ) {
    *pc = i ; *pf = j ;
  } else {
    *pc = j ; *pf = i ;
  }
  
  return ;
}

static void insert_point(int *p, int np, int idx, int pi)

{
  int i ;

  for ( i = np ; i >= idx ; i -- ) p[i] = p[i-1] ;

  p[idx] = pi ;
  
  return ;
}

static void delete_point(int *p, int np, int idx)

{
  int i ;

  for ( i = idx ; i < np ; i ++ ) p[i] = p[i+1] ;
  
  return ;
}

static void uv_midpoint(double *uv1, double *uv2, double *uv, int *fix)

{
  int i ;
  double dmax ;

  dmax = 0 ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    uv[i] = 0.5*(uv1[i] + uv2[i]) ;
    if ( ABS(uv1[i] - uv2[i]) > dmax ) {
      dmax = ABS(uv1[i] - uv2[i]) ;
      *fix = i ;
    }
  }
  
  return ;
}

int ssi_curve_join_points(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
			  double *c, int cstr, int *nc, int *dup, int nd,
			  int *p, int *brks, int *nb, ssi_crawl_workspace_t *ws)

{
  int np, i, j, k, n, pc, pf, off, fix ;
  double *x1, *x2, d1, d2, *uv, work[27] ;
  
  brks[0] = 0 ; *nb = 0 ;

  d1 = 1.0 / (1 << bvh1->depth ) ;
  d2 = 1.0 / (1 << bvh2->depth ) ;

  ws->dlim = 2.0*M_SQRT2*d1 ;
  ssi_crawl(ws, &(p[brks[(*nb)]]), &np) ;
  /* brks[1] = np ; (*nb) ++ ; */
  /* return 0 ; */
  
  while ( np > 0 ) {
    off = brks[(*nb)] ;
    for ( i = 0 ; i < np ; i ++ ) {
      if ( point_is_duplicate(dup, nd, p[off+i], &j) ) {
	/*find which point is closer, insert it and add a break*/
	closest_point(c, cstr, p[off+i-1], &(dup[j*2]), &pc, &pf) ;
	p[off+i] = pc ;
	insert_point(p, brks[(*nb)]+np, off+i+1, pf) ; np ++ ;
	dup[2*j+0] = -dup[2*j+0] ; dup[2*j+1] = -dup[2*j+1] ;
	brks[(*nb) + 1] = off + i + 1 ; (*nb) ++ ;
      }      
    }

    brks[(*nb) + 1] = brks[(*nb)] + np ; (*nb) ++ ;
    ssi_crawl(ws, &(p[brks[(*nb)]]), &np) ;
  }

  /* return 0 ; */

  for ( i = 0 ; i < (*nb) ; i ++ ) {
    for ( j = brks[i] ; j < brks[i+1] - 1 ; j ++ ) {
      x1 = &(c[p[j+0]*cstr]) ;
      x2 = &(c[p[j+1]*cstr]) ;
      /* if ( ssi_point_is_boundary(x1) && ssi_point_is_boundary(x2) ) { */
	/*check distance between boundary points*/
	if ( ABS(x1[0] - x2[0]) > d1 ||
	     ABS(x1[1] - x2[1]) > d1
	     /* || */
	     /* ABS(x1[2] - x2[2]) > d2 || */
	     /* ABS(x1[3] - x2[3]) > d2 */
	     ) {
	  uv = &(c[(*nc)*cstr]) ;
	  uv_midpoint(x1, x2, uv, &fix) ;
	  /* if ( refine_intersection(bvh1, bvh2, uv, 1e-12, */
	  /* 			   1 << bvh2->depth, fix, work) == 0 ) { */
	  /*   insert_point(brks, (*nb) + 1, i + 1, j + 1) ; (*nb) ++ ; */
	  /*   insert_point(p, brks[(*nb)], j+1, (*nc)) ; np ++ ; */
	  /*   (*nc) ++ ; */
	  /*   for ( k = i + 1 ; k <= (*nb) ; k ++ ) brks[k] ++ ; */
	  /* } else { */
	  /*   /\*split the curve*\/ */
	  /*   insert_point(brks, (*nb) + 1, i + 1, j + 1) ; (*nb) ++ ; */
	  /* } */
	  /* fix = 1 ; */
	  if ( refine_intersection(bvh1, bvh2, uv, 1e-14,
				   1 << bvh2->depth, fix, work) != 0 ) {
	    /*split the curve*/
	    insert_point(brks, (*nb) + 1, i + 1, j + 1) ; (*nb) ++ ;
	  }
	}
    }
    /* } */
  }
  
  return 0 ;
  
  /*remove single point curves*/
  for ( i = 0 ; i < (*nb) ; i ++ ) {
    if ( brks[i+1] == brks[i] + 1 ) {
      /*remove the single point from the global list of points*/
      delete_point(   p,    brks[(*nb)], brks[i]) ; brks[(*nb)] -- ;
      /*remove the breakpoints for the single point*/
      delete_point(brks, (*nb),       i) ; (*nb) -- ;
      for ( j = i ; j <= (*nb) ; j ++ ) brks[j] -- ;
    }
  }
  
  return 0 ;
}

double ssi_bvh_u_global(ssi_bvh_t *bvh, double ul)

{
  double ug ;
  ssi_mapping_func_t mfunc ;

  mfunc = ssi_bvh_mapping_u(bvh) ;
  if ( mfunc == NULL ) {
    fprintf(stderr, "%s: no mapping defined\n", __FUNCTION__) ;
    exit(1) ;
  }

  /* if ( ul < 0 || ul > 1 ) { */
  /*   fprintf(stderr, "%s: ul (%lg) out of range (0,1)\n", __FUNCTION__, ul) ; */
  /*   exit(1) ; */
  /* } */

  ug = mfunc(ul, ssi_bvh_mapping_data_u(bvh)) ;

  return ug ;
}

double ssi_bvh_v_global(ssi_bvh_t *bvh, double vl)

{
  double vg ;
  ssi_mapping_func_t mfunc ;

  mfunc = ssi_bvh_mapping_v(bvh) ;
  if ( mfunc == NULL ) {
    fprintf(stderr, "%s: no mapping defined\n", __FUNCTION__) ;
    exit(1) ;
  }

  /* if ( vl < 0 || vl > 1 ) { */
  /*   fprintf(stderr, "%s: vl (%lg) out of range (0,1)\n", __FUNCTION__, vl) ; */
  /*   exit(1) ; */
  /* } */
  
  vg = mfunc(vl, ssi_bvh_mapping_data_v(bvh)) ;
  
  return vg ;
}

int ssi_point_is_boundary(double *v)

{
  if ( v[0] == 0 || v[0] == 1 || v[1] == 0 || v[1] == 1 ) return 1 ;
  if ( v[2] == 0 || v[2] == 1 || v[3] == 0 || v[3] == 1 ) return 2 ;
  
  return 0 ;
}

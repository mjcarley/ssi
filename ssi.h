#ifndef __SSI_H_INCLUDED__
#define __SSI_H_INCLUDED__

#include <stdio.h>

#define MIN(_a,_b) ((_a) < (_b) ? (_a) : (_b))
#define MAX(_a,_b) ((_a) > (_b) ? (_a) : (_b))
#define ABS(_a)    (((_a) >= 0) ? (_a) : (-(_a)))

#define SSI_BVH_DEPTH_MAX 16

#define ssi_vector_distance2(_A,_B)		\
  (((_A)[0]-(_B)[0])*((_A)[0]-(_B)[0]) +	\
   ((_A)[1]-(_B)[1])*((_A)[1]-(_B)[1]) +	\
   ((_A)[2]-(_B)[2])*((_A)[2]-(_B)[2]))

#define ssi_vector_distance(_A,_B) (sqrt(ssi_vector_distance2((_A),(_B))))

#define ssi_vector_diff(_A,_B,_C)			\
  do {							\
  (_A)[0] = (_B)[0] - (_C)[0] ;				\
  (_A)[1] = (_B)[1] - (_C)[1] ;				\
  (_A)[2] = (_B)[2] - (_C)[2] ;				\
  } while (0)

#define ssi_vector_shift(_A,_B)			\
  do {						\
  (_A)[0] += (_B)[0] ;				\
  (_A)[1] += (_B)[1] ;				\
  (_A)[2] += (_B)[2] ;				\
  } while (0)

#define ssi_vector_scale(_A,_s)					\
  do {								\
    (_A)[0] *= (_s) ; (_A)[1] *= (_s) ; (_A)[2] *= (_s) ;	\
  } while (0)

#define ssi_vector_scalar(_A,_B)				\
  (((_A)[0])*((_B)[0])+((_A)[1])*((_B)[1])+((_A)[2])*((_B)[2]))

#define ssi_vector_length(_A)		\
  (sqrt(ssi_vector_scalar((_A),(_A))))

#define ssi_vector_cross(_C,_A,_B)				\
  ((_C)[0] = (_A)[1]*(_B)[2] - (_A)[2]*(_B)[1],			\
   (_C)[1] = (_A)[2]*(_B)[0] - (_A)[0]*(_B)[2],			\
   (_C)[2] = (_A)[0]*(_B)[1] - (_A)[1]*(_B)[0])

#define ssi_vector_extend(_x,_o,_r,_t)		\
  do {						\
    (_x)[0] = (_o)[0] + (_t)*(_r[0]) ;		\
    (_x)[1] = (_o)[1] + (_t)*(_r[1]) ;		\
    (_x)[2] = (_o)[1] + (_t)*(_r[2]) ;		\
  } while (0)

typedef int (* ssi_surface_func_t)(double u, double v, double *x,
				   double *xu, double *xv, void *data) ;

typedef double (* ssi_mapping_func_t)(double u, void *data) ;

typedef struct _ssi_bvh_t ssi_bvh_t ;

struct _ssi_bvh_t {
  int depth, off[SSI_BVH_DEPTH_MAX] ;
  double *x, *xb, uvlim[4] ;
  ssi_mapping_func_t umfunc, vmfunc ;
  ssi_surface_func_t sfunc ;
  void *sdata, *umdata, *vmdata ;
} ;

#define ssi_bvh_umin(_b)      ((_b)->uvlim[0])
#define ssi_bvh_umax(_b)      ((_b)->uvlim[1])
#define ssi_bvh_vmin(_b)      ((_b)->uvlim[2])
#define ssi_bvh_vmax(_b)      ((_b)->uvlim[3])
#define ssi_bvh_depth(_b)     ((_b)->depth)
#define ssi_bvh_function(_b)  ((_b)->sfunc)
#define ssi_bvh_data(_b)      ((_b)->sdata)
#define ssi_bvh_mapping_u(_b) ((_b)->umfunc)
#define ssi_bvh_mapping_v(_b) ((_b)->vmfunc)
#define ssi_bvh_mapping_data_u(_b) ((_b)->umdata)
#define ssi_bvh_mapping_data_v(_b) ((_b)->vmdata)

double ssi_mapping_linear(double u, void *data) ;
double ssi_mapping_cosine(double u, void *data) ;
double ssi_mapping_split_cosine(double u, void *data) ;

ssi_bvh_t *ssi_bvh_new(int depth) ;
double ssi_bvh_u_global(ssi_bvh_t *bvh, double ul) ;
double ssi_bvh_v_global(ssi_bvh_t *bvh, double vl) ;
int ssi_bvh_fill_grid(ssi_bvh_t *bvh, ssi_surface_func_t sfunc, void *data) ;
int ssi_bvh_make_bboxes(ssi_bvh_t *bvh, int level, double del) ;
int ssi_bvh_patch_bbox(double *x1, double *x2, double *x3, double *x4,
		       double del, double *b1, double *b2) ;
int ssi_bvh_write_grid(FILE *f, ssi_bvh_t *bvh) ;
int ssi_surface_write_grid(FILE *f, 
			   ssi_surface_func_t sfunc2, void *data,
			   double umin, double umax,
			   double vmin, double vmax,
			   int nu, int nv) ;
int ssi_bvh_write_bboxes(FILE *f, ssi_bvh_t *bvh, int level) ;
int ssi_bbox_ray_slab_test(double *p0, double *p1,
			   double *orig, double *invr,
			   double rtmin, double rtmax) ;
int ssi_bvh_ray_intersection(ssi_bvh_t *bvh,
			     double *orig, double *r,
			     double rtmin, double rtmax,
			     double *uv, int uvstr) ;
int ssi_bvh_ray_intersection_slow(ssi_bvh_t *bvh,
				  double *orig, double *r,
				  double rtmin, double rtmax,
				  double *uv, int uvstr, int *ni) ;
int ssi_surface_surface_intersection(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
				     double *c, int cstr, int *nc,
				     int ncmax, double tol, double *work) ;
int ssi_curve_sort_points(double *c, int cstr, int nc, int *dup, int *nd) ;
int ssi_intersection_curve_order(double *uv, int uvstr, int n,
				 int *brks, int *nb, int nbmax) ;
int ssi_curve_join_points(ssi_bvh_t *bvh1, ssi_bvh_t *bvh2,
			  double *c, int cstr, int *nc, int *dup, int nd,
			  int *p, int *brks, int *nb,
			  ssi_crawl_workspace_t *ws) ;

int ssi_delaunay_triangle(double *pts, int npts, int *e, int *ne) ;
int ssi_delaunay(double *c, int cstr, int nc, int *e, int *ne, double *dwork) ;
int ssi_point_is_boundary(double *v) ;

#endif /*__SSI_H_INCLUDED__*/

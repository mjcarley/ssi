#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>

#include "bilinear.h"
#include "crawl.h"
#include "ssi.h"

int sphere_func(double u, double v, double *x, double *xu, double *xv,
		void *data) ;
int cylinder_func(double u, double v, double *x, double *xu, double *xv,
		  void *data) ;
int wing_func(double u, double v, double *x, double *xu, double *xv,
	      void *data) ;

char *progname ;
char *surface_fnames[] = {"sphere", "cylinder", "wing", NULL} ;
void *surface_funcs[]  = {sphere_func, cylinder_func, wing_func, NULL} ;
int   surface_nargs[]  = {4, 5, 5, -1} ;

int sphere_func(double u, double v, double *x, double *xu, double *xv,
		void *data)

{
  double sgn, *c = data ;
  double V, Vv ;
  
  sgn = ( u >= 0 ? 1 : -1 ) ;
  u = fabs(u) ;
  V = sqrt(v)*sqrt(1-v) ;
  x[0] =     2*(u-0.5)*V + 0.5 ;
  x[1] = sgn*2*sqrt(u)*sqrt(1-u)*V ;
  x[2] =     v ;

  if ( c == NULL ) return 0 ;

  x[0] -= 0.5 ; x[2] -= 0.5 ;

  ssi_vector_scale(x,c[3]) ;
  ssi_vector_shift(x, c) ;

  if ( xu == NULL || xv == NULL ) return 0 ;

  Vv = 0.5*sqrt(1-v)/sqrt(v) - 0.5*sqrt(v)/sqrt(1-v) ;

  xv[0] =     2*(u-0.5)*Vv ;
  xv[1] = sgn*2*sqrt(u)*sqrt(1-u)*Vv ;
  xv[2] =     1.0 ;

  xv[0] *= c[3] ; xv[1] *= c[3] ; xv[2] *= c[3] ; 

  xu[0] = sgn*2*V ;
  if ( sgn >= 0 ) {
    xu[1] =  (sqrt(1-u)/sqrt(u) - sqrt(u)/sqrt(1-u))*V ;
  } else {
    xu[1] = -(sqrt(1-u)/sqrt(u) + sqrt(u)/sqrt(1-u))*V ;
  }
  xu[2] = 0.0 ;
  
  xu[0] *= c[3] ; xu[1] *= c[3] ; xu[2] *= c[3] ; 

  return 0 ;
}

int cylinder_func(double u, double v, double *x, double *xu, double *xv,
		  void *data)

{
  double sgn, *c = data ;

  sgn = ( u >= 0 ? 1 : -1 ) ;
  u = fabs(u) ;
  x[0] =     2*(u-0.5) ;
  x[1] = sgn*2*sqrt(u)*sqrt(1-u) ;
  x[2] =     v ;

  if ( c == NULL ) return 0 ;

  x[0] *= c[3] ; x[1] *= c[3] ; x[2] *= c[4] ;

  ssi_vector_shift(x, c) ;
  
  if ( xu == NULL || xv == NULL ) return 0 ;
  
  return 0 ;
}

int wing_func(double u, double v, double *x, double *xu, double *xv,
	      void *data)

{
  double sgn, *c = data ;

  sgn = ( u >= 0 ? 1 : -1 ) ;
  u = fabs(u) ;
  x[2] =     2*(u-0.5) ;
  x[0] = sgn*2*sqrt(u)*(1-u)*0.2 + 0.0*u ;
  x[1] =     v ;

  if ( c == NULL ) return 0 ;

  x[2] *= c[3] ; x[0] *= c[3] ; x[1] *= c[4] ;

  ssi_vector_shift(x, c) ;
  
  return 0 ;
}

static void parse_func(void *func, char **name, int *na)

{
  int i ;
  *name = NULL ; *na = 0 ;

  for ( i = 0 ; surface_fnames[i] != NULL ; i ++ ) {
    if ( surface_funcs[i] == func ) {
      *name = surface_fnames[i] ; *na = surface_nargs[i] ;
      return ;
    }
  }
  
  return ;
}

static void parse_surface(char *arg, void **func, double *c)

{
  char *tok ;
  int i, na ;
  
  tok = strtok(arg, "(), ") ;

  *func = NULL ;
  for ( i = 0 ; surface_fnames[i] != NULL ; i ++ ) {
    if ( strcmp(tok, surface_fnames[i]) == 0 ) {
      *func = surface_funcs[i] ; na = surface_nargs[i] ;
      break ;
    }
  }

  if ( *func == NULL ) {
    fprintf(stderr, "%s: unrecognised surface `%s'\n", progname, arg) ;
    exit(1) ;
  }

  i = 0 ;
  tok = strtok(NULL, "(), ") ;    
  while ( tok != NULL ) {
    /* fprintf(stderr, "%s\n", tok) ; */
    c[i] = atof(tok) ; i ++ ;
    tok = strtok(NULL, "(), ") ;    
  } 

  if ( i != na ) {
    fprintf(stderr,
	    "%s: incorrect number (%d) of parameters for surface\n",
	    progname, i) ;
    exit(1) ;
  }
  
  return ;
}

static double distance4(double *x, double *y)

{
  double d ;

  d =
    (x[0] - y[0])*(x[0] - y[0]) +
    (x[1] - y[1])*(x[1] - y[1]) +
    (x[2] - y[2])*(x[2] - y[2]) +
    (x[3] - y[3])*(x[3] - y[3]) ;    
  
  return d ;
}

static void list_surfaces(void)

{
  fprintf(stderr,
	  "Surface functions:\n"
	  "  sphere(x,y,z,D)\n"
	  "  cylinder(x,y,z,r,L)\n"
	  "  wing(x,y,z,c,L)\n\n"
	  "  (x,y,z) coordinates of centre of sphere, centre of cylinder end,\n"
	  "          or wing section leading edge;\n"
	  "  D       sphere diameter\n"
	  "  L       cylinder or wing length\n"
	  "  c       wing chord\n") ;
  
  return ;
}

int main(int argc, char **argv)

{
  ssi_bvh_t *bvh1, *bvh2 ;
  ssi_crawl_workspace_t *ws ;
  int depth1, depth2, i, k, cstr, nc, ncmax, nd, p[8192], na ;
  int dup[64], brks[64], nb, j0min, j1min, k0min, k1min ;
  double del, u1, v1, u2, v2, c1[8], c2[8], *p0, *p1, *q0, *q1 ;
  double *c, x1[3], x2[3], work[8192], tol, *pm, d0min, d1min, d0, d1 ;
  void *sfunc1, *sfunc2 ;
  FILE *output ;
  char ch, *sname ;
  
  depth1 = depth2 = 6 ; del = 1e-6 ; tol = 1e-12 ;

  ncmax = 8192 ; cstr = 4 ;
  progname = argv[0] ;

  /*defaults for basic testing*/
  /*sphere parameters (centre, radius)*/
  sfunc1 = sphere_func ;
  c1[0] =  0.3 ; c1[1] = 0.4 ; c1[2] = -0.1 ; c1[3] = 1.4 ;
  /*cylinder parameters (origin,radius,height)*/
  sfunc2 = cylinder_func ;
  c2[0] = 0.4 ; c2[1] = 0.1 ; c2[2] = -1.1 ; c2[3] = 0.5 ; c2[4] = 2.0 ;

  while ( (ch = getopt(argc, argv, "D:d:lS:s:t:")) != EOF ) {
    switch ( ch ) {
    default: fprintf(stderr,
		     "%s: unrecognised option `%c'\n", progname, ch) ;
      exit(1) ;
      break ;
    case 'D': depth1 = atoi(optarg) ; break ;
    case 'd': depth2 = atoi(optarg) ; break ;
    case 'l': list_surfaces() ; return 0 ; break ;
    case 'S': parse_surface(optarg, &sfunc1, c1) ; break ;
    case 's': parse_surface(optarg, &sfunc2, c2) ; break ;
    case 't': tol = atof(optarg) ; break ;
    }
  }

  c = (double *)malloc(ncmax*cstr*sizeof(double)) ;

  ws = ssi_crawl_workspace_alloc(ncmax, 2*ncmax) ;
    
  bvh1 = ssi_bvh_new(depth1) ; bvh2 = ssi_bvh_new(depth2) ;
  
  ssi_bvh_umin(bvh1) = -1 ; ssi_bvh_umax(bvh1) =  1 ;
  ssi_bvh_mapping_u(bvh1) = ssi_mapping_split_cosine ;
  ssi_bvh_mapping_v(bvh1) = ssi_mapping_cosine ;

  ssi_bvh_umin(bvh2) = -1 ; ssi_bvh_umax(bvh2) =  1 ;
  ssi_bvh_mapping_u(bvh2) = ssi_mapping_split_cosine ;
  ssi_bvh_mapping_v(bvh2) = ssi_mapping_cosine ;

  ssi_bvh_fill_grid(bvh1, sfunc1, c1) ;
  ssi_bvh_fill_grid(bvh2, sfunc2, c2) ;
  
  parse_func(ssi_bvh_function(bvh1), &sname, &na) ; pm = ssi_bvh_data(bvh1) ;
  fprintf(stderr, "surface 1: %s(", sname) ;
  for ( i = 0 ; i < na-1 ; i ++ ) fprintf(stderr, "%lg,", pm[i]) ;
  fprintf(stderr, "%lg)\n", pm[na-1]) ;
  fprintf(stderr, "    depth: %d\n", ssi_bvh_depth(bvh1)) ;
  parse_func(ssi_bvh_function(bvh2), &sname, &na) ; pm = ssi_bvh_data(bvh2) ;
  fprintf(stderr, "surface 2: %s(", sname) ;
  for ( i = 0 ; i < na-1 ; i ++ ) fprintf(stderr, "%lg,", pm[i]) ;
  fprintf(stderr, "%lg)\n", pm[na-1]) ;
  fprintf(stderr, "    depth: %d\n", ssi_bvh_depth(bvh2)) ;

  fprintf(stderr, "intersection tolerance %lg\n", tol) ;
  
  for ( i = ssi_bvh_depth(bvh1) ; i >= 0 ; i -- ) {
    ssi_bvh_make_bboxes(bvh1, i, del) ;
  }
  for ( i = ssi_bvh_depth(bvh2) ; i >= 0 ; i -- ) {
    ssi_bvh_make_bboxes(bvh2, i, del) ;
  }

  output = fopen("bvh1.dat", "w") ;
  ssi_bvh_write_grid(output, bvh1) ;
  fclose(output) ;

  output = fopen("bvh2.dat", "w") ;
  ssi_bvh_write_grid(output, bvh2) ;
  fclose(output) ;

  output = fopen("boxes.dat", "w") ;
  ssi_bvh_write_bboxes(output, bvh1, ssi_bvh_depth(bvh1)-1) ;
  fclose(output) ;

  nc = 0 ; nd = 0 ;
  ssi_surface_surface_intersection(bvh1, bvh2, c, cstr, &nc, ncmax, tol, work) ;
  ssi_curve_sort_points(c, cstr, nc, dup, &nd) ;

  /* for ( i = 0 ; i < nd ; i ++ ) { */
  /*   fprintf(stderr, "%d %d\n", dup[i*2+0], dup[i*2+1]) ; */
  /* } */
  
  fprintf(stderr, "%d intersection points\n", nc) ;

  if ( nc == 0 ) return 0 ;
  
  /* fprintf(stderr, "%d duplicate points\n", nd) ; */
  
  ssi_delaunay(c, cstr, nc - nd , ws->e, &(ws->ne), work) ;
  ssi_crawl_init(c, cstr, nc - nd , ws) ;

  ssi_curve_join_points(bvh1, bvh2, c, cstr, &nc, dup, nd, p, brks, &nb, ws) ;

  output = fopen("points.dat", "w") ;

  for ( i = 0 ; i < nc ; i ++ ) {
    u1 = ssi_bvh_u_global(bvh1,c[i*cstr+0]) ;
    v1 = ssi_bvh_v_global(bvh1,c[i*cstr+1]) ;
    ssi_bvh_function(bvh1)(u1, v1, x1, NULL, NULL, ssi_bvh_data(bvh1)) ;
    u1 = c[i*cstr+0] ; v1 = c[i*cstr+1] ; 
    u2 = ssi_bvh_u_global(bvh2,c[i*cstr+2]) ;
    v2 = ssi_bvh_v_global(bvh2,c[i*cstr+3]) ;
    ssi_bvh_function(bvh2)(u2, v2, x2, NULL, NULL, ssi_bvh_data(bvh2)) ;
    u2 = c[i*cstr+2] ; v2 = c[i*cstr+3] ; 
    fprintf(output, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e "
	    "%1.16e %1.16e %1.16e %1.16e\n",
	    x1[0], x1[1], x1[2], x2[0], x2[1], x2[2], u1, v1, u2, v2) ;
  }
  
  fclose(output) ;

  /*nearest points to curves*/
  for ( k = 0 ; k < nb ; k ++ ) {
    int i0, i1 ;
    i0 = p[brks[k]] ; i1 = p[brks[k+1]-1] ;
    p0 = &(c[i0*cstr]) ; p1 = &(c[i1*cstr]) ;
    d0min = d1min = DBL_MAX ;
    for ( i = 0 ; i < k ; i ++ ) {
      int j0, j1 ;
      j0 = p[brks[i]] ; j1 = p[brks[i+1]-1] ;
      q0 = &(c[j0*cstr]) ; q1 = &(c[j1*cstr]) ;
      d0 = distance4(p0,q0) ; if ( d0 < d0min ) {
	d0min = d0 ; j0min = j0 ; k0min = i ;
      }
      d0 = distance4(p0,q1) ; if ( d0 < d0min ) {
	d0min = d0 ; j0min = j1 ; k0min = i ;
      }
      d1 = distance4(p1,q0) ; if ( d1 < d1min ) {
	d1min = d0 ; j1min = j0 ; k1min = i ;
      }
      d1 = distance4(p1,q1) ; if ( d1 < d1min ) {
	d1min = d0 ; j1min = j1 ; k1min = i ;
      }
    }
    for ( i = k+1 ; i < nb ; i ++ ) {
      int j0, j1 ;
      j0 = p[brks[i]] ; j1 = p[brks[i+1]-1] ;
      q0 = &(c[j0*cstr]) ; q1 = &(c[j1*cstr]) ;
      d0 = distance4(p0,q0) ; if ( d0 < d0min ) {
	d0min = d0 ; j0min = j0 ; k0min = i ;
      }
      d0 = distance4(p0,q1) ; if ( d0 < d0min ) {
	d0min = d0 ; j0min = j1 ; k0min = i ;
      }
      d1 = distance4(p1,q0) ; if ( d1 < d1min ) {
	d1min = d0 ; j1min = j0 ; k1min = i ;
      }
      d1 = distance4(p1,q1) ; if ( d1 < d1min ) {
	d1min = d0 ; j1min = j1 ; k1min = i ;
      }
    }

    fprintf(stderr, "curve %d: %d %d (%lg) %d %d (%lg)\n",
	    k, k0min, j0min, sqrt(d0min), k1min, j1min, sqrt(d1min)) ;    
  }
  
  fprintf(stderr, "%d distinct curves\n", nb) ;
  
  output = fopen("curves.dat", "w") ;
  for ( k = 0 ; k < nb ; k ++ ) {
    for ( i = brks[k] ; i < brks[k+1] ; i ++ ) {
      fprintf(output, "%d %d\n", p[i], k) ;
    }
  }  
  fclose(output) ;
  
  return 0 ;
}

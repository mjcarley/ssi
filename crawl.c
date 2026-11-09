#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "crawl.h"

#define MIN(_a,_b) ((_a) < (_b) ? (_a) : (_b))
#define MAX(_a,_b) ((_a) > (_b) ? (_a) : (_b))
#define ABS(_a)    (((_a) >= 0) ? (_a) : (-(_a)))

static int point_is_boundary(double *v, int vstr, int i)

{
  double *x = &(v[i*vstr]) ;

  if ( x[0] == 0 || x[0] == 1 ) return 1 ;
  if ( x[1] == 0 || x[1] == 1 ) return 1 ;
  
  return 0 ;
}

ssi_crawl_workspace_t *ssi_crawl_workspace_alloc(int nvmax, int nemax)

{
  ssi_crawl_workspace_t *w ;

  w = (ssi_crawl_workspace_t *)malloc(sizeof(ssi_crawl_workspace_t)) ;

  w->e  = (int *)malloc((2*nemax + nvmax)*sizeof(int)) ;
  w->vt = &(w->e[2*nemax]) ;

  ssi_crawl_workspace_edge_number_max(w) = nemax ;
  ssi_crawl_workspace_vertex_number_max(w) = nvmax ;
  
  w->elen = (double *)malloc(nemax*sizeof(double)) ;

  w->dlim = 0.0 ;
  
  return w ;
}

int ssi_crawl_init(double *v, int vstr, int nv, ssi_crawl_workspace_t *w)

{
  int i, *e, j ;
  double *x1, *x2 ;

  if ( nv > ssi_crawl_workspace_vertex_number_max(w) ) {
    fprintf(stderr, "%s: not enough space (%d) for %d vertices\n",
	    __FUNCTION__, ssi_crawl_workspace_vertex_number_max(w), nv) ;
    exit(1) ;
  }
  
  for ( i = 0 ; i < ssi_crawl_workspace_edge_number(w) ; i ++ ) {
    e = ssi_crawl_workspace_edge(w,i) ;
    if ( e[0] > e[1] ) {
      j = e[0] ; e[0] = e[1] ; e[1] = j ;
    }
    x1 = &(v[e[0]*vstr]) ; x2 = &(v[e[1]*vstr]) ;
    w->elen[i] = sqrt((x1[0] - x2[0])*(x1[0] - x2[0]) +
		      (x1[1] - x2[1])*(x1[1] - x2[1])) ;
  }

  ssi_crawl_workspace_vertex_number(w) = nv ;
  for ( i = 0 ; i < ssi_crawl_workspace_vertex_number(w) ; i ++ ) {
    ssi_crawl_workspace_vertex_tag(w,i) = 0 ;
  }

  w->v = v ; w->vstr = vstr ;
  
  return 0 ;
}  

static void insert_vertex(int *p, int *np, int v, int idx, int *tags)

{
  int i ;
  
  tags[v] = 1 ;
  if ( idx == (*np) ) {
    p[(*np)] = v ; (*np) ++ ;
    return ;
  }

  for ( i = (*np) ; i > idx ; i -- ) {
    p[i] = p[i-1] ;
  }

  p[idx] = v ; (*np) ++ ;
  
  return ;
}

int ssi_crawl(ssi_crawl_workspace_t *w, int *p, int *np)

{
  int i, imin, *e, added, vstr ;
  double dmin, d1, *x1, *x2, *v ;

  (*np) = 0 ; v = w->v ; vstr = w->vstr ;
  /*find shortest untagged edge*/
  dmin = DBL_MAX ; imin = -1 ;
  for ( i = 0 ; i < ssi_crawl_workspace_edge_number(w) ; i ++ ) {
    if ( !ssi_crawl_workspace_edge_tag(w,i) && w->elen[i] < dmin ) {
      dmin = w->elen[i] ; imin = i ;
    }
  }

  /*no untagged edges left*/
  if ( imin == -1 ) return 0 ;

  e = ssi_crawl_workspace_edge(w,imin) ;
  /*add the first edge to the list*/
  p[(*np)] = e[0] ; ssi_crawl_workspace_vertex_tag(w,p[(*np)]) = 1 ; (*np) ++ ;
  p[(*np)] = e[1] ; ssi_crawl_workspace_vertex_tag(w,p[(*np)]) = 1 ; (*np) ++ ;

  added = 1 ; d1 = DBL_MAX ;

  while ( added == 1 ) {
    added = 0 ;
    if ( (*np) > 4 ) {
      x1 = &(v[p[0]*vstr]) ; x2 = &(v[p[(*np)-1]*vstr]) ;
      d1 = sqrt((x1[0] - x2[0])*(x1[0] - x2[0]) +
		(x1[1] - x2[1])*(x1[1] - x2[1])) ;
    }

    /* if ( w->dlim != 0 ) d1 = MIN(d1, w->dlim) ; */
    
    /* d1 = MIN(d1, (1.42/(1 << 5))) ; */
    
    /*find the shortest edge connected to an endpoint and not already
      tagged*/
    dmin = DBL_MAX ; imin = -1 ;
    for ( i = 0 ; i < ssi_crawl_workspace_edge_number(w) ; i ++ ) {
      e = ssi_crawl_workspace_edge(w,i) ;
      if ( (e[0] == p[0]       || e[1] == p[0] ||
	    e[0] == p[(*np)-1] || e[1] == p[(*np)-1]) &&
	   !ssi_crawl_workspace_edge_tag(w,i) &&
	   !(point_is_boundary(v,vstr,e[0]) && w->vt[e[0]]) &&
	   !(point_is_boundary(v,vstr,e[1]) && w->vt[e[1]]) &&
	   w->elen[i] < dmin ) {
	dmin = w->elen[i] ; imin = i ;
      }
    }
    if ( imin == -1 ) {
      /*no candidate edge, terminate*/
      if ( *np < 4 ) return 0 ;
      
      /*find second longest edge connected to an extremity point*/
      dmin = DBL_MAX ; 
      for ( i = 0 ; i < ssi_crawl_workspace_edge_number(w) ; i ++ ) {
	if ( !ssi_crawl_workspace_edge_tag(w,i) ) {
	  /*untagged edge: check it*/
	  e = ssi_crawl_workspace_edge(w,i) ;
	  if ( e[0] == p[0] || e[0] == p[(*np)-1] ||
	       e[1] == p[0] || e[1] == p[(*np)-1] ) {
	    /*edge connected to an extremity: is it shorter than d1?*/
	    if ( d1 > w->elen[i] ) return 0 ;
	  }
	}
      }

      /*check for boundary hit*/
      x1 = &(v[p[(*np)-1]*vstr]) ;
      if ( x1[0] == 0 || x1[0] == 1.0 || x1[1] == 0.0 || x1[1] == 1.0 ||
	   x1[2] == 0 || x1[2] == 1.0 || x1[3] == 0.0 || x1[3] == 1.0 )
	return 0 ;
      x1 = &(v[p[0]*vstr]) ;
      if ( x1[0] == 0 || x1[0] == 1.0 || x1[1] == 0.0 || x1[1] == 1.0 ||
	   x1[2] == 0 || x1[2] == 1.0 || x1[3] == 0.0 || x1[3] == 1.0 )
	return 0 ;
      
      p[(*np)] = p[0] ; (*np) ++ ;
      
      return 0 ;
    }

    e = ssi_crawl_workspace_edge(w,imin) ;
    if ( w->elen[imin] < d1 ) {
      if ( e[0] == p[0] ) {
	insert_vertex(p, np, e[1], 0, w->vt) ;
      } else 
	if ( e[1] == p[0] ) {
	  insert_vertex(p, np, e[0], 0, w->vt) ;
	} else 
	  if ( e[0] == p[(*np)-1] ) {
	    insert_vertex(p, np, e[1], (*np), w->vt) ;
	  } else 
	    if ( e[1] == p[(*np)-1] ) {
	      insert_vertex(p, np, e[0], (*np), w->vt) ;
	    }
      added = 1 ;
    }    
  }
  
  return 0 ;
}

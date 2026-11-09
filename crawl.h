#ifndef __CRAWL_H_INCLUDED__
#define __CRAWL_H_INCLUDED__

typedef struct _ssi_crawl_workspace_t ssi_crawl_workspace_t ;

struct _ssi_crawl_workspace_t {
  int *e, ne, nemax, *vt, nv, nvmax, vstr ;
  double *v, *elen, dlim ;
} ;

#define ssi_crawl_workspace_edge(_w,_i)            (&((_w)->e[2*(_i)]))
#define ssi_crawl_workspace_vertex_tag(_w,_i)      ((_w)->vt[(_i)])
#define ssi_crawl_workspace_edge_number(_w)        ((_w)->ne)
#define ssi_crawl_workspace_edge_number_max(_w)    ((_w)->nemax)
#define ssi_crawl_workspace_vertex_number(_w)      ((_w)->nv)
#define ssi_crawl_workspace_vertex_number_max(_w)  ((_w)->nvmax)
#define ssi_crawl_workspace_edge_tag(_w,_i)				\
  ((ssi_crawl_workspace_vertex_tag((_w),(_w)->e[2*(_i)+0]) == 1) &&	\
   (ssi_crawl_workspace_vertex_tag((_w),(_w)->e[2*(_i)+1]) == 1))

ssi_crawl_workspace_t *ssi_crawl_workspace_alloc(int nvmax, int nemax) ;
int ssi_crawl_init(double *v, int vstr, int nv, ssi_crawl_workspace_t *w) ;
int ssi_crawl(ssi_crawl_workspace_t *w, int *p, int *np) ;

#endif /*__CRAWL_H_INCLUDED__*/

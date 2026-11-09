/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2023 Michael Carley
 *
 * AGG is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  AGG is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AGG.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <triangle.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include <triangle_api.h>
#pragma GCC diagnostic pop

extern void traversalinit(struct memorypool *pool);
extern void *traverse(struct memorypool *pool);
extern vertex vertextraverse(mesh *m);
extern triangle *triangletraverse(mesh *m);
extern int plus1mod3[] ;
extern int minus1mod3[] ;

static void reset_triangleio(triangleio *io)

{
  io->pointlist = (double *) NULL;
  io->pointattributelist = (double *) NULL;
  io->pointmarkerlist = (int *) NULL;
  io->numberofpoints = 0;
  io->numberofpointattributes = 0;
  
  io->trianglelist = (int *) NULL;
  io->triangleattributelist = (double *) NULL;
  io->trianglearealist = (double *) NULL;
  io->neighborlist = (int *) NULL;
  io->numberoftriangles = 0;
  io->numberofcorners = 0;
  io->numberoftriangleattributes = 0;
  
  io->segmentlist = (int *) NULL;
  io->segmentmarkerlist = (int *) NULL;
  io->numberofsegments = 0;
  
  io->holelist = (double *) NULL;
  io->numberofholes = 0;
  io->regionlist = (double *) NULL;
  io->numberofregions = 0;
  
  io->edgelist = (int *) NULL;
  io->edgemarkerlist = (int *) NULL;
  io->numberofedges = 0;

  return ;
}

int ssi_delaunay_triangle(double *pts, int npts, int *e, int *ne)

{
  char *args = "ze" ;
  context *ctx ;
  triangleio in ;
  /*it has to be called this because of how the macros in triangle.h
    are defined*/
  mesh *m ;
  behavior *b ;
  statistics s ;
  int i0, i1, vertexnumber ;
  struct otri triangleloop, trisym;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
  struct osub checkmark;
#pragma GCC diagnostic pop
  vertex p1, p2, vertexloop ;
  long edgenumber ;
  triangle ptr;
  subseg sptr;
  /* int tristatus ; */

  reset_triangleio(&in);
  in.pointlist = pts ;
  in.numberofpoints = npts ;

  ctx = triangle_context_create() ;
  triangle_context_options(ctx, args) ;
  
  triangle_mesh_create(ctx, &in) ;
  triangle_mesh_statistics(ctx, &s);
  m = ctx->m ;
  b = ctx->b ;

  /*
   * code adapted from Christian Woltering's Triangle API library:
   * https://github.com/wo80/Triangle/
   */
  traversalinit(&m->vertices);
  vertexnumber = b->firstnumber;
  vertexloop = vertextraverse(m);
  while (vertexloop != (vertex) NULL) {
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      setvertexmark(vertexloop, vertexnumber) ;
      vertexnumber++ ;
    }
    vertexloop = vertextraverse(m);
  }

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  edgenumber = b->firstnumber;
  *ne = 0 ;
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each triangle.  If there isn't another triangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent triangle, operate on the edge only if the current triangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (triangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
	 triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) {
	tspivot(triangleloop, checkmark);
	org(triangleloop, p1);
	dest(triangleloop, p2);
	i0 = vertexmark(p1) ; i1 = vertexmark(p2) ;
	e[2*(*ne)+0] = i0 ; e[2*(*ne)+1] = i1 ; 
	edgenumber++; (*ne) ++ ;
      }
    }
    triangleloop.tri = triangletraverse(m);
  }
  
  triangle_context_destroy(ctx) ;

  return 0 ;
}
			  
int ssi_delaunay(double *c, int cstr, int nc, int *e, int *ne, double *dwork)

{
  int i ;
  double *pts ;

  pts = dwork ;
  for ( i = 0 ; i < nc ; i ++ ) {
    pts[2*i+0] = c[cstr*i+0] ; pts[2*i+1] = c[cstr*i+1] ;
  }

  ssi_delaunay_triangle(pts, nc, e, ne) ;
  
  return 0 ;
}

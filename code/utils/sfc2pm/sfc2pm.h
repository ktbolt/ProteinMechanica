
#include <stdio.h>

typedef struct GeomPolygonEdge {
  int count;
  int poly;
  int node1, node2;
  struct GeomPolygonEdge *next;
  } GeomPolygonEdge;

typedef  float  (DmPoint3)[3];

#define dm_MaxVal(vmax, val)    \
          ((val) > (vmax) ? (val) : (vmax))

#define dm_MinVal(vmin, val)    \
          ((val) < (vmin) ? (val) : (vmin))

void
GeomIndexPolygonsReOrient (int p, int num_polys, int *conn, int hgeom, int n, 
                           int *poly_table, GeomPolygonEdge **edge_table);

void 
GeomIndexPolygonsOrientTrav (int p, int num_polys, int *conn, 
                             GeomPolygonEdge **edge_table,
                             int hgeom, int n, int *poly_table, int poly_proc[]);

void
PolyReorient (int num_polys, int *conn, int num_verts, DmPoint3 *verts);


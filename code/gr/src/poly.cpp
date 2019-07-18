
/* -------------------------------------------------------------------------- *
 *                              Protein Mechanica                             *
 * -------------------------------------------------------------------------- *
 * This is part of the Protein Mechanica coarse-grained molecular motor       *
 * modeling application originating from Simbios, the NIH National Center for *
 * Physics-Based Simulation of Biological Structures at Stanford, funded      *
 * under the NIH Roadmap for Medical Research, grant U54 GM072970.            *
 * See https://simtk.org.                                                     *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: David Parker                                                      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//*============================================================*
//* poly:                 p o l y g o n                        *
//*============================================================*

#include "pm/gr/poly.h"
#include "gc.h"

namespace PmGraphics {

#define ndbg_GrPolygon

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pick                      ==========*
//*============================================================*
// pick a polygon geometry.

void
GrPolygon::pick (GrPickResult& pick)
  {
  #ifdef dbg_GrPolygon
  fprintf (stderr, "\n>>>>>> GrPolygon::pick: \n");
  #endif

  if (!pick_enabled) {
    return;
    }

  if (connectivity) {
    pick_indexed (pick);
    return;
    }

  GrVector3 v, ldir, ipt, min_ipt, p1, p2, lp1, lp2;
  GrVector3 *norms; 
  int n;

  p1 = pick.line[0];
  p2 = pick.line[1];
  ldir = p1 - p2;

  bool ipt_found = false, intersect;
  float tol = pick.atts.tol;
  float min_dist = 1e6;
  float r, dist;
  this->picked_entity = -1;

  // compute face normals //
  getFaceNormals (&norms);

  // compute intersection of polygon with line //
  gc_PolyLineInt (num_vertices, vertices, norms[0], pick.line, tol, intersect, ipt);

  if (!intersect) {
    return;
    }

  picked_entity = 1;

  GrPickGeom pg;
  pg.geom = this;
  pg.int_pt = min_ipt;
  pg.dist = min_dist;
  pg.intersected = false;
  pick.geom_list.push_back(pg);
  }

//*============================================================*
//*==========              pick_indexed              ==========*
//*============================================================*
// pick indexed polygons.

void
GrPolygon::pick_indexed (GrPickResult& pick)
  {
  #ifdef dbg_GrPolygon
  fprintf (stderr, "\n>>>>>> GrPolygon::pick_indexed: \n");
  #endif

  GrVector3 ipt, min_ipt, pts[20], *norms, line[2]; 
  int i, j, k;

  bool ipt_found = false, intersect;
  float tol = pick.atts.tol;
  float min_dist = 1e6;
  float dist;
  this->picked_entity = -1;
  getFaceNormals (&norms);
  int *conn = connectivity;
  int n = hgeom;

  if (xform.set) {
    PmMatrix3x3 mat = xform.matrix;
    mat.transpose();
    line[0] = mat*(pick.line[0]-xform.center - xform.translation) + xform.center;
    line[1] = mat*(pick.line[1]-xform.center - xform.translation) + xform.center;
    }
  else {
    line[0] = pick.line[0];
    line[1] = pick.line[1];
    }

  for (i = 0; i < number; i++) {
    if (!hgeom) {
      n = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      k = *conn;
      pts[j] = vertices[k];
      conn++;
      }

    gc_PolyLineInt (n, pts, norms[i], line, tol, intersect, ipt);

    if (intersect) {
      picked_entity = i+1;

      if (xform.set) {
        PmMatrix3x3 mat = xform.matrix;
        PmVector3 tpt;
        tpt = mat*(ipt-xform.center) + xform.translation + xform.center;
        ipt = tpt;
        }

      GrVector3 v = pick.line[0] - ipt; 
      GrPickGeom pg;
      pg.geom = this;

      if (pick.atts.type == GR_PICK_VERTEX) {
        fprintf (stderr, "   >>> pick vertex \n");
        v = pts[0] - ipt;
        float dmin = v.length();
        float d;
        int jmin = 0;

        for (j = 1; j < n; j++) {
          v = pts[j] - ipt;
          d = v.length();

          if (d < dmin) {
            dmin = d;
            jmin = j;
            }
          }

        pg.int_pt = pts[jmin];
        }
      else {
        pg.int_pt = ipt;
        }

      pg.dist = v.length();
      pg.intersected = true;
      pick.geom_list.push_back(pg);
      }
    }
  }

//*============================================================*
//*==========              render                    ==========*
//*============================================================*
// render a polygon 

void
GrPolygon::render() {

  #ifdef dbg_GrPolygon
  fprintf (stderr, ">>>>>>> GrPolygon::render:  conn [%x] \n", this->connectivity);
  fprintf (stderr, "   >>>> hgeom [%d] \n", hgeom);
  #endif

  if (!active) {
    return;
    }

  api_setup();

  if (connectivity) {
    api_render_indexed ();
    }

  else if (outline) {
    api_render_outline ();
    }

  else {
    api_render ();
    }

  api_reset_state();
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========                init                    ==========*
//*============================================================*
// initialize attributes.

void
GrPolygon::init(int has_conn)
  {
  #ifdef dbg_GrPolygon
  fprintf (stderr, ">>>>>>> GrPolygon::init   has conn [%d] \n", has_conn);
  fprintf (stderr, "   >>>> hgeom [%d] \n", hgeom);
  #endif
  active = true;
  outline = false;
  face_normals = NULL;
  vertex_normals = NULL;
  line_width = 1.0;
  marker_size = 1.0;
  line_style = GR_LINE_STYLE_SOLID;
  display = GR_GEOMETRY_DISPLAY_SOLID;
  shading = GR_GEOMETRY_SHADING_FLAT;
  pick_enabled = true;
  picked_entity = -1;
  }

//*============================================================*
//*==========              getVertexNormals          ==========*
//*============================================================*
// get vertex normals.             

void
GrPolygon::getVertexNormals(GrVector3 **norms)
  {

  *norms = NULL;

  if (vertex_normals) {
    *norms = vertex_normals;
    return;
    }

  GrVector3 *fnorms, norm;
  int i, j;
  int ptr, next_ptr; 
  int node, face;
  int *conn = connectivity;

  // get polygon face normals //

  getFaceNormals (&fnorms);

  vector<int> edge_table(num_vertices, -1);

  // compute face adjacency list //

  int n = hgeom;
  size = number;
  vector<int> face_list(size);
  vector<int> link(size); 
  next_ptr = 0;

  for (i = 0; i < number; i++) {
    if (!hgeom) {
      n = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      node = *conn;
      conn++;

      if (edge_table[node] == -1) {
        ptr = next_ptr++;

        if (ptr == size) {
          size += 1000;
          face_list.resize(size); 
          link.resize(size); 
          }

        edge_table[node] = ptr;
        face_list[ptr] = i;
        link[ptr] = -1;
        }
      else {
        ptr = next_ptr++;

        if (ptr == size) {
          size += 1000;
          face_list.resize(size); 
          link.resize(size); 
          }

        link[ptr] = edge_table[node];
        edge_table[node] = ptr;
        face_list[ptr] = i;
        }
      }
    }


  // average face normals

  vertex_normals = new GrVector3[num_vertices];

  for (i = 0; i < num_vertices; i++) {
    ptr = edge_table[i];
    n = 0;
    norm.set(0.0, 0.0, 0.0);

    while (ptr != -1) {
      face = face_list[ptr];
      norm = norm + fnorms[face];
      ptr = link[ptr];
      n += 1;
      }

    vertex_normals[i][0] = norm[0] / n;
    vertex_normals[i][1] = norm[1] / n;
    vertex_normals[i][2] = norm[2] / n;
    }

  *norms = vertex_normals;
  }

//*============================================================*
//*==========              getFaceNormals            ==========*
//*============================================================*
// get vertex normals.

void
GrPolygon::getFaceNormals(GrVector3 **norms)
  {
  *norms = NULL;

  if (face_normals) {
    *norms = face_normals;
    return;
    }

  GrVector3 pverts[40], v1, v2, norm;
  float nx, ny, nz, mag;
  int i, j, k;

  // if a single polygon //

  if (!connectivity) {
    face_normals = new GrVector3[1];
    gc_PolyNormComp (num_vertices, vertices, face_normals[0]);
    *norms = face_normals;
    return;
    }

  // indexed polygon //

  int n = hgeom;
  int *conn = connectivity;
  face_normals = new GrVector3[number];

  for (i = 0; i < number; i++) {
    if (!hgeom) {
      n = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      k = *conn;
      pverts[j][0] = vertices[k][0];
      pverts[j][1] = vertices[k][1];
      pverts[j][2] = vertices[k][2];
      conn++;
      }

    v1 = pverts[1]   - pverts[0];
    v2 = pverts[n-1] - pverts[0];
    nx = v1[1]*v2[2] - v1[2]*v2[1];
    ny = v1[2]*v2[0] - v1[0]*v2[2];
    nz = v1[0]*v2[1] - v1[1]*v2[0];
    mag = sqrt(nx*nx + ny*ny + nz*nz);

    if (mag != 0.0) {
      //if (rev) mag = -mag;
      nx = nx / mag;
      ny = ny / mag;
      nz = nz / mag;
      }

    face_normals[i][0] = nx;
    face_normals[i][1] = ny;
    face_normals[i][2] = nz;
    }

  *norms = face_normals;
  }

}


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
//* cyl:                 c y l i n d e r                       *
//*============================================================*

#include "pm/gr/cyl.h"
#include "gc.h"

namespace PmGraphics {

GrCylinder::~GrCylinder ()
  {
  fprintf (stderr, "\n>>>>>> GrCylinder::  dtor \n");
  }

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////


//*============================================================*
//*==========              setAxis                   ==========*
//*============================================================*
// set axes for cylinders.

void
GrCylinder::setAxis (GrVector3 axes[]) {

  if (!num_vertices) {
    return;
    }

  this->axes = new GrVector3[num_vertices];

  for (int i = 0; i < num_vertices; i++) {
    this->axes[i] = axes[i];
    this->axes[i].normalize();
    }
  }

void 
GrCylinder::setAxis(GrVector3& axis) 
  { 
  this->axis = axis; 
  this->axis.normalize(); 
  }

//*============================================================*
//*==========              setLength                 ==========*
//*============================================================*
// set lengths for cylinders.

void
GrCylinder::setLength (float lengths[]) {

  if (!num_vertices) {
    return;
    }

  this->lengths = new float[num_vertices];

  for (int i = 0; i < num_vertices; i++) {
    this->lengths[i] = lengths[i];
    }
  }

//*============================================================*
//*==========              pick                      ==========*
//*============================================================*
// pick a cylinder.

void
GrCylinder::pick (GrPickResult& pick)
  {
  //fprintf (stderr, "\n>>>>>> GrCylinder::pick: \n");

  if (!pick_enabled) {
    return;
    }

  GrVector3 line[2], dir, ipt, min_ipt;
  line[0] = pick.line[0];
  line[1] = pick.line[1];
  /*
  fprintf (stderr, "   >>> line[0] (%g %g %g) \n", line[0][0], line[0][1], line[0][2] );
  fprintf (stderr, "   >>> line[1] (%g %g %g) \n", line[1][0], line[1][1], line[1][2] );
  */

  bool ipt_found = false, intersect;
  float tol = pick.atts.tol;
  float min_dist = 1e6;
  float r, l, dist;
  this->picked_entity = -1;

  for (int i = 0; i < num_vertices; i++) {
    if (radii) {
      r = radii[i];
      }
    else {
      r = radius;
      }

    if (lengths) {
      l = lengths[i];
      }
    else {
      l = length;
      }

    if (axes) {
      dir = axes[i];
      }
    else {
      dir = axis;
      }

    gc_CylLineInt (r, l, vertices[i], dir, line, tol, intersect, ipt);

    if (intersect) {
      ipt_found = true;
      this->picked_entity = i+1;
      break;
      }
    }

  if (!ipt_found) {
    return;
    }

  GrPickGeom pg;
  pg.geom = this;
  pg.int_pt = ipt;
  pg.dist = 0.0;
  pg.intersected = true;
  pick.geom_list.push_back(pg);
  }

//*============================================================*
//*==========              setRadius                 ==========*
//*============================================================*
// set radii for cylinders.

void
GrCylinder::setRadius (float rads[]) {

  if (!num_vertices) {
    return;
    }

  this->radii = new float[num_vertices];

  for (int i = 0; i < num_vertices; i++) {
    this->radii[i] = rads[i]; 
    }
  }

//*============================================================*
//*==========              render                    ==========*
//*============================================================*
// render cylinders.

void
GrCylinder::render() {

  if (!active) {
    return;
    }

  api_setup();

  api_render ();

  api_reset_state();
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              initData                  ==========*
//*============================================================*
// initialize cylinder data.

void
GrCylinder::initData() 
  {
  fast_render = true;
  axis.set(1,0,0); 
  axes = NULL; 
  radius = 0.5;
  radii = NULL;
  length = 1.0;
  lengths = NULL;
  prv_data = NULL;
  }

//*============================================================*
//*==========              genCylinderGeom           ==========*
//*============================================================*
// generate cylinder geometry.  

void
GrCylinder::genCylinderGeom (float rad, float length, GrVector3 origin, GrVector3 axis,
                             int *p_num_polys, int **p_conn, int *p_num_verts, 
                             GrVector3 **p_verts, GrVector3 **p_vnorms,
                             GrVector3 **p_fnorms) 
  {

  int num_tri, *conn;

  int num_verts;

  GrVector3 *verts, *fnorms, *vnorms; 

  int i, j;

  float x, y, z;

  float cx, cy, cz; 

  GrVector3 u, v, w, p;

  float h, t, dt, ct, st, mag;

  int n, i1, i2, i3, j1, j2, j3;

  #ifdef dbg_genCylinderGeom
  fprintf (stderr, ">>>>>> GrCylinder::genCylinderGeom \n");
  fprintf (stderr, "   >>> radius [%g] \n", radius);
  fprintf (stderr, "   >>> length [%g] \n", length);
  #endif

  h = 0.5*length;
  n = 16;
  dt = 2.0*M_PI / n;
  t = 0.0;

  u = axis;
  gc_BasisComp (u, v, w);
  #ifdef dbg_genCylinderGeom
  fprintf (stderr, "   >>> u (%g %g %g) \n", u[0], u[1], u[2]);
  fprintf (stderr, "   >>> v (%g %g %g) \n", v[0], v[1], v[2]);
  fprintf (stderr, "   >>> w (%g %g %g) \n", w[0], w[1], w[2]);
  #endif

  num_verts = 2*n + 2;
  verts = new GrVector3[num_verts];
  vnorms = new GrVector3[num_verts];

  for (i = 0; i < n; i++) {
    ct = rad*cos(t);
    st = rad*sin(t);

    for (j = 0; j < 3; j++) {
      verts[i][j]   = origin[j] + ct*v[j] + st*w[j] - h*u[j];
      verts[i+n][j] = origin[j] + ct*v[j] + st*w[j] + h*u[j];
      vnorms[i][j]   = ct*v[j] + st*w[j];
      vnorms[i+n][j] = ct*v[j] + st*w[j];
      }
 
    t += dt;
    }

  for (j = 0; j < 3; j++) {
    verts[2*n][j]   = origin[j] - h*u[j];
    verts[2*n+1][j] = origin[j] + h*u[j];
    vnorms[2*n][j]  = -axis[j];
    vnorms[2*n+1][j] = axis[j];
    }

  num_tri = 4*n;
  conn = new int[3*num_tri];
  fnorms = new GrVector3[num_tri];
  num_tri = 0;

  for (i = 0; i < n; i++) {
    i1 = i; i2 = i+1; i3 = i+n;

    if (i2 == n) {
      i2 = 0;
      }

    conn[3*num_tri+0] = i1;
    conn[3*num_tri+1] = i2;
    conn[3*num_tri+2] = i3;
    fnorms[num_tri] = (1/3.0)*(vnorms[i1] + vnorms[i2] + vnorms[i3]);
    num_tri += 1;

    j1 = i2;
    j2 = i2+n;
    j3 = i3;

    conn[3*num_tri+0] = j1;
    conn[3*num_tri+1] = j2;
    conn[3*num_tri+2] = j3;
    fnorms[num_tri] = (1.0/3.)*(vnorms[j1] + vnorms[j2] + vnorms[j3]);
    num_tri += 1;
    }

  for (i = 0; i < n; i++) {
    i1 = i+1; i2 = i; i3 = 2*n;
    if (i1 == n) i1 = 0;
    conn[3*num_tri+0] = i1;
    conn[3*num_tri+1] = i2;
    conn[3*num_tri+2] = i3;
    fnorms[num_tri] = -axis;
    num_tri += 1;

    j2 = i+1+n; j1 = i+n; j3 = 2*n+1;
    if (j2 == 2*n) j2 = n;

    conn[3*num_tri+0] = j1;
    conn[3*num_tri+1] = j2;
    conn[3*num_tri+2] = j3;
    fnorms[num_tri] = axis;
    num_tri += 1;
    }

  *p_num_verts = num_verts;
  *p_verts = verts;
  *p_num_polys = num_tri;
  *p_conn = conn;
  *p_vnorms = vnorms;
  *p_fnorms = fnorms;
  }

}


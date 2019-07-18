
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
//* ellipsoid:              e l l i p s o i d                  *
//*============================================================*

#include "pm/gr/ellipsoid.h"
#include "gc.h"

namespace PmGraphics {


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

GrEllipsoid::~GrEllipsoid()
  {
  }

//*============================================================*
//*==========              pick                      ==========*
//*============================================================*
// pick a sphere.

void
GrEllipsoid::pick (GrPickResult& pick)
  {
  //fprintf (stderr, "\n>>>>>> GrEllipsoid::pick: \n");

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
  float r, dist;
  this->picked_entity = -1;

  for (int i = 0; i < num_vertices; i++) {

    compLineInt (line, vertices[i], r, intersect, dist, ipt);

    if (intersect) {
      fprintf (stderr, "   >>> intersect [%d]   dist [%g] \n", i, dist);
      ipt_found = true;
      dist = fabs (dist);

      if (dist < min_dist) {
        min_dist = dist;
        min_ipt = ipt;
        this->picked_entity = i+1;
        }
      }
    }

  if (!ipt_found) {
    return;
    }

  GrPickGeom pg;
  pg.geom = this;
  pg.int_pt = min_ipt;
  pg.dist = min_dist;
  pg.intersected = true;
  pick.geom_list.push_back(pg);
  }


//*============================================================*
//*==========              render                    ==========*
//*============================================================*
// render ellipsoids.

void
GrEllipsoid::render() {

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
//*==========              initAtts                  ==========*
//*============================================================*
// initialize line attributes.

void
GrEllipsoid::initAtts(GrEllipsoidData *data) 
  {
  fast_render = true;
  prv_data = NULL;
  axes_data = new GrEllipsoidData[num_vertices];

  for (int i = 0; i < num_vertices; i++) {
    axes_data[i].rad1 = data[i].rad1;
    axes_data[i].axis1 = data[i].axis1;
    axes_data[i].rad2 = data[i].rad2;
    axes_data[i].axis2 = data[i].axis2;
    axes_data[i].rad3 = data[i].rad3;
    axes_data[i].axis3 = data[i].axis3;
    }
  }

//*============================================================*
//*==========              genEllipsoidGeom          ==========*
//*============================================================*
// generate sphere geometry.  

void
GrEllipsoid::genEllipsoidGeom (GrEllipsoidData& data, GrVector3 center, 
                               int *p_num_polys, int **p_conn, int *p_num_verts, 
                               GrVector3 **p_verts, GrVector3 **p_fnorms)
  {

  GrGcTri tri[1000];

  int num_tri, num_verts, num_iverts, *conn; 

  GrVector3 *verts, *iverts, *fnorms, poly[3], *vnorms; 


  // get sphere geometry //

  gc_SphereGen (tri, 3, &num_tri);

  num_verts = 3 * num_tri;
  verts = new GrVector3[num_verts];
  fnorms = new GrVector3[num_tri];
  int n = 0;

  for (int i = 0; i < num_tri; i++) {
    verts[n++] = tri[i].p1;
    verts[n++] = tri[i].p2;
    verts[n++] = tri[i].p3;
    poly[0] = tri[i].p1;
    poly[1] = tri[i].p2;
    poly[2] = tri[i].p3;
    gc_PolyNormComp(3, poly, fnorms[i]);
    }

  gc_TriGetIndexes (num_tri, verts, &conn, &num_iverts, &iverts);

  // transform sphere geometry into ellipsoid //

  float a  = data.rad1;
  float b  = data.rad2;
  float c  = data.rad3;

  GrVector3 u = data.axis1;
  GrVector3 v = data.axis2;
  GrVector3 w = data.axis3;
  GrVector3 p;
  float cx, cy, cz, dpu, dpv, dpw;
  cx = center[0];
  cy = center[1];
  cz = center[2];

  for (int i = 0; i < num_iverts; i++) {
    p = iverts[i];
    dpu = p * u;
    dpv = p * v;
    dpw = p * w;
    iverts[i][0] = a*dpu*u[0] + b*dpv*v[0] + c*dpw*w[0] + cx;
    iverts[i][1] = a*dpu*u[1] + b*dpv*v[1] + c*dpw*w[1] + cy;
    iverts[i][2] = a*dpu*u[2] + b*dpv*v[2] + c*dpw*w[2] + cz;
    }

  *p_num_polys = num_tri;
  *p_conn = conn;
  *p_num_verts = num_iverts;
  *p_verts = iverts;
  *p_fnorms = fnorms;
  }

//*============================================================*
//*==========              compLineInt               ==========*
//*============================================================*
// intersect sphere with a line.

void 
GrEllipsoid::compLineInt (GrVector3 line[2], GrVector3 center, float radius, 
                          bool& intersect, float& dist, GrVector3 ipt)
  {

  float xc, yc, zc, x0, y0, z0, xd, yd, zd;

  float a, b, c, d; 

  float t, t0, t1; 

  GrVector3 v;

 /**************
  ***  body  ***
  **************/

  xc = center[0];  yc = center[1];  zc = center[2];
  x0 = line[0][0]; y0 = line[0][1]; z0 = line[0][2];
  v = line[0] - line[1];
  xd = v[0];  yd = v[1];  zd = v[2];
  v = line[0] - center;
  xc = v[0];  yc = v[1];  zc = v[2];

  a = xd*xd + yd*yd + zd*zd;
  b = 2.0 * (xd*xc  +  yd*yc  +  zd*zc);
  c = xc*xc + yc*yc + zc*zc - radius*radius;
  d = b*b - 4.0*a*c;

  if (d < 0.0) {
    intersect = false;
    return;                  
    }

  d = sqrt(d);
  t0 = (-b - d) / (2*a);

  if (t0 > 0.0) {
    t = t0;
    }
  else {
    t = (-b + d) / (2*a);
    }

  intersect = true;
  ipt[0] = x0 + t*xd;
  ipt[1] = y0 + t*yd;
  ipt[2] = z0 + t*zd;
  dist = t * sqrt(a);
  }


}


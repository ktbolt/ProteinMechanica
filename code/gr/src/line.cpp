
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
//* line:                 l i n e                              *
//*============================================================*

#include "pm/gr/line.h"

namespace PmGraphics {

typedef struct LineIterState {
  int m, n;
  int i, j, *conn;
  int num;
  int num_verts;
  GrVector3 *verts;
  int closed, disjoint;
  int inc;
  int entity;
  } LineIterState;

bool getNextLineSeg (bool& init, LineIterState& state, GrVector3& p1,
                     GrVector3& p2);

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pick                      ==========*
//*============================================================*
// pick a line geometry.

void
GrLine::pick (GrPickResult& pick)
  {
  #ifdef dbg_GrLine_pick 
  fprintf (stderr, "\n>>>>>> GrLine::pick:   \"%s\"\n", name.c_str());
  #endif

  if (!pick_enabled) {
    #ifdef dbg_GrLine_pick 
    fprintf (stderr, "   >>> pick not enabled \n");
    #endif
    return;
    }

  GrVector3 v, ldir, ipt, min_ipt, p1, p2, lp1, lp2;
  int n;
  bool ninit = true;
  LineIterState istate;

  istate.num = number;
  istate.num_verts = num_vertices;
  istate.verts = vertices;
  istate.conn = connectivity;
  istate.closed = closed;
  istate.disjoint = disjoint;

  #ifdef dbg_GrLine_pick 
  fprintf (stderr, "\n>>>>>> GrLine::pick:   \"%s\"\n", name.c_str());
  fprintf (stderr, "   >>> disjoint = %d \n", disjoint); 
  fprintf (stderr, "   >>> closed   = %d \n", closed); 
  fprintf (stderr, "   >>> num_vertices = %d \n", num_vertices); 
  fprintf (stderr, "   >>> number = %d \n", number); 
  fprintf (stderr, "   >>> size   = %d \n", size); 
  fprintf (stderr, "   >>> conn = %x \n", connectivity); 

  if (connectivity) {
    fprintf (stderr, "   >>> conn = "); 
    for (int i = 0; i < size; i++) {
      fprintf (stderr, " %d", istate.conn[i]); 
      }
    fprintf (stderr, "\n"); 
    }
  #endif

  // NOTE: xform.matrix must be set for this to work //

  if (xform.set) {
    #ifdef dbg_GrLine_pick 
    fprintf (stderr, "   >>> xform set \n");
    #endif
    PmMatrix3x3 mat = xform.matrix;
    mat.transpose();
    p1 = mat*(pick.line[0]-xform.center - xform.translation) + xform.center;
    p2 = mat*(pick.line[1]-xform.center - xform.translation) + xform.center;
    }
  else {
    p1 = pick.line[0];
    p2 = pick.line[1];
    }

  #ifdef dbg_GrLine_pick 
  fprintf (stderr, "   >>> number %d \n", number);
  #endif

  ldir = p2 - p1;
  bool ipt_found = false, intersect;
  float tol = pick.atts.tol;
  float min_dist = 1e6;
  float r, dist;
  this->picked_entity = -1;

  // extract the line segments for the object //

  while (getNextLineSeg(ninit, istate, lp1, lp2)) {
    v = lp2 - lp1; 
    compLineLineDist (lp1, v, p1, ldir, ipt, dist);
    #ifdef dbg_GrLine_pick 
    fprintf (stderr, "   >>> dist [%g] \n", dist);
    #endif

    if ((dist < min_dist) && (dist >= 0.0)) {
      min_dist = dist;
      min_ipt = ipt;
      ipt_found = true;
      //picked_entity = istate.i / istate.inc;
      picked_entity = istate.entity;

      if (dist == 0.0) {
        break;
        }
      }
    }

  if (!ipt_found || (min_dist > tol)) {
    return;
    }

  GrPickGeom pg;
  pg.geom = this;
  pg.int_pt = min_ipt;
  pg.dist = min_dist;
  pg.intersected = false;
  pick.geom_list.push_back(pg);
  }

//*============================================================*
//*==========              render                    ==========*
//*============================================================*
// render a line.

void
GrLine::render() {
  /*
  fprintf (stderr, "\n>>>>>> GrLine::render  name[%s] \n", name.c_str());
  fprintf (stderr, "\n>>>>>> GrLine::render  active [%d] \n", active);
  fprintf (stderr, "   >>> disjoint [%d] \n", disjoint);
  fprintf (stderr, "   >>> conn     [%x] \n", connectivity);
  */

  if (!active) {
    return;
    }

  api_setup();

  if (disjoint) {
    api_render_disjoint ();
    }

  else if (connectivity) {
    api_render_indexed ();
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
//*==========                 init                   ==========*
//*============================================================*
// initialize line attributes.

void
GrLine::init() 
  {
  //fprintf (stderr, "\n>>>>>> GrLine::init   this [%x] \n", this);
  //fprintf (stderr, "  >>>>>> GrLine::init   conn [%x] \n", connectivity);
  disjoint = false;
  closed = false;
  line_style = GR_LINE_STYLE_SOLID;
  line_width = 1.0;
  }

//*============================================================*
//*==========          compLineLineDist              ==========*
//*============================================================*
// compute the distance between two lines.

void
GrLine::compLineLineDist (GrVector3& p1, GrVector3& v1, GrVector3& p2, GrVector3& v2,
                          GrVector3& ipt, float& p_dist)
  {

  float u, t, dist; 
  float s1, s2;
  int i;

  float a =  v1 * v1;
  float b = -v1 * v2;
  float c = -b;
  float d = -v2 * v2;
  float e = v1*p2 - v1*p1;
  float f = v2*p2 - v2*p1;
  float det = a*d - c*b;

  if (det == 0.0) {
    p_dist = -1.0;
    return;
    }

  t = (e*d - f*b) / det;
  u = (a*f - e*c) / det;
  dist = 0.0;

  if ((t < 0.0) || (t > 1.0)) {
    p_dist = -1.0;
    return;
    }
  t = (e*d - f*b) / det;
  u = (a*f - e*c) / det;
  dist = 0.0;

  if ((t < 0.0) || (t > 1.0)) {
    p_dist = -1.0;
    return;
    }

  if ((u < 0.0) || (u > 1.0)) {
    p_dist = -1.0;
    return;
    }

  for (i = 0; i < 3; i++) {
    s1 = p1[i] + t*v1[i];
    s2 = p2[i] + u*v2[i];
    dist += (s1 - s2) * (s1 - s2);
    ipt[i] = s1;
    }

  p_dist = sqrt (dist);
  }

//*============================================================*
//*==========            getNextLineSeg              ==========*
//*============================================================*
// get a series of line segments.                

bool 
getNextLineSeg (bool& init, LineIterState& state, GrVector3& p1, GrVector3& p2)
  {
  #ifdef dbg_getNextLineSeg 
  fprintf (stderr, "------ getNextLineSeg ------ \n");
  #endif
  int i1, i2;

  if (init) {
    state.inc = 1;
    state.entity = 0;

    if (state.conn) {
      state.n = state.num;
      state.m = *state.conn++;
      state.i = 1;
      state.j = 1;
      //fprintf (stderr, ">>>>>> getNextLineSeg:  has conn \n");
      //fprintf (stderr, "                        state m [%d]  \n", state.m);
      }
    else {
      if (state.closed) {
        state.n = state.num_verts;
        }
      else {
        state.n = state.num_verts-1;
        }

      if (state.disjoint) {
        state.inc = 2;
        }

      state.i = 0;
      }

    init = false;
    }

  if (state.conn) {
    #ifdef dbg_getNextLineSeg 
    fprintf (stderr, "   >>> state i [%d]  state j [%d] \n", state.i, state.j);
    fprintf (stderr, "       state m [%d]  state n [%d] \n", state.m, state.n);
    #endif

    if (state.j == state.m) {
      if (state.i == state.n) {
        return false;
        }

      state.conn++;
      state.m = *state.conn++;
      state.j = 1;
      state.i += 1;
      #ifdef dbg_getNextLineSeg 
      fprintf (stderr, "   >>> state i [%d]  \n", state.i);
      fprintf (stderr, "       state m [%d]  \n", state.m);
      #endif

      if (state.m == 1) {
        return false;
        }
      }

    i1 = *state.conn++;
    i2 = *state.conn;
    #ifdef dbg_getNextLineSeg 
    fprintf (stderr, "   >>> i1 [%d]  i2 [%d] \n", i1, i2);
    #endif
    state.j += 1;
    }
  else if (state.i >= state.n) {
    return false;
    }
  else {
    i1 = state.i;
    i2 = (state.i+1) % state.num_verts;
    state.i += state.inc;
    }

  //fprintf (stderr, "i1 i2  %d %d \n", i1, i2);
  p1 = state.verts[i1];
  p2 = state.verts[i2];
  state.entity += 1;
  return true;
  }


}


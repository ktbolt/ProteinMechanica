
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
//* point:                 p o i n t                           *
//*============================================================*

#include "pm/gr/point.h"
#include "gc.h"

namespace PmGraphics {

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pick                      ==========*
//*============================================================*
// pick a point object.

void
GrPoint::pick (GrPickResult& pick)
  {
  //fprintf (stderr, "\n>>>>>> GrPoint::pick: \n");

  if (!pick_enabled) {
    return;
    }

  GrVector3 line[2], dir, ipt, min_ipt, proj_pt;
  line[0] = pick.line[0];
  line[1] = pick.line[1];
  //fprintf (stderr, "   >>> line[0] (%g %g %g) \n", line[0][0], line[0][1], line[0][2] );
  //fprintf (stderr, "   >>> line[1] (%g %g %g) \n", line[1][0], line[1][1], line[1][2] );

  bool ipt_found = false, intersect = false;
  float tol = pick.atts.tol;
  float min_dist = 1e6;
  float r, dist;
  this->picked_entity = -1;

  for (int i = 0; i < num_vertices; i++) {
    gc_LinePointProj (line, vertices[i], dist, proj_pt);
    //fprintf (stderr, "   >>> dist [%g] \n", dist);

    if (dist < min_dist) { 
      min_dist = dist;
      min_ipt = vertices[i];
      ipt_found = true;
      picked_entity = i+1;

      if (dist == 0.0) {
        intersect = true;
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
  pg.intersected = intersect;
  pick.geom_list.push_back(pg);

  }

//*============================================================*
//*==========              render                    ==========*
//*============================================================*
// render points.

void
GrPoint::render() {

  if (!active) {
    return;
    }

  api_setup();

  if (marker_type != GR_GEOMETRY_MARKER_NONE) {
    api_render_markers ();
    }
  else {
    api_render ();
    }

  api_reset_state();
  }

//*============================================================*
//*==========              setMarkerType             ==========*
//*============================================================*
// set the marker type.

void 
GrPoint::setMarkerType (const GrGeometryMarkerType type) {
  this->marker_type = type;
  }

//*============================================================*
//*==========              setMarkerSize             ==========*
//*============================================================*
// set the marker size.

void
GrPoint::setMarkerSize (const float size) {
  this->marker_size = size;
  }

//*============================================================*
//*==========                    setSize             ==========*
//*============================================================*
// set the point size.

void
GrPoint::setSize (const float size) {
  this->marker_size = size;
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              initAtts                  ==========*
//*============================================================*
// initialize point data.

void
GrPoint::initAtts () 
  {
  prv_data = NULL;
  marker_type = GR_GEOMETRY_MARKER_POINT;
  marker_size = 1.0;
  }


}



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
//* ellipsoid:             e l l i p s o i d                   *
//*============================================================*

#ifndef _ELLIPSOID_GR_H_
#define _ELLIPSOID_GR_H_

#include "pm/gr/gr.h"
#include "pm/gr/vgeom.h"

namespace PmGraphics {

typedef struct GrEllipsoidData {
  float rad1, rad2, rad3;
  GrVector3 axis1, axis2, axis3;
  } GrEllipsoidData;


// GrEllipsoid
// -----------

 class PM_EXPORT GrEllipsoid : public GrVertexGeometry {

   public:

      GrEllipsoid(const string name, int num_verts, GrVector3 *verts, 
                  GrEllipsoidData *data) :
         GrVertexGeometry(name, num_verts, verts) { initAtts(data); }

      ~GrEllipsoid();

      void render();

   private:
     GrEllipsoidData *axes_data;
     bool fast_render;
     void *prv_data;

     void genEllipsoidGeom (GrEllipsoidData& data, GrVector3 center, int *p_num_polys, 
                            int **p_conn, int *p_num_verts, GrVector3 **p_verts, 
                            GrVector3 **p_norms);
     void compLineInt (GrVector3 line[2], GrVector3 center, float radius, bool& intersect,
                       float& dist, GrVector3 ipt);

     void pick (GrPickResult& pick);

     void initAtts(GrEllipsoidData *data);
     void api_render();
     void api_reset_state();
     void api_setup();
   };

}

#endif




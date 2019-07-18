
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
//* sphere:                s p h e r e                         *
//*============================================================*

#ifndef _SPHERE_GR_H_
#define _SPHERE_GR_H_

#include "pm/gr/gr.h"
#include "pm/gr/vgeom.h"

namespace PmGraphics {


// GrSphere
// --------

 class PM_EXPORT GrSphere : public GrVertexGeometry {

   public:

      GrSphere (const string name, int num_verts, GrVector3 *verts) :
         GrVertexGeometry(name, num_verts, verts) { initAtts(); }

      ~GrSphere();

      void setOptimized(bool flag) { optimized = flag; }
      void setRadius (float rads[]);
      void setRadius (float radius);
      void setResolution(int res);
      void render();

   private:

     float *radii, radius;
     bool optimized;
     int resolution;
     void *prv_data;

     void genSphereGeom (float rad, GrVector3 center, int *p_num_polys, int **p_conn,
                         int *p_num_verts, GrVector3 **p_verts, GrVector3 **p_fnorms,
                         GrVector3 **p_vnorms);
     void compLineInt (GrVector3 line[2], GrVector3 center, float radius, bool& intersect,
                       float& dist, GrVector3& ipt);

     void pick (GrPickResult& pick);

     void initAtts();
     void api_render();
     void api_render_optimized();
     void api_reset_state();
     void api_setup();

     void api_gen_half_sphere(GrGeometryDisplayType dtype, float radius, int slices,
                              int stacks);

     void api_init_sphere_dlist();

   };

}

#endif




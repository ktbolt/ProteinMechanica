
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

#ifndef _POLYGON_GR_H_
#define _POLYGON_GR_H_

#include "pm/gr/gr.h"
#include "pm/gr/geom.h"

namespace PmGraphics {


// GrPolygon
// ---------

 class PM_EXPORT GrPolygon: public GrGeometry {

   public:

      GrPolygon (const string name, int num_verts, GrVector3 *verts) : 
         GrGeometry(name, num_verts, verts) { init(false); } 

      GrPolygon (const string name, GrGeometryType type, int num_verts, 
                 GrVector3 *verts) : 
         GrGeometry(name, num_verts, verts) { geom_type = type; init(false); } 

      GrPolygon(const string name, int num, GrIndex& conn, int num_verts, 
                GrVector3 *verts) :
         GrGeometry(name, num, conn, num_verts, verts) { init(true); } 

      GrPolygon(const string name, int num, GrIndex& conn, int num_verts, 
                GrVector3 *verts, int hgeom) :
         GrGeometry(name, num, conn, num_verts, verts, hgeom) { init(true); } 

      ~GrPolygon(); 

      void render();
      void setOutline(bool flag) { this->outline = flag; }

   private:

     GrVector3 *face_normals, *vertex_normals;
     bool outline;
     float line_width;
     GrLineStyleType line_style;

     void init(int has_conn);
     void getVertexNormals(GrVector3 **norms);
     void getFaceNormals(GrVector3 **norms);

     void pick (GrPickResult& pick);
     void pick_indexed (GrPickResult& pick);

     void api_render();
     void api_render_outline();
     void api_render_indexed();
     void api_setup();
     void api_reset_state();

   };


}

#endif




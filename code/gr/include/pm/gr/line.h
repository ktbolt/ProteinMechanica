
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

#ifndef _LINE_GR_H_
#define _LINE_GR_H_

#include "pm/gr/gr.h"
#include "pm/gr/geom.h"

namespace PmGraphics {


// GrLine
// ---------

 class PM_EXPORT GrLine : public GrGeometry {

   public:

      GrLine (const string name, int num_verts, GrVector3 *verts) :
         GrGeometry(name, num_verts, verts) { connectivity = NULL;
                                              vertex_colors = NULL;
                                              init(); }

      GrLine (const string name, int num_verts, GrVector3 *verts, bool copy) :
         GrGeometry(name, num_verts, verts, copy) { connectivity = NULL;
                                                    vertex_colors = NULL;
                                                    init(); }


      GrLine (const string name, int dflag, int num_verts, GrVector3 *verts) :
         GrGeometry(name, num_verts, verts) { connectivity = NULL;
                                              vertex_colors = NULL;
                                              init(); 
                                              disjoint = dflag; }

      GrLine (const string name, int num, GrIndex& conn, int num_verts,GrVector3 *verts) :
         GrGeometry(name, num, conn, num_verts, verts) { init (); }

      bool isClosed() { return closed; }
      void setClosed (const bool flag) { closed = flag; }

      bool isDisjoint () { return disjoint; }
      void setDisjoint (const bool flag) { disjoint = flag; }

      void setWidth (const float width ) { line_width = width; }

      void render();

   private:

     // attributes
     float line_width;
     GrLineStyleType line_style;
     
     bool disjoint;
     bool closed;

     void pick (GrPickResult& pick);
     void pick_indexed (GrPickResult& pick);

     void compLineLineDist (GrVector3& p1, GrVector3& v1, GrVector3& p2, GrVector3& v2,
                            GrVector3& ipt, float& p_dist);

     void api_render();
     void api_render_disjoint();
     void api_render_indexed();
     void api_reset_state();
     void api_setup();

     void init(); 

   };

}

#endif





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
//* geom:                 g e o m e t r y                      *
//*============================================================*

#ifndef _GEOMETRY_GR_H_
#define _GEOMETRY_GR_H_

#include "pm/gr/gr.h"
#include "pm/gr/vgeom.h"

namespace PmGraphics {


// GrGeometry
// ------------

 class PM_EXPORT GrGeometry : public GrVertexGeometry {

   public:

     GrGeometry (const string name, int num_verts, GrVector3 *verts) :
         GrVertexGeometry(name, num_verts, verts) { };

     GrGeometry (const string name, int num_verts, GrVector3 *verts, bool copy) :
         GrVertexGeometry(name, num_verts, verts, copy) { };

     GrGeometry (const string name, int num, GrIndex& conn, int num_verts, 
                 GrVector3 *verts);

     GrGeometry (const string name, int num, GrIndex& conn, int num_verts, 
                 GrVector3 *verts, int hgeom);

     ~GrGeometry(); 

     void setVertexColors (int num, const GrColor *colors);
     void setVertexColors (int num, const GrVector3 *colors);
     void setVertexColors (GrColorVector& colors);

   protected:

     int size, *connectivity;
     int hgeom;
     GrColor *vertex_colors;
     void setConnData(int num, GrIndex& conn);

   private:

     virtual void api_render() = 0;

   };

}


#endif




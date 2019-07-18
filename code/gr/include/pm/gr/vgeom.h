
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
//* vgeom:        v e r t e x      g e o m e t r y             *
//*============================================================*

#ifndef _VERTEX_GEOMETRY_GR_H_
#define _VERTEX_GEOMETRY_GR_H_

#include "pm/gr/gr.h"

namespace PmGraphics {


// GrVertexGeometry
// ----------------

 class PM_EXPORT GrVertexGeometry : private GrNoCopy {

   public:

     GrVertexGeometry ();
     GrVertexGeometry (const string name, int num_verts, GrVector3 *verts, 
                       bool copy=true);
     ~GrVertexGeometry();

     bool isActive();
     void setActive(const bool flag);

     void setColor (const GrColor& color) { this->color = color; }
     void setColor (const float red, const float green, const float blue) { 
                    this->color.set(red,green,blue); }
     void setColor (const GrVector3& color) { this->color.set(color); }
     void setLighting(bool flag) { this->twosided_lighting = flag; }

     void setColors (int num, GrColor *colors);
     void setColors (int num, GrVector3 *colors);
     void setColors (GrColorVector& colors);
     void setLineWidth(const float width);

     void setDisplayType (const GrGeometryDisplayType type) { this->display = type; }

     void getName (string& str) { str = name; }
     void getPickedEntity (int& id) { id = picked_entity; }

     void setScene(GrScene *scene) { this->scene = scene; }
     void setShadingType (const GrGeometryShadingType type) { this->shading = type; }
     void setMarkerSize (const float size) { this->marker_size = size; }

     void translate (const GrVector3& vec);
     void rotate (const GrVector3& vec);
     void scale (const GrVector3& vec);
     void scale (const float val);

     void setVertices (int num_verts, GrVector3 *verts);

     void applyXform();
     void setXform(GrXform& xform) { this->xform = xform; }
     bool hasXform() { return this->xform.set; }

     virtual void render() = 0;
     virtual void pick(GrPickResult& pick) = 0;

   protected:
     string name;
     bool active, twosided_lighting;
     GrXform xform;
     GrColor color;
     GrGeometryShadingType shading;
     GrGeometryDisplayType display;
     float line_width, marker_size;
     GrGeometryType geom_type;
     GrScene *scene;

     int number;
     int num_colors;
     GrColor *colors;

     bool copy_vertices;
     int num_vertices;
     GrVector3 *vertices;

     bool pick_enabled;
     int picked_entity;

     void setVertexData (int num_verts, GrVector3 *verts);
     void initAtts();
     void api_applyXform(GrXform& xform);

   private:

     virtual void api_render() = 0;

   };

}


#endif





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

#include "pm/gr/vgeom.h"

namespace PmGraphics {

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / detructor          ==========*
//*============================================================*

GrVertexGeometry::GrVertexGeometry()
  {
  this->number = 0;
  this->num_colors = 0;
  this->colors = NULL;
  this->copy_vertices = true;
  //fprintf (stderr, "\n>>>>>> GrVertexGeometry:: ctor  \n"); 
  }

GrVertexGeometry::GrVertexGeometry (const string name, int num_verts, GrVector3 *verts,
                                    bool copy) 
  {
  this->number = 0;
  this->num_colors = 0;
  this->colors = NULL;
  this->copy_vertices = copy;

  if (!num_verts || !verts) {
    return;
    }
  
  /*
  fprintf (stderr, "\n>>>>>> GrVertexGeometry:: ctor  name[%s]  copy = %d \n", 
           name.c_str(), copy);
  */

  this->name = name;
  setVertexData (num_verts, verts);
  initAtts();
  }

GrVertexGeometry::~GrVertexGeometry() {
  //fprintf (stderr, "\n>>>>>> GrVertexGeometry:: dtor \n");
  }

//*============================================================*
//*==========              isActive                  ==========*
//*============================================================*
// set if the geometry is active (visible).

bool 
GrVertexGeometry::isActive () { 
  return this->active; 
  }

void 
GrVertexGeometry::setActive (const bool flag) { 
  this->active = flag; 
  }

//*============================================================*
//*==========              setColors                 ==========*
//*============================================================*
// set the colors for each geometry entitiy. 

void 
GrVertexGeometry::setColors (int num, GrColor *colors) 
  {
  //fprintf (stderr, "\n>>>>>> GrVertexGeometry:: set colors \n");
  this->num_colors = num;
  this->colors = new GrColor[num]; 

  for (int i = 0; i < num; i++) {
    this->colors[i] = colors[i]; 
    }
  }

//*============================================================*
//*==========              setColors                 ==========*
//*============================================================*
// use color vector.

void
GrVertexGeometry::setColors (GrColorVector& colors)
  {
  //fprintf (stderr, "\n>>>>>> GrVertexGeometry:: set colors \n");
  int n = colors.size();
  this->num_colors = n;
  this->colors = new GrColor[n];

  for (int i = 0; i < n; i++) {
    //fprintf (stderr, "%g %g %g \n", colors[i].red, colors[i].green, colors[i].blue);
    this->colors[i] = colors[i];
    }
  }

//*============================================================*
//*==========              setColors                 ==========*
//*============================================================*
// use vectors.

void
GrVertexGeometry::setColors (int num, GrVector3 colors[])
  {

  this->num_colors = num;
  this->colors = new GrColor[num];

  for (int i = 0; i < num; i++) {
    //fprintf (stderr, "%g %g %g \n", colors[i][0], colors[i][1], colors[i][2]);
    this->colors[i].set(colors[i]);
    }
  }

//*============================================================*
//*==========               setLineWidth             ==========*
//*============================================================*
// set line width.

void 
GrVertexGeometry::setLineWidth(const float width) { 
  this->line_width = width; 
  }

//*============================================================*
//*==========               rotate                   ==========*
//*============================================================*
// set geometry's local rotation.

void
GrVertexGeometry::rotate (const GrVector3& vec) {
  this->xform.setRotation(vec[0], vec[1], vec[2]);
  this->xform.set = true;
  }

//*============================================================*
//*==========               scale                    ==========*
//*============================================================*
// set geometry's local scaling.  

void
GrVertexGeometry::scale (const GrVector3& vec) {
  this->xform.scale = vec;
  this->xform.set = true;
  }

void
GrVertexGeometry::scale (const float val) {
  this->xform.scale[0] = val;
  this->xform.scale[1] = val;
  this->xform.scale[2] = val;
  this->xform.set = true;
  }

//*============================================================*
//*==========            translate                   ==========*
//*============================================================*
// set geometry's local translation.

void 
GrVertexGeometry::translate (const GrVector3& vec) { 
  this->xform.translation = vec;
  this->xform.set = true;
  }

//*============================================================*
//*==========              setVertices               ==========*
//*============================================================*
// set the vertex data for a previously defined object.

void
GrVertexGeometry::setVertices (int num_verts, GrVector3 *verts)
  {
  /*
  fprintf (stderr, "\n>>>>>> GrVertexGeometry::setVertices\n");
  fprintf (stderr, "   >>> num verts[%d]\n", num_verts);
  fprintf (stderr, "   >>> num_vertices[%d]\n", num_vertices);
  */

  // Note [2Jan2010] this is being used by trace. //

  if (!copy_vertices) {
    //fprintf (stderr, "\n>>>>>> GrVertexGeometry::setVertices don't copy \n");
    //fprintf (stderr, "   >>> num verts[%d]  verts[%x] \n", num_verts, verts);
    this->num_vertices = num_verts;
    this->vertices = verts;

    /*
    for (int i = 0; i < num_verts; i++) {
      fprintf (stderr, "%d  %f %f %f \n", i, verts[i][0], verts[i][1], verts[i][2]);
      }
    */
    return;
    }

  if (num_verts != num_vertices) {
    return;
    }

  for (int i = 0; i < num_verts; i++) {
    this->vertices[i] = verts[i];
    }
  }

//*============================================================*
//*==========                applyXform              ==========*
//*============================================================*
// apply a geometry's local transformation.

void 
GrVertexGeometry::applyXform()
  {
  api_applyXform (xform); 
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              setVertexData             ==========*
//*============================================================*
// set vertex data.

void
GrVertexGeometry::setVertexData (int num_verts, GrVector3 *verts)
  {

  /*
  fprintf (stderr, "\n>>>>>> GrVertexGeometry:: setVertexData \n");
  fprintf (stderr, "   >>> num verts [%d] \n", num_verts);
  */

  if (!num_verts) {
    return;
    }

  this->num_vertices = num_verts;
  this->vertices = new GrVector3 [num_verts];

  for (int i = 0; i < num_verts; i++) {
    //fprintf (stderr, "%d  (%g %g %g) \n", i, verts[i][0], verts[i][1], verts[i][2]);
    this->vertices[i] = verts[i];
    }
  }

//*============================================================*
//*==========              initAtts                  ==========*
//*============================================================*
// initialze attributes.

void
GrVertexGeometry::initAtts() 
  {
  active = true;
  color.set(1, 1, 1);
  display = GR_GEOMETRY_DISPLAY_SOLID;
  shading = GR_GEOMETRY_SHADING_NONE;
  geom_type = GR_GEOMETRY_UNKNOWN;
  line_width = 1.0;
  pick_enabled = true;
  picked_entity = -1;
  //fprintf (stderr, "\n>>>>>> GrVertexGeometry::initAtts \n");
  twosided_lighting = false;
  }

}



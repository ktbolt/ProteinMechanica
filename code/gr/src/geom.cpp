
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

#include "pm/gr/geom.h"

namespace PmGraphics {

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / detructor          ==========*
//*============================================================*

GrGeometry::GrGeometry (const string name, int num, GrIndex& conn, int num_verts, 
                        GrVector3 *verts) 
  {
  this->number = 0;
  this->size = 0;
  this->connectivity = NULL;
  this->num_colors = 0;
  this->colors = NULL;
  this->vertex_colors = NULL;
  this->hgeom = 0;
  this->active = true;
  this->pick_enabled = true;

  if (!num_verts) {
    return;
    }

  //fprintf (stderr, "\n>>>>>> GrGeometry: conn ctor \n");
  this->name = name;
  setVertexData (num_verts, verts);
  setConnData (num, conn);
  }

GrGeometry::GrGeometry (const string name, int num, GrIndex& conn, int num_verts,
                        GrVector3 *verts, int hgeom) {
  this->number = 0;
  this->size = 0;
  this->connectivity = NULL;
  this->num_colors = 0;
  this->colors = NULL;
  this->vertex_colors = NULL;
  this->hgeom = hgeom;
  this->active = true;
  this->pick_enabled = true;

  if (!num_verts) {
    return;
    }

  //fprintf (stderr, "\n>>>>>> GrGeometry: conn ctor \n");

  this->name = name;
  setVertexData (num_verts, verts);
  setConnData (num, conn);
  }

GrGeometry::~GrGeometry ()
  {
  //fprintf (stderr, "\n>>>>>> GrGeometry: dtor  name [%s] \n", name.c_str());
  }

//*============================================================*
//*==========              setVertexColors           ==========*
//*============================================================*
// set the colors for each vertex of a geometry entitiy.

void 
GrGeometry::setVertexColors (int num, const GrColor *colors) 
  {

  if (num != num_vertices) { 
    return;
    }

  this->vertex_colors = new GrColor[num]; 

  for (int i = 0; i < num; i++) {
    this->vertex_colors[i] = colors[i]; 
    }
  }

//*============================================================*
//*==========              setVertexColors           ==========*
//*============================================================*
// set the colors for each vertex of a geometry entitiy.

void
GrGeometry::setVertexColors (int num, const GrVector3 *colors)
  {

  if (num != num_vertices) {
    return;
    }

  this->vertex_colors = new GrColor[num];

  for (int i = 0; i < num; i++) {
    this->vertex_colors[i] = colors[i];
    }
  }

//*============================================================*
//*==========              setVertexColors           ==========*
//*============================================================*
// set the colors for each vertex of a geometry entitiy using a
// color vector.

void
GrGeometry::setVertexColors (GrColorVector& colors)
  {

  if (!num_vertices) {
    return;
    }

 int n = colors.size();

 if ((n == 0) || (n != num_vertices)) {
   fprintf (stderr,
       "\n**** Error [GrGeometry::setVertexColors] color data wrong size.\n");
   return;
   }

  this->vertex_colors = new GrColor[num_vertices];

  for (int i = 0; i < num_vertices; i++) {
    this->vertex_colors[i] = colors[i];
    }
  }


////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              setConnData               ==========*
//*============================================================*
// set connectivity data.

void
GrGeometry::setConnData (int num, GrIndex& conn)
  {
  this->number = num;
  this->size = conn.size; 

  if (hgeom) {
    for (int i = 0; i < conn.size; i++) {
      int n = conn.vals[i];

      if (n >= this->num_vertices) {
        fprintf (stderr, 
          "\n**** Error [GrGeometry::setConnData] conn index [%d] out of range.\n", n); 
        return;
        }
      }
    }

  else { 
    int i = 0;

    while (i < conn.size) {
      int n = conn.vals[i++];
      //fprintf (stderr, "%d: ",  n);

      for (int j = 0; j < n; j++, i++) {
        if (conn[i] >= this->num_vertices) {
          fprintf (stderr, 
            "\n**** Error [GrGeometry::setConnData] conn index [%d] out of range.\n", i); 
          fprintf (stderr, "    num vertices[%d]  conn[%d] = %d\n", this->num_vertices,
                   i, conn[i]); 
          return;
          }

        //fprintf (stderr, "%d ", conn[i]);
        }
      //fprintf (stderr, "\n");
      }
    }

  this->connectivity = new int[conn.size];

  for (int i = 0; i < conn.size; i++) {
    this->connectivity[i] = conn.vals[i];
    }
  }

}



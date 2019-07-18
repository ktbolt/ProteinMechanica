
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
//* poly_opengl:        o p e n g l   p o l y g o n            *
//*============================================================*

#include "pm/gr/poly.h"
#include "gr_opengl.h"

namespace PmGraphics {

#define ndbg_GrPolygon

//*============================================================*
//*==========              api_setup                 ==========*
//*============================================================*

void
GrPolygon::api_setup()
  {
  glDisable (GL_LIGHTING);
  glLineWidth (line_width);
  glColor3f (color[0], color[1], color[2]);
  opengl_setLinePattern (line_style);
  glPointSize (marker_size);
  }

//*============================================================*
//*==========              api_reset_state           ==========*
//*============================================================*

void
GrPolygon::api_reset_state()
  {
  glEnable (GL_LIGHTING);
  glLineWidth (1.0);
  glDisable (GL_LINE_STIPPLE);
  }

//*============================================================*
//*==========              api_render                ==========*
//*============================================================*
// render a polygon. 

void
GrPolygon::api_render()
  {

  glBegin (GL_POLYGON);

  for (int i = 0; i < num_vertices; i++) {
    glVertex3f (vertices[i][0], vertices[i][1], vertices[i][2]);
    }

  glEnd ();
  }

//*============================================================*
//*==========              api_render_outline        ==========*
//*============================================================*
// render the outline of a polygon.

void
GrPolygon::api_render_outline()
  {

  glBegin (GL_LINE_STRIP);

  for (int i = 0; i < num_vertices; i++) {
    glVertex3f (vertices[i][0], vertices[i][1], vertices[i][2]);
    }

  glVertex3f (vertices[0][0], vertices[0][1], vertices[0][2]);
  glEnd ();
  }

//*============================================================*
//*==========          api_render_indexed            ==========*
//*============================================================*
// render a set of polygons indexed into its vertex array.

void
GrPolygon::api_render_indexed()
  {
  #define ndbg_GrPolygon
  #ifdef dbg_GrPolygon
  fprintf (stderr, "\n>>>>>>> GrPolygon::api_render_indexed: \n");
  fprintf (stderr, "    >>> display [%d] \n", display);
  fprintf (stderr, "    >>> hgeom   [%d] \n", hgeom);
  #endif
  int n, k;
  int *conn = connectivity;
  bool vcolors;
  GLenum mode;
  GrVector3 *norms = NULL;

  // get the display mode: point, line, fill  //

  opengl_getPolygonMode (display, mode);

  // get normals for shading //

  if (shading == GR_GEOMETRY_SHADING_COLOR) {
    getVertexNormals (&norms);
    glEnable (GL_LIGHTING);
    glShadeModel (GL_SMOOTH);
    }
  else if (shading == GR_GEOMETRY_SHADING_FLAT) { 
    getFaceNormals (&norms);
    glEnable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    }

  if (twosided_lighting) {
    glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    }

  if (vertex_colors) {
    vcolors = true;
    //fprintf (stderr, "    >>> has vertex colors \n");

    if (norms) {
      opengl_renderPolygonsColorsShading (mode, shading, number, conn, vertices, vcolors, 
                                          vertex_colors, norms, hgeom);
      }
    else {
      opengl_renderPolygonsColors (mode, number, conn, vertices, vcolors, vertex_colors, 
                                   hgeom);
      }
    }

  else if (colors) {
    vcolors = false;

    if (norms) {
      opengl_renderPolygonsColorsShading (mode, shading, number, conn, vertices, vcolors, 
                                          colors, norms, hgeom);
      }
    else {
      opengl_renderPolygonsColors (mode, number, conn, vertices, vcolors, colors, hgeom);
      }
    }

  else {

    if (norms) {
      opengl_renderPolygonsShading (mode, shading, number, conn, vertices, norms, hgeom);
      }
    else {
      opengl_renderPolygons (mode, number, connectivity, vertices, hgeom);
      }
    }

  if (twosided_lighting) {
    glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    }
  }


}   // PmGraphics


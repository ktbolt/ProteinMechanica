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
//* point_opengl:        o p e n g l   p o i n t               *
//*============================================================*

#include "pm/gr/point.h"
#include "gr_opengl.h"

namespace PmGraphics {

typedef struct PointPrvData {
  GrVector3 *verts;
  } PointPrvData;

//*============================================================*
//*==========              api_setup                 ==========*
//*============================================================*

void
GrPoint::api_setup()
  {

  if (!prv_data) {
    PointPrvData *prvd = new PointPrvData;
    prv_data = prvd;
    }

  glDisable (GL_LIGHTING);
  glPointSize (marker_size);
  //glPointSize (1.0);
  //printf ("##### %f \n", marker_size);
  glColor3f (color[0], color[1], color[2]);
  }

//*============================================================*
//*==========              api_reset_state           ==========*
//*============================================================*

void
GrPoint::api_reset_state()
  {
  glEnable (GL_LIGHTING);
  glLineWidth (1.0);
  glDisable (GL_LINE_STIPPLE);
  }

//*============================================================*
//*==========              api_render                ==========*
//*============================================================*
// render points.   

void
GrPoint::api_render()
  {

  PointPrvData *prvd = (PointPrvData*)prv_data;

  if (colors) {
    glBegin (GL_POINTS);

    for (int i = 0; i < num_vertices; i++) {
      glColor3fv ((GLfloat*)(colors+i));
      glVertex3fv ((GLfloat*)(vertices+i));
      }

    glEnd ();
    }
  else {
    glBegin (GL_POINTS);

    for (int i = 0; i < num_vertices; i++) {
      glVertex3fv ((GLfloat *)(vertices+i));
      }

    glEnd ();
    }
  }

//*============================================================*
//*==========              api_render_markers        ==========*
//*============================================================*
// render markers.

void
GrPoint::api_render_markers()
  {

  PointPrvData *prvd = (PointPrvData*)prv_data;
  GrColor *mcolors = colors;
  float mdim, x, y, z;
  mdim = marker_size;
  glLineWidth (2.0);

  switch (marker_type) {
    case GR_GEOMETRY_MARKER_CROSS:
      glBegin (GL_LINES);

      for (int i = 0; i < num_vertices; i++) {
        if (mcolors) {
          glColor3fv ((GLfloat*)(colors+i));
          mcolors++;
          }

        //if (scales) mdim = scales[i];
        x = vertices[i][0];
        y = vertices[i][1];
        z = vertices[i][2];
        glVertex3f (x-mdim, y, z);
        glVertex3f (x+mdim, y, z);

        glVertex3f (x, y-mdim, z);
        glVertex3f (x, y+mdim, z);

        glVertex3f (x, y, z-mdim);
        glVertex3f (x, y, z+mdim);
        }

      glEnd();
    break;

    default:
      api_render();
    }
  }



}  // PmGraphics


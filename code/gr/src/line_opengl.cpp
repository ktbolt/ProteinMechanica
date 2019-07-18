
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
//* line_opengl:        o p e n g l   l i n e                  *
//*============================================================*

#include "pm/gr/line.h"
#include "gr_opengl.h"

namespace PmGraphics {

//*============================================================*
//*==========              api_setup                 ==========*
//*============================================================*

void
GrLine::api_setup()
  {
  glDisable (GL_LIGHTING);
  glLineWidth (line_width);
  glColor3f (color[0], color[1], color[2]);
  opengl_setLinePattern (line_style);
  }

//*============================================================*
//*==========              api_reset_state           ==========*
//*============================================================*

void
GrLine::api_reset_state()
  {
  glEnable (GL_LIGHTING);
  glLineWidth (1.0);
  glDisable (GL_LINE_STIPPLE);
  }

//*============================================================*
//*==========              api_render                ==========*
//*============================================================*
// render a single connected line represented by a set of 
// vertices. 

void
GrLine::api_render()
  {
  /*
  fprintf (stderr, "\n>>>>>> GrLine::api_render \n");
  fprintf (stderr, "   >>> colors [%x] \n", colors);
  fprintf (stderr, "   >>> vertex_colors [%x] \n", vertex_colors);
  fprintf (stderr, "   >>> color (%f %f %f) \n", color[0], color[1], color[2]);
  */

  if (vertex_colors) { 
    glShadeModel (GL_SMOOTH);
    }

  if (display == GR_GEOMETRY_DISPLAY_POINT) {
    glBegin (GL_POINTS);
    }
  else { 
    if (closed) {
      glBegin (GL_LINE_LOOP);
      }
    else {
      glBegin (GL_LINE_STRIP);
      }
    }

  if ((colors) && (num_colors >= num_vertices-1)) {
    for (int i = 0; i < num_vertices-1; i++) {
      glColor3fv ((GLfloat*)(colors+i));
      glVertex3f (vertices[i][0], vertices[i][1], vertices[i][2]);
      glVertex3f (vertices[i+1][0], vertices[i+1][1], vertices[i+1][2]);
      }
    }

  else if (vertex_colors) { 
    //fprintf (stderr, "\n>>>>>> GrLine::api_render: vertex colors \n");

    for (int i = 0; i < num_vertices; i++) {
      glColor3fv ((GLfloat*)(vertex_colors+i));
      glVertex3f (vertices[i][0], vertices[i][1], vertices[i][2]);
      }
    }

  else {
    for (int i = 0; i < num_vertices; i++) {
      glVertex3f (vertices[i][0], vertices[i][1], vertices[i][2]);
      }
    }

  glEnd ();
  }

//*============================================================*
//*==========              api_render_disjoint       ==========*
//*============================================================*
// render disjoint lines.

void
GrLine::api_render_disjoint()
  {
  //fprintf (stderr, "\n>>>>>> GrLine::api_render_disjoint: \n");

  if ((colors) && (num_colors == num_vertices/2)) {
    int j;
    glBegin (GL_LINES);

    for (int i = 0, j = 0; i < num_vertices; i += 2, j++) {
      glColor3fv  ((GLfloat*)(colors+j));
      //fprintf (stderr, "%g %g %g\n", colors[j].red, colors[j].green, colors[j].blue);
      glVertex3fv ((GLfloat*)(vertices+i));
      glVertex3fv ((GLfloat*)(vertices+i+1));
      }

    glEnd ();
    }

  else if (vertex_colors) { 
    fprintf (stderr, "\n>>>>>> GrLine::api_render_disjoint: vertex colors \n");
    glShadeModel (GL_SMOOTH);
    glBegin (GL_LINES);

    for (int i = 0; i < num_vertices; i += 2) {
      glColor3fv  ((GLfloat*)(vertex_colors+i));
      glVertex3fv ((GLfloat*)(vertices+i));
      glColor3fv  ((GLfloat*)(vertex_colors+i+1));
      glVertex3fv ((GLfloat*)(vertices+i+1));
      }
    }

  else {
    glBegin (GL_LINES);

    for (int i = 0; i < num_vertices; i += 2) {
      glVertex3fv ((GLfloat *)(vertices+i));
      glVertex3fv ((GLfloat *)(vertices+i+1));
      }

    glEnd ();
    }
  }

//*============================================================*
//*==========         api_render_indexed             ==========*
//*============================================================*
// render indexed lines.

void
GrLine::api_render_indexed()
  {

  /*
  fprintf (stderr, "\n>>>>>> GrLine::api_render_indexed \n");
  fprintf (stderr, "   >>> colors [%x] \n", colors);
  fprintf (stderr, "   >>> vertex_colors [%x] \n", vertex_colors);
  fprintf (stderr, "   >>> color (%f %f %f) \n", color[0], color[1], color[2]);
  */

  int i, j, k, n;
  int *conn = connectivity;

  if (colors && (num_colors == number)) {
    int j;

    for (i = 0; i < number; i++) {
      glBegin (GL_LINE_STRIP);
      n = *conn;
      conn++;
      glColor3fv  ((GLfloat*)(colors+i));

      for (j = 0; j < n; j++) {
        k = *conn;
        glVertex3fv ((GLfloat*)(vertices+k));
        conn++;
        }

      glEnd ();
      }
    }

  else if (vertex_colors) { 
    glShadeModel (GL_SMOOTH);

    for (i = 0; i < number; i++) {
      glBegin (GL_LINE_STRIP);
      n = *conn;
      conn++;

      for (j = 0; j < n; j++) {
        k = *conn;
        glColor3fv  ((GLfloat*)(vertex_colors+k));
        glVertex3fv ((GLfloat*)(vertices+k));
        conn++;
        }

      glEnd ();
      }
    }

  else {
    for (i = 0; i < number; i++) {
      glBegin (GL_LINE_STRIP);
      n = *conn;
      //fprintf (stderr, "%d: ", n);
      conn++;

      for (j = 0; j < n; j++) {
        k = *conn;
        //fprintf (stderr, "%d ", k);
        glVertex3fv ((GLfloat*)(vertices+k));
        conn++;
        }

      //fprintf (stderr, "\n");
      glEnd ();
      }
    }
  }


}  // PmGraphics


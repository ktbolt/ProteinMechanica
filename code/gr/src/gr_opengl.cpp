
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
//* gr_opengl:        o p e n g l   u t i l i t i e s          *
//*============================================================*

#include "pm/gr/gr.h"
#include "gr_opengl.h"

namespace PmGraphics {

#define ndbg_polygons
#define ndbg_xform    


//*============================================================*
//*==========         opengl_setLinePattern          ==========*
//*============================================================*
// set line pattern. 

void
opengl_setLinePattern (GrLineStyleType type)
  {

  if (type != GR_LINE_STYLE_SOLID) {
    glEnable (GL_LINE_STIPPLE);
    glLineStipple (1, 0x0101);
    }
  }

//*============================================================*
//*==========            opengl_getPolygonMode       ==========*
//*============================================================*
// get the rendering mode for polygons.

void
opengl_getPolygonMode (GrGeometryDisplayType type, GLenum& mode)
  {
  switch (type) {
    case GR_GEOMETRY_DISPLAY_POINT:
      mode = GL_POINTS;
    break;

    case GR_GEOMETRY_DISPLAY_LINE:
      mode = GL_LINE_LOOP;
    break;

    case GR_GEOMETRY_DISPLAY_SOLID:
      mode = GL_POLYGON;
    break;
    }
  }

//*============================================================*
//*==========         opengl_renderPolygons          ==========*
//*============================================================*
// render polygons.  

void
opengl_renderPolygons (GLenum mode, int num, int *conn, GrVector3 *verts, int hgeom)
  {

  int i, j, k, n;
  n = hgeom;

  for (i = 0; i < num; i++) {
    glBegin (mode);
    if (!hgeom) n = *conn++;

    for (j = 0; j < n; j++) {
      k = *conn;
      glVertex3fv ((GLfloat*)(verts+k));
      conn++;
      }

    glEnd ();
    }
  }

//*============================================================*
//*==========         opengl_renderPolygonsColors    ==========*
//*============================================================*
// render polygons with colors defined at vertices or for each 
// polygon.

void
opengl_renderPolygonsColors (GLenum mode, int num, int *conn, GrVector3 *verts, 
                             bool vcolors, GrColor *colors, int hgeom)
  {

  //fprintf (stderr, "\n>>>>>> opengl_renderPolygonsColors \n");

  int i, j, k, n = hgeom;

  // colors for each vertex //

  if (vcolors) {
    glShadeModel (GL_SMOOTH);

    for (i = 0; i < num; i++) {
      glBegin (mode);
      if (!hgeom) n = *conn++;

      for (j = 0; j < n; j++) {
        k = *conn;
        glColor3fv ((GLfloat*)(colors+k));
        glVertex3fv ((GLfloat*)(verts+k));
        conn++;
        }

      glEnd ();
      }
    }


  // a single color for each face. //

  else {

    for (i = 0; i < num; i++) {
      glBegin (mode);
      if (!hgeom) n = *conn++;
      glColor3fv ((GLfloat*)(colors+i));

      for (j = 0; j < n; j++) {
        k = *conn;
        glVertex3fv ((GLfloat*)(verts+k));
        conn++;
        }

      glEnd ();
      }
    }
  }

//*============================================================*
//*==========    opengl_renderPolygonsColorsShading  ==========*
//*============================================================*
// render polygons with color and normal data.

void
opengl_renderPolygonsColorsShading (GLenum mode, GrGeometryShadingType shading, 
                                    int num, int *conn, GrVector3 *verts, bool vcolors, 
                                    GrColor *colors, GrVector3 *norms, int hgeom)
  {
  /*
  fprintf (stderr, "\n>>>>>> opengl_renderPolygonsColorsShading \n");
  fprintf (stderr, ">>> hgeom = %d \n", hgeom);
  fprintf (stderr, ">>> vcolors = %d \n", vcolors);
  */

  int i, j, k, n = hgeom;

  // single normal per polygon

  if (shading == GR_GEOMETRY_SHADING_FLAT) {
    if (vcolors) {
      glShadeModel (GL_SMOOTH);

      for (i = 0; i < num; i++) {
        glBegin (mode);
        if (!hgeom) n = *conn++;
        glNormal3fv ((GLfloat*)(norms));

        for (j = 0; j < n; j++) {
          k = *conn;
          glColor3fv ((GLfloat*)(colors+k));
          glVertex3fv ((GLfloat*)(verts+k));
          //fprintf(stderr,"%g %g %g\n",colors[k].red,colors[k].green,colors[k].blue);
          conn++;
          }

        glEnd ();
        norms++;
        }
      }

    else {
      for (i = 0; i < num; i++) {
        glBegin (mode);
        if (!hgeom) n = *conn++;
        glNormal3fv ((GLfloat*)(norms));
        glColor3fv ((GLfloat*)(colors));

        for (j = 0; j < n; j++) {
          k = *conn;
          glVertex3fv ((GLfloat*)(verts+k));
          conn++;
          }

        glEnd ();
        norms++;
        colors++;
        }
      }
    }

  // normals defined at the polygon vertices

  else if (shading == GR_GEOMETRY_SHADING_COLOR) {
    if (vcolors) {
      for (i = 0; i < num; i++) {
        glBegin (mode);
        if (!hgeom) n = *conn++;

        for (j = 0; j < n; j++) {
          k = *conn;
          glNormal3fv ((GLfloat *)(norms+k));
          glColor3fv ((GLfloat*)(colors+k));
          glVertex3fv ((GLfloat*)(verts+k));
          conn++;
          }

        glEnd ();
        }
      }

    else {
      for (i = 0; i < num; i++) {
        glBegin (mode);
        if (!hgeom) n = *conn++;

        for (j = 0; j < n; j++) {
          k = *conn;
          glNormal3fv ((GLfloat *)(norms+k));
          glVertex3fv ((GLfloat*)(verts+k));
          conn++;
          }

        glColor3fv ((GLfloat*)(colors+i));
        glEnd ();
        }
      }
    }
  }

//*============================================================*
//*==========      opengl_renderPolygonsShading      ==========*
//*============================================================*
// render polygons with shading.

void
opengl_renderPolygonsShading (GLenum mode, GrGeometryShadingType shading, int num, 
                              int *conn, GrVector3 *verts, GrVector3 *norms, int hgeom)
  {
  #ifdef dbg_polygons
  fprintf (stderr,">>>>>> opengl_renderPolygonsShading: \n");
  fprintf (stderr,"   >>> hgeom [%d] \n", hgeom);
  #endif 

  int i, j, k, n = hgeom;

  // single normal per polygon //

  if (shading == GR_GEOMETRY_SHADING_FLAT) {
    for (i = 0; i < num; i++) {
      glBegin (mode);
      if (!hgeom) n = *conn++;
      glNormal3fv ((GLfloat *)(norms));

      for (j = 0; j < n; j++) {
        k = *conn;
        glVertex3fv ((GLfloat*)(verts+k));
        conn++;
        }

      glEnd ();
      norms++;
      }
    }

  // normals defined at the polygon vertices //

  else if (shading == GR_GEOMETRY_SHADING_COLOR) {

    for (i = 0; i < num; i++) {
      glBegin (mode);
      if (!hgeom) n = *conn++;

      for (j = 0; j < n; j++) {
        k = *conn;
        glNormal3fv ((GLfloat*)(norms+k));
        glVertex3fv ((GLfloat*)(verts+k));
        conn++;
        }

      glEnd ();
      }
    }
  }

//*============================================================*
//*==========         opengl_applyXform              ==========*
//*============================================================*
// apply a transformation.             

void
opengl_applyXform (GrXform& xform)
  {

  #ifdef dbg_xform    
  fprintf (stderr, "  >>>> opengl_applyXform: is set [%d] \n", xform.set);
  #endif

  if (!xform.set) {
    return;
    }

  float tx, ty, tz, mat[16];
  GrMatrix3x3 xmat;

  tx = xform.center[0] + xform.translation[0];
  ty = xform.center[1] + xform.translation[1];
  tz = xform.center[2] + xform.translation[2];
  glTranslatef (tx, ty, tz);

  xmat = xform.matrix;
  opengl_getXformMatrix (xmat, mat);
  glMultMatrixf (mat);

  glScalef (xform.scale[0], xform.scale[0], xform.scale[0]);
  glTranslatef (-xform.center[0], -xform.center[1], -xform.center[2]);
  }

//*============================================================*
//*==========         opengl_applyInvXform           ==========*
//*============================================================*
// apply a inverse transformation.             

void
opengl_applyInvXform (GrXform& xform)
  {

  if (!xform.set) {
    return;
    }

  float tx, ty, tz, mat[16];
  PmMatrix3x3 xmat;

  tx = xform.center[0] + xform.translation[0];
  ty = xform.center[1] + xform.translation[1];
  tz = xform.center[2] + xform.translation[2];

  glTranslatef (xform.center[0], xform.center[1], xform.center[2]);
  glScalef (xform.scale[0], xform.scale[0], xform.scale[0]);

  xmat = xform.matrix;
  xmat.transpose();
  opengl_getXformMatrix (xmat, mat);
  glMultMatrixf(mat);

  glTranslatef (-tx, -ty, -tz);
  }

//*============================================================*
//*==========         opengl_GetXformMatrix          ==========*
//*============================================================*
// get the OpenGL transformation matrix.

void
opengl_getXformMatrix (PmMatrix3x3& matrix, float mat[16])
  {
  mat[0] = matrix(0,0);
  mat[1] = matrix(1,0);
  mat[2] = matrix(2,0);
  mat[3] = 0.0;

  mat[4] = matrix(0,1);
  mat[5] = matrix(1,1);
  mat[6] = matrix(2,1);
  mat[7] = 0.0;

  mat[8] = matrix(0,2);
  mat[9] = matrix(1,2);
  mat[10] = matrix(2,2);
  mat[11] = 0.0; 

  mat[12] = 0.0; 
  mat[13] = 0.0; 
  mat[14] = 0.0; 
  mat[15] = 1.0;
  }


}  // PmGraphics


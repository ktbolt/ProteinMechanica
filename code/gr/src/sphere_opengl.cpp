
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
//* sphere_opengl:        o p e n g l   s p h e r e            *
//*============================================================*

#include "pm/gr/scene.h"
#include "pm/gr/sphere.h"
#include "gr_opengl.h"

namespace PmGraphics {

typedef struct SphereGeom {
  int num_polys, *conn, num_verts; 
  GrVector3 *verts, *fnorms, *vnorms, *tmp_verts;
  } SphereGeom;

typedef struct SphereOptimized {
  int point_dlist;
  int line_dlist;
  int outline_dlist;
  int surf_dlist;
  SphereGeom geom;
  } SphereOptimized;

typedef struct SpherePrvData {
  SphereOptimized optimized;
  SphereGeom geom;
  GrVector3 *verts;
  } SpherePrvData;

//*============================================================*
//*==========              api_setup                 ==========*
//*============================================================*

void
GrSphere::api_setup()
  {

  // create private data. this will include a sphere of radius 
  // 1.0 centered at (0,0,0).

  if (!prv_data) {
    SpherePrvData *prvd = new SpherePrvData;
    prv_data = prvd;

    if (optimized) {
      api_init_sphere_dlist();
      }
    else {
      GrVector3 center(0,0,0);

      genSphereGeom (1.0, center, &prvd->geom.num_polys, &prvd->geom.conn,
                     &prvd->geom.num_verts, &prvd->geom.verts, &prvd->geom.fnorms,
                     &prvd->geom.vnorms);

      prvd->geom.tmp_verts = new GrVector3[prvd->geom.num_verts];
      }
    }

  glDisable (GL_LIGHTING);
  glLineWidth (line_width);
  glColor3f (color[0], color[1], color[2]);
  }

//*============================================================*
//*==========              api_reset_state           ==========*
//*============================================================*

void
GrSphere::api_reset_state()
  {
  glEnable (GL_LIGHTING);
  glLineWidth (1.0);
  glDisable (GL_LINE_STIPPLE);
  }

//*============================================================*
//*==========              api_render                ==========*
//*============================================================*
// render spheres.   

void
GrSphere::api_render()
  {

  //fprintf (stderr, "\n>>>>>>> GrSphere::api_render optimized = %d \n", optimized);

  if (optimized) {
    api_render_optimized();
    return;
    }

  //fprintf (stderr, "\n>>>>>>> GrSphere::api_render \n");

  GrVector3 v; 
  SpherePrvData *prvd = (SpherePrvData*)prv_data;

  int num_polys = prvd->geom.num_polys; 
  int *conn = prvd->geom.conn;
  int num_verts = prvd->geom.num_verts; 
  GrVector3 *verts = prvd->geom.verts; 
  GrVector3 *tmp_verts = prvd->geom.tmp_verts; 
  GrVector3 *fnorms = prvd->geom.fnorms;
  GrVector3 *vnorms = prvd->geom.vnorms;
  GrVector3 *norms;
  bool vcolors;
  GLenum mode;
  float r = radius;

  // get the display mode: point, line, fill.
  opengl_getPolygonMode (display, mode);

  // get normals for shading
  if (shading == GR_GEOMETRY_SHADING_COLOR) {
    glEnable (GL_LIGHTING);
    glShadeModel (GL_SMOOTH);
    norms = vnorms;
    }
  else if (shading == GR_GEOMETRY_SHADING_FLAT) {
    glEnable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    norms = fnorms;
    }
  else if (shading == GR_GEOMETRY_SHADING_NONE) { 
    norms = NULL;
    }

  //fprintf (stderr, "    >>> colors [%x] \n", colors); 
  //fprintf (stderr, "    >>> norms  [%x] \n", norms); 
  int hgeom = 0;

  for (int i = 0; i < num_vertices; i++) {
    v = vertices[i];
    //fprintf (stderr, "    >>> v (%g %g %g) \n", v[0], v[1], v[2]);

    if (radii) {
      r = radii[i];
      }

    for (int j = 0; j < num_verts; j++) {
      tmp_verts[j] = r*verts[j] + v;
      }

    if (colors) {
      glColor3fv ((GLfloat*)(colors+i));
      }

    if (norms) {
      opengl_renderPolygonsShading (mode, shading, num_polys, conn, tmp_verts, norms, 
                                    hgeom);
      }

    else {
      opengl_renderPolygons (mode, num_polys, conn, tmp_verts, hgeom);
      }
    }
  }

//*============================================================*
//*==========         api_render_optimized           ==========*
//*============================================================*
// render spheres using display lists.

void
GrSphere::api_render_optimized()
  {
  //fprintf (stderr, "\n>>>>>>> GrSphere::api_render_optimized \n");
  //fprintf (stderr, "   >>>> name %s \n", name.c_str());
  SpherePrvData *prvd = (SpherePrvData*)prv_data;
  GLenum mode;
  float r = radius;
  float cx, cy, cz; 
  GrVector3 v, angles;
  int dlist;
  float mat[16];
  GrMatrix3x3 xmat;

  // get scene xform //
  GrXform scene_xform;
  this->scene->getRotationAngles(angles);

  // get the display mode: point, line, fill.
  opengl_getPolygonMode (display, mode);

  // get normals for shading
  if (shading == GR_GEOMETRY_SHADING_COLOR) {
    glEnable (GL_LIGHTING);
    glShadeModel (GL_SMOOTH);
    }
  else if (shading == GR_GEOMETRY_SHADING_FLAT) {
    glEnable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    }
  else if (shading == GR_GEOMETRY_SHADING_NONE) {
    }

  if (display == GR_GEOMETRY_DISPLAY_POINT) {
    glDisable (GL_LIGHTING);
    glPointSize (1.0);
    dlist = prvd->optimized.point_dlist;
    }
  else if (display == GR_GEOMETRY_DISPLAY_LINE) {
    glDisable (GL_LIGHTING);
    dlist = prvd->optimized.line_dlist;
    }
  else if (display == GR_GEOMETRY_DISPLAY_OUTLINE) {
    glDisable (GL_LIGHTING);
    dlist = prvd->optimized.outline_dlist;
    }
  else {
    glEnable (GL_LIGHTING);
    dlist = prvd->optimized.surf_dlist;
    }

  // if a local xform is set then back transform //

  if (xform.set) {
    xmat = xform.matrix;
    xmat.transpose();
    opengl_getXformMatrix (xmat, mat);
    }

  for (int i = 0; i < num_vertices; i++) {
    v = vertices[i];

    if (radii) {
      r = radii[i];
      }

    if (colors) {
      glColor3fv ((GLfloat*)(colors+i));
      }

    glPushMatrix();

    // translate to sphere center //
    glTranslatef (v[0], v[1], v[2]);

    // scale by sphere radius //
    glScalef (r, r, r);

    // if a local xform is set then back transform //
    if (xform.set) {
      glMultMatrixf (mat);
      }

    // rotate to face eye //
    glRotatef (-angles[2], 0.0, 0.0, 1.0);
    glRotatef (-angles[1], 0.0, 1.0, 0.0);
    glRotatef (-angles[0], 1.0, 0.0, 0.0);

    glCallList (dlist);

    glPopMatrix();
    }
  }

//*============================================================*
//*==========         api_init_sphere_dlist          ==========*
//*============================================================*
// initialize display lists for sphere.

void
GrSphere::api_init_sphere_dlist()
  {
  //fprintf (stderr, "\n>>>>>>> GrSphere::api_init_sphere_dlist \n");
  SpherePrvData *prvd = (SpherePrvData*)prv_data;

  float r; 
  int dlist;
  int res = resolution;
  int num_stack = (res - 6 > 2 ? res : 2); 
  //fprintf (stderr, "   >>>> res[%d] \n", res);
  //fprintf (stderr, "   >>>> num stack[%d] \n", num_stack);

  // generate display lists for various display types // 

  // points //
  r = 1.0;
  dlist = glGenLists (1);
  glNewList (dlist, GL_COMPILE);
  api_gen_half_sphere (GR_GEOMETRY_DISPLAY_POINT, r, res, num_stack);
  glEndList();
  prvd->optimized.point_dlist = dlist; 

  // lines //
  r = 1.0;
  dlist = glGenLists (1);
  glNewList (dlist, GL_COMPILE);
  api_gen_half_sphere (GR_GEOMETRY_DISPLAY_LINE, r, res, num_stack);
  glEndList();
  prvd->optimized.line_dlist = dlist; 

  // outline //
  r = 1.0;
  dlist = glGenLists (1);
  glNewList (dlist, GL_COMPILE);
  api_gen_half_sphere (GR_GEOMETRY_DISPLAY_OUTLINE, r, res, num_stack);
  glEndList();
  prvd->optimized.outline_dlist = dlist;

  // surface //
  r = 1.0;
  dlist = glGenLists (1);
  glNewList (dlist, GL_COMPILE);
  api_gen_half_sphere (GR_GEOMETRY_DISPLAY_SOLID, r, res, num_stack);
  glEndList();
  prvd->optimized.surf_dlist = dlist;
  }

//*============================================================*
//*==========         api_gen_half_sphere            ==========*
//*============================================================*
// generate a half sphere for various display types.

void 
GrSphere::api_gen_half_sphere(GrGeometryDisplayType dtype, float radius, int slices, 
                              int stacks)
  {

  float rho, drho, theta, dtheta, x, y, z;
  int i, j, imin, imax;

  drho = M_PI / (float) stacks;
  dtheta = 2.0 * M_PI / (float) slices;

  imin = 1;
  imax = stacks / 2;

  if (dtype == GR_GEOMETRY_DISPLAY_POINT) {
    glBegin (GL_POINTS);

    for (j = 0; j <= slices; j++) {
      theta = (j == slices) ? 0.0 : j * dtheta;
      x = (float)(-sin(theta) * sin(drho));
      y = (float)(cos(theta) * sin(drho));
      z = (float)(cos(drho));
      glVertex3f (x*radius, y*radius, z*radius);
      }

    for (i = imin; i < imax; i++) {
      rho = i * drho;

      for (j = 0; j <= slices; j++) {
        theta = (j == slices) ? 0.0f : j * dtheta;
        x = (float)(-sin(theta) * sin(rho));
        y = (float)(cos(theta) * sin(rho));
        z = (float)(cos(rho));
        glVertex3f (x * radius, y * radius, z * radius);
        }
      }

    glEnd();
    return;
    }

  if (dtype == GR_GEOMETRY_DISPLAY_LINE) {
    for (i = imin; i < imax; i++) {
      rho = i * drho;
      glBegin (GL_LINE_LOOP);

      for (j = 0; j <= slices; j++) {
        theta = (j == slices) ? 0.0f : j * dtheta;
        x = (float)(-sin(theta) * sin(rho));
        y = (float)(cos(theta) * sin(rho));
        z = (float)(cos(rho));
        glVertex3f (x*radius, y*radius, z*radius);
        }

      glEnd();
      }

    glBegin (GL_LINES);

    for (j = 0; j <= slices; j++) {
      theta = (j == slices) ? 0.0 : j * dtheta;
      x = (float)(-sin(theta) * sin(rho));
      y = (float)(cos(theta) * sin(rho));
      z = (float)(cos(rho));
      glVertex3f (0.0, 0.0, 0.0);
      glVertex3f (x*radius, y*radius, z*radius);
      }

    glEnd();
    return;
    }

  if (dtype == GR_GEOMETRY_DISPLAY_OUTLINE) {
    glBegin (GL_LINE_LOOP);
    rho = (imax-1) * drho;

    for (j = 0; j <= slices; j++) {
      theta = (j == slices) ? 0.0f : j * dtheta;
      x = (float)(-sin(theta) * sin(rho));
      y = (float)(cos(theta) * sin(rho));
      z = (float)(cos(rho));
      glVertex3f (x*radius, y*radius, z*radius);
      }

    glEnd();

    return;
    }

  // draw top as a triangle fan // 

  glBegin (GL_TRIANGLE_FAN);
  glNormal3f (0.0, 0.0, 1.0);
  glVertex3f (0.0, 0.0, radius);

  for (j = 0; j <= slices; j++) {
    theta = (j == slices) ? 0.0 : j * dtheta;
    x = (float)(-sin(theta) * sin(drho));
    y = (float)(cos(theta) * sin(drho));
    z = (float)(cos(drho));

    glNormal3f (x, y, z);
    glVertex3f (x * radius, y * radius, z * radius);
    }

  glEnd();

  // draw intermediate stacks as quad strips // 

  for (i = imin; i < imax; i++) {
    rho = i * drho;

    glBegin (GL_QUAD_STRIP);

    for (j = 0; j <= slices; j++) {
      theta = (j == slices) ? 0.0f : j * dtheta;
      x = (float)(-sin(theta) * sin(rho));
      y = (float)(cos(theta) * sin(rho));
      z = (float)(cos(rho));
      glNormal3f(x, y, z);
      glVertex3f(x * radius, y * radius, z * radius);

      x = (GLfloat)(-sin(theta) * sin(rho + drho));
      y = (GLfloat)(cos(theta) * sin(rho + drho));
      z = (GLfloat)(cos(rho + drho));
      glNormal3f(x, y, z);
      glVertex3f(x * radius, y * radius, z * radius);
      }

    glEnd();
    }
  }



}

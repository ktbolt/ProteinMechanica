
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
// * render_opengl:       o p e n g l    i n e r f a c e       *
//*============================================================*

#include "pm/gr/gr.h"
#include "pm/gr/win.h"
#include "pm/gr/scene.h"
#include "scene_prv.h"
#include "gr_opengl.h"

namespace PmGraphics {

const int PICK_BUFFER_SIZE = 1000;

//*============================================================*
//*==========              api_createRenderer        ==========*
//*============================================================*
// create a renderer.

void
GrScene::api_createRenderer() 
  {
  //fprintf (stderr, "\n>>>>>> GrScene::create_renderer \n");
  ScenePrvData *scene_data = (ScenePrvData*)prv_data;

  static float ambient[] = {0.3, 0.3, 0.3, 0.0};

  glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
  glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, ambient);
  glEnable (GL_COLOR_MATERIAL);

  //glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glEnable (GL_NORMALIZE);
  glEnable (GL_LIGHTING);
  glShadeModel (GL_FLAT);

  glDepthFunc (GL_LESS);
  glEnable (GL_DEPTH_TEST);

  glMatrixMode (GL_MODELVIEW);
  }

//*============================================================*
//*==========              api_render                ==========*
//*============================================================*
// render the geometry defined for a scene.

void
GrScene::api_render()
  {
  //fprintf (stderr, "\n>>>>>> GrScene::api_render \n");
  ScenePrvData *prvd = (ScenePrvData*)prv_data;
  float tx, ty, tz, cx, cy, cz;
  GrColor color;

  if (!prvd->init_renderer) {
    api_createRenderer();
    prvd->init_renderer = true;
    }

  float x = -xform.center[0];

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  window->getColor(color);
  glClearColor (color[0], color[1], color[2], 1.0);

  glLoadIdentity ();
  glPushMatrix ();

  // set xform matrix for scene //

  tx = xform.center[0] + xform.translation[0];
  ty = xform.center[1] + xform.translation[1];
  tz = xform.center[2] + xform.translation[2];
  glTranslatef (tx, ty, tz);

  glRotatef (xform.angles[0], 1.0, 0.0, 0.0);
  glRotatef (xform.angles[1], 0.0, 1.0, 0.0);
  glRotatef (xform.angles[2], 0.0, 0.0, 1.0);

  glScalef (xform.scale[0], xform.scale[0], xform.scale[0]);

  cx = xform.center[0];
  cy = xform.center[1];
  cz = xform.center[2];
  glTranslatef (-cx, -cy, -cz);
  //fprintf (stderr, ">>> center = %f %f %f \n", cx, cy, cz);

  // traverse geometries for the scene //

  //fprintf (stderr, "  >>>> num geoms [%d] \n", this->geometry.size());

  for (int i = 0; i < this->geometry.size(); i++) {
    GrVertexGeometry *geom = geometry[i];
    string name;
    geom->getName(name);
    //fprintf (stderr, "  >>>> %s geom active [%d] \n", name.c_str(), geom->isActive()); 

    if (!geom->isActive()) {
      continue;
      }

    glPushMatrix ();

    if (geom->hasXform()) {
      geom->applyXform();
      }

    if (prvd->pick.perform) {
      int pick_id = (uintptr_t)geom;  //PROBLEM 
      //#int pick_id = 0; 
      glLoadName (pick_id);
      }

    geom->render();
    glPopMatrix ();
    }

  glPopMatrix ();

  if (!prvd->pick.perform) {
    this->window->swapBuffers();
    }
  }

//*============================================================*
//*==========              api_setViewport           ==========*
//*============================================================*
// set viewport.                                  

void
GrScene::api_setViewport()
  {
  //fprintf (stderr, "\n>>>>>> GrScene::api_setViewport \n");
  ScenePrvData *prvd = (ScenePrvData*)prv_data;

  int w = prvd->viewport.win_width;
  int h = prvd->viewport.win_height;
  GrExtent extent = prvd->viewport.extent;

  /*
  fprintf (stderr, "extent           \n"); 
  fprintf (stderr, "       x = %f %f \n", extent.min[0], extent.max[0]);
  fprintf (stderr, "       y = %f %f \n", extent.min[1], extent.max[1]);
  fprintf (stderr, "       z = %f %f \n", extent.min[2], extent.max[2]);
  */

  glViewport (0, 0, w, h);

  glMatrixMode (GL_PROJECTION);

  glLoadIdentity();

  GrVector3 min, max;

  extent.get (min, max);

  /*
  fprintf (stderr, "  >>>> viewport   w [%d]  h [%d] \n", w, h);
  fprintf (stderr, "  >>>> extent   min (%g %g %g) \n", min[0], min[1], min[2]);
  fprintf (stderr, "                max (%g %g %g) \n", max[0], max[1], max[2]);
  */

  if (prvd->pick.perform) {
    GLint viewport[4];
    int x, y;
    x = prvd->pick.screen_pt[0];
    y = prvd->pick.screen_pt[1];
    glGetIntegerv (GL_VIEWPORT, viewport);
    gluPickMatrix((GLdouble)x, (GLdouble)(y), 5.0, 5.0, viewport);
    }

  glOrtho (min[0], max[0], min[1], max[1], min[2], max[2]);
  glMatrixMode (GL_MODELVIEW);
  }

//*============================================================*
//*==========              api_unprojectPoint        ==========*
//*============================================================*
// transform screem coordinates to world coordinates.

void
GrScene::api_unprojectPoint (int screen_pt[2], GrVector3 world_pts[2]) 
  {

  GLint viewport[4];

  GLdouble mvmatrix[16], projmatrix[16];

  GLdouble x, y, z;

  GLdouble wx, wy, wz;

  float tx, ty, tz;

  // load transformation

  glLoadIdentity ();
  glPushMatrix ();

  tx = xform.center[0] + xform.translation[0];
  ty = xform.center[1] + xform.translation[1];
  tz = xform.center[2] + xform.translation[2];

  glTranslatef (tx, ty, tz);

  glRotatef (xform.angles[0], 1.0, 0.0, 0.0);
  glRotatef (xform.angles[1], 0.0, 1.0, 0.0);
  glRotatef (xform.angles[2], 0.0, 0.0, 1.0);

  glScalef (xform.scale[0], xform.scale[0], xform.scale[0]);

  glTranslatef (-xform.center[0], -xform.center[1], -xform.center[2]);

  //  get view matrices

  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  //  unproject screen points

  x = screen_pt[0];
  y = screen_pt[1];
  z = 0.0;

  gluUnProject (x, y, z, mvmatrix, projmatrix, viewport, &wx, &wy, &wz);
  world_pts[0][0] = wx;
  world_pts[0][1] = wy;
  world_pts[0][2] = wz;

  z = 1.0;
  gluUnProject (x, y, z, mvmatrix, projmatrix, viewport, &wx, &wy, &wz);
  world_pts[1][0] = wx;
  world_pts[1][1] = wy;
  world_pts[1][2] = wz;

  glPopMatrix ();
  }

//*============================================================*
//*==========              api_pick                  ==========*
//*============================================================*
// perform a pick.                                     

void
GrScene::api_pick (int screen_pt[2], GrPickResult& pick)
  {

  GrVector3 world_pts[2];

  GLuint select_buf[PICK_BUFFER_SIZE];

  api_unprojectPoint (screen_pt, world_pts);

  pick.line[0] = world_pts[0];
  pick.line[1] = world_pts[1];

  setPerformPick (true);

  (void)glRenderMode(GL_SELECT);

  glInitNames();

  glPushName((GLuint) ~0);

  render();

  setPerformPick (false);

  glRenderMode (GL_RENDER);
  }

}  // PmGraphics



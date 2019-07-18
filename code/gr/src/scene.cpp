
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
//* scene:               s c e n e                             *
//*============================================================*

#include "pm/gr/scene.h"
#include "scene_prv.h"

namespace PmGraphics {


//*============================================================*
//*==========      constructors / detructor          ==========*
//*============================================================*

GrScene::GrScene() : name(NULL)
  {
  init();
  }


GrScene::GrScene(const string name)
  {
  this->name = name;
  //fprintf (stderr, "\n>>>>>>> GrScene::GrScene:  name [%s] \n", name.c_str());
  init();
  }

//*============================================================*
//*==========                init                    ==========*
//*============================================================*
// initialize scene data.

void
GrScene::init() 
  {
  ScenePrvData *scene_data = new ScenePrvData;
  prv_data = scene_data; 
  ScenePrvData *prvd = scene_data; 

  GrVector3 min(0.0, 0.0, 0.0);
  GrVector3 max(1.0, 1.0, 1.0);
  extent.set(min, max);

  window = NULL;

  // initialize pick data
  prvd->pick.perform = false;
  prvd->pick.active = false;
  prvd->pick.atts.type = GR_PICK_UNKNOWN;
  prvd->pick.return_func = NULL; 
  prvd->pick.result = NULL; 
  prvd->pick.requested = false; 
  prvd->pick.num = 0; 
  prvd->pick.num_requested = 0; 

  prvd->init_renderer = false;
  //api_createRenderer ();

  } 

GrScene::~GrScene() { 
  fprintf (stderr, "\n>>>>>>> GrScene::  dtor  name [%s] \n", name.c_str());
  }

//////////////////////////////////////////////////////////////
//                    p u b l i c                          //
////////////////////////////////////////////////////////////

//*============================================================*
//*==========                 render                 ==========*
//*============================================================*
// render the scene into its window.

void
GrScene::render () {

  if (!this->window) {  
    return;
    }

  // set the current graphics window //

  this->window->setCurrentWindow();  

  // set lights //

  this->positionLights();  

  // render the geometries in the scene //

  api_render ();

  // write window image to file //

  this->window->record();  

  }

//*============================================================*
//*==========                 setExtent              ==========*
//*============================================================*
// set scene geometric extent.

void
GrScene::setExtent (GrExtent& extent) 
  {
  /*
  fprintf (stderr, "\n>>>>>> GrScene::setExtent \n");
  fprintf (stderr, "       x = %f %f \n", extent.min[0], extent.max[0]);
  fprintf (stderr, "       y = %f %f \n", extent.min[1], extent.max[1]);
  fprintf (stderr, "       z = %f %f \n", extent.min[2], extent.max[2]);
  */
  this->extent = extent; 

  for (int i = 0; i < 3; i++) {
   this->xform.center[i] = (extent.max[i] + extent.min[i]) / 2.0;
   }

  /*
  fprintf (stderr, ">>> center  %f %f %f \n", this->xform.center[0], 
                                              this->xform.center[1],  
                                              this->xform.center[2]); 
  */

  ((ScenePrvData*)prv_data)->viewport.extent = extent;
  api_setViewport ();
  }

//*============================================================*
//*==========                 addGeometry            ==========*
//*============================================================*
// add geometry to the scene.

//GrScene::addGeometry (GrGeometry *geom) {

void
GrScene::addGeometry (GrVertexGeometry *geom) {
  this->geometry.push_back (geom);
  geom->setScene(this);
  }

//*============================================================*
//*==========                 addLight               ==========*
//*============================================================*
// add a light to the scene.

void
GrScene::addLight (GrLight *light) {
  this->lights.push_back (light);
  }

//*============================================================*
//*==========            positionLights              ==========*
//*============================================================*
// position lights.                

void
GrScene::positionLights () {
  int n = this->lights.size();

  for (int i = 0; i < n; i++) {
    GrLight *light = lights[i];
    light->enable();
    }
  }

//*============================================================*
//*==========            initPick                    ==========*
//*============================================================*
// initialize picking operation foar a scene.

void
GrScene::initPick ()
  {
  ScenePrvData *prvd = (ScenePrvData*)prv_data;
  prvd->pick.requested = false;
  prvd->pick.num = 0;
  prvd->pick.num_requested = 0;
  }

//*============================================================*
//*==========            setPickActive               ==========*
//*============================================================*
// set picking to be active.                

void
GrScene::setPickActive (bool flag)
  {
  ScenePrvData *prvd = (ScenePrvData*)prv_data;
  prvd->pick.active = flag;
  }

//*============================================================*
//*==========            setPickAtts                 ==========*
//*============================================================*
// initialize picking operation for a scene.

void
GrScene::setPickAtts (GrPickAtts& atts)
  {
  ScenePrvData *prvd = (ScenePrvData*)prv_data;
  prvd->pick.atts.tol = atts.tol;
  prvd->pick.atts.type = atts.type;
  }

//*============================================================*
//*==========            finishPick                  ==========*
//*============================================================*
// finish a picking operation for a scene.

void
GrScene::finishPick (GrPickResult& pick_res)
  {
  ScenePrvData *prvd = (ScenePrvData*)prv_data;

  // reduce pick list to the one closest to the eye

  GrVector3 pt1 = pick_res.line[0];
  GrPickGeom *min_pick = NULL; 

  if (pick_res.geom_list.size() == 0) {
    initPick();
    return;
    }

  if (prvd->pick.return_func) {
    (*prvd->pick.return_func)(this, pick_res);
    }

  initPick();
  }

//*============================================================*
//*==========            performPick                 ==========*
//*============================================================*
// perform a picking operation on the scene.

void
GrScene::performPick (int screen_pt[2], GrPickResult& pick) 
  {

  ScenePrvData *prvd = (ScenePrvData*)prv_data;

  #ifdef dbg_GrScene_performPick 
  fprintf (stderr, "\n>>>>>> GrScene::performPick: \n");
  fprintf (stderr, "   >>> requested [%d] \n", prvd->pick.requested);
  fprintf (stderr, "   >>> active    [%d] \n", prvd->pick.active);
  fprintf (stderr, "   >>> pick.atts.type [%d] \n", prvd->pick.atts.type); 
  #endif

  if (!prvd->pick.requested && !prvd->pick.active) {
    return;
    }

  prvd->pick.screen_pt[0] = screen_pt[0]; 
  prvd->pick.screen_pt[1] = screen_pt[1]; 

  api_pick (screen_pt, pick);
  pick.geom_list.clear();
  pick.atts.type = prvd->pick.atts.type;
  pick.atts.tol = prvd->pick.atts.tol;

  for (int i = 0; i < this->geometry.size(); i++) {
    GrVertexGeometry *geom = geometry[i];
    geom->pick(pick);
    }

  if (prvd->pick.requested) { 
    prvd->pick.num += 1;

   if (prvd->pick.num == prvd->pick.num_requested) {
      finishPick (pick);
      }
    }
  else {
    finishPick (pick);
    }
  }

//*============================================================*
//*==========            requestPick                 ==========*
//*============================================================*
// request a pick operation on the scene.

void
GrScene::requestPick (int num, GrPickType type, GrPickReturnFunc func)
  {
  ScenePrvData *prvd = (ScenePrvData*)prv_data;

  prvd->pick.requested = true;
  prvd->pick.num_requested = num;
  prvd->pick.num = 0;
  prvd->pick.atts.type = type;
  prvd->pick.return_func = func;
  }

//*============================================================*
//*==========            setPerformPick              ==========*
//*============================================================*
// set the scene perform pick flag.

void
GrScene::setPerformPick (const bool flag) { 
  ((ScenePrvData*)prv_data)->pick.perform = flag;
  }

//*============================================================*
//*==========            getPerformPick              ==========*
//*============================================================*
// get the scene perform pick flag.

bool 
GrScene::getPerformPick() {
  return ((ScenePrvData*)prv_data)->pick.perform;
  }

//*============================================================*
//*==========                 setViewport            ==========*
//*============================================================*
// set the viewport for the scene.

void
GrScene::setViewport (const int& win_width, const int& win_height, const GrExtent& extent)
  {
  ((ScenePrvData*)prv_data)->viewport.win_width = win_width;
  ((ScenePrvData*)prv_data)->viewport.win_height = win_height;
  ((ScenePrvData*)prv_data)->viewport.extent = extent;

  api_setViewport ();
  }

//*============================================================*
//*==========                 setWindow              ==========*
//*============================================================*
// set the scene's window.

void
GrScene::setWindow(GrWindow *win) {
  this->window = win; 
  }


//*============================================================*
//*==========                 backXformPoint         ==========*
//*============================================================*
// back transform a screen point into a couple of points in
// world coordinates.

void 
GrScene::backXformPoint (int screen_pts[2], GrVector3 world_pts[2]) {
  api_unprojectPoint (screen_pts, world_pts);
  }

//*============================================================*
//*==========                 rotate                 ==========*
//*============================================================*
// rotate the scene.

void
GrScene::rotate (bool inc, GrVector3& rot)
  {

  //fprintf (stderr, "\n---------- GrScene::rotate ---------- \n");
  //fprintf (stderr, ">>>>>> scene rot (%g %g %g) \n", rot[0], rot[1], rot[2]); 

  if (inc) {
    this->xform.angles[0] += rot[0];
    this->xform.angles[1] += rot[1];
    this->xform.angles[2] += rot[2];
    }
  else {
    this->xform.angles[0] = rot[0];
    this->xform.angles[1] = rot[1];
    this->xform.angles[2] = rot[2];
    }

  this->render();
  }

//*============================================================*
//*==========                 scale                  ==========*
//*============================================================*
// scale the scene.

void
GrScene::scale (bool inc, float scale)
  {

  //fprintf (stderr, "\n---------- GrScene::scale ---------- \n");
  //fprintf (stderr, ">>>>>> scene scale [%g] \n", scale); 

  if (inc) {
    this->xform.scale[0] *= scale;
    }

  else {
    this->xform.scale[0] = scale;
    }

  this->render();
  }

//*============================================================*
//*==========                 translate              ==========*
//*============================================================*
// translate the scene.

void
GrScene::translate (bool inc, GrVector3& trans)
  {

  /*
  fprintf (stderr, "\n---------- GrScene::translate ---------- \n");
  fprintf (stderr, ">>>>>> inc [%d] \n", inc);
  fprintf (stderr, ">>>>>> trans (%g %g %g) \n", trans[0], trans[1], trans[2]);
  fprintf (stderr, ">>>>>> scene trans (%g %g %g) \n",
     this->xform.translation[0], this->xform.translation[1] , this->xform.translation[2]);
  */

  if (inc) {
    this->xform.translation[0] += trans[0];
    this->xform.translation[1] += trans[1];
    this->xform.translation[2] += trans[2];
    }
  else {
    this->xform.translation[0] = trans[0];
    this->xform.translation[1] = trans[1];
    this->xform.translation[2] = trans[2];
    }

  this->render();
  }

//*============================================================*
//*==========                 getXform               ==========*
//*============================================================*
// get the xform for a scene.

void 
GrScene::getXform(GrXform& grxform)
  {
  fprintf (stderr, "\n---------- GrScene::getXform ---------- \n");
  grxform.center = this->xform.center;
  grxform.translation = this->xform.translation;
  grxform.angles = this->xform.angles;
  grxform.scale = this->xform.scale;
  }


}


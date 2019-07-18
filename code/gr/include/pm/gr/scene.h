
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
//* scene:                  s c e n e                          *
//*============================================================*

#ifndef _SCENE_GR_H_
#define _SCENE_GR_H_

#include "pm/gr/gr.h"
#include "pm/gr/win.h"
#include "pm/gr/geom.h"
#include "pm/gr/vgeom.h"
#include "pm/gr/light.h"

namespace PmGraphics {

class GrSceneXform {
  public:
    GrSceneXform(){ xinit(); };
    bool set;
    GrVector3 angles;
    GrVector3 center;
    GrVector3 scale;
    GrVector3 translation;
    GrMatrix3x3 matrix;

   void xinit() {
     set = false;
     scale.set (1, 1, 1);
     angles.set (0, 0, 0);
     translation.set (0, 0, 0);
     center.set (0, 0, 0);
     matrix.diag(1.0);
     }
  };


// GrScene
// ------------

 class PM_EXPORT GrScene {

   public:

     GrScene();
     GrScene(const string name); 
     ~GrScene();

     void addGeometry (GrVertexGeometry *geom);
     //void addGeometry (GrGeometry *geom);
     void addLight(GrLight *light);

     void setViewport(const int& win_width, const int& win_height, 
                      const GrExtent& extent);
     void setWindow(GrWindow *win);
     void setExtent (GrExtent& extent);

     void rotate(bool inc, GrVector3& rot);
     void getRotationAngles(GrVector3& angles) { angles = this->xform.angles; }
     void scale(bool inc, float scale);
     void translate(bool inc, GrVector3& trans);

     void performPick (int screen_pts[2], GrPickResult& pick); 
     void requestPick (int num, GrPickType type, GrPickReturnFunc func);
     void setPickActive (bool flag);
     void setPickAtts (GrPickAtts& atts); 
     void setPerformPick (const bool flag); 
     bool getPerformPick (); 

     void getXform(GrXform& xform);
     void setXformCenter(GrVector3& center) { this->xform.center = center; }
     void resetXform() { this->xform.xinit(); this->render(); }
     void backXformPoint (int screen_pts[2], GrVector3 world_pts[2]);

     void render();

   private:

     string name;
     GrExtent extent;
     GrSceneXform xform;
     int active;
     GrRenderMode render_mode;
     GrWindow *window;

     vector<class GrLight*> lights;
     vector<class GrVertexGeometry*> geometry;
     //vector<class GrGeometry*> geometry;

     void *prv_data;

     void init();
     void positionLights();
     void initPick();
     void finishPick(GrPickResult& pick);

     // graphics api specific functions
     void api_createRenderer();
     void api_pick (int screen_pt[2], GrPickResult& pick);
     void api_render();
     void api_setViewport();
     void api_unprojectPoint (int *screen_pts, GrVector3 world_pts[2]);

   };

}


#endif




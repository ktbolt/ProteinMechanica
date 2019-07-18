
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
//* gr:                 g r a p h i c s                        *
//*============================================================*

#ifndef _GRAPHICS_PRV_PM_H_ 
#define _GRAPHICS_PRV_PM_H_

#include "pm/pm.h"

#ifdef USE_GRAPHICS

#include "pm/gr/gr.h"
#include "pm/gr/win.h"
#include "pm/gr/scene.h"
#include "pm/gr/cyl.h"
#include "pm/gr/ellipsoid.h"
#include "pm/gr/line.h"
#include "pm/gr/poly.h"
#include "pm/gr/point.h"
#include "pm/gr/sphere.h"

using namespace PmGraphics;

typedef struct GraphicsContext {
  GrWindow *window;
  GrScene *scene;
  //GrGeometry *geom;
  GrVertexGeometry *geom;
  GrVertexGeometry *geom1;
  PmGraphicsGeometry *pmgr_geom;
  vector<GrVertexGeometry*> geom_list;
  } GraphicsContext;

typedef struct GraphicsPrvData {
  GrWindow *window;
  GrScene *scene;
  } GraphicsPrvData;

#else

typedef struct GraphicsContext {
  void *window;
  } GraphicsContext;

typedef struct GraphicsPrvData {
  void *window;
  } GraphicsPrvData;

#endif

#endif




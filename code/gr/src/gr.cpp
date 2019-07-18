
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
//* gr:                 s y s t e m                            *
//*============================================================*

#include "pm/gr/gr.h"

namespace PmGraphics {

GrSystem grSystem;

//*============================================================*
//*==========      constructors / detructor          ==========*
//*============================================================*

GrSystem::GrSystem() {
    //fprintf (stderr, ">>>>>> GrSystem: init \n");
    }

//*============================================================*
//*==========              setCommandCallback        ==========*
//*============================================================*

void 
GrSystem::setCommandCallback1 (void (*func)(const string cmd)) 
  {
  cmd_callback_func1 = func; 
  }

void 
GrSystem::setCommandCallback (void (*func)(FILE *fp, const string cmd)) 
  {
  cmd_callback_func = func; 
  }

//*============================================================*
//*==========              processCommand            ==========*
//*============================================================*
// process a command from the keyboard.

void
GrSystem::processCommand (const string cmd)
  {
  //fprintf (stderr, ">>>>>> GrSystem::processCommand \n");

  if (cmd_callback_func) { 
    //fprintf (stderr, ">>> cmd_callback_func \n");
    (*cmd_callback_func)(0, cmd); 
    }

  if (cmd_callback_func1) { 
    //fprintf (stderr, ">>> cmd_callback_func1 \n");
    (*cmd_callback_func1)(cmd); 
    }
  }

//*============================================================*
//*==========              addWindow                 ==========*
//*============================================================*

void
GrSystem::addWindow(GrWindow *win) 
  {
  windows.push_back (win);
  }

//*============================================================*
//*==========        gr_RotationFromAngles           ==========*
//*============================================================*

void
gr_RotationFromAngles (float xr, float yr, float zr, PmMatrix3x3& mat)
  {

  float sx, sy, sz, cx, cy, cz, f;

  f = (4.0 * atan(1.0)) / 180.0;
  sx = sin(f*xr);
  sy = sin(f*yr);
  sz = sin(f*zr);

  cx = cos(f*xr);
  cy = cos(f*yr);
  cz = cos(f*zr);

  mat(0,0) = cy * cz;
  mat(0,1) = -cy * sz;
  mat(0,2) = sy;

  mat(1,0) = sx * sy * cz + cx * sz;
  mat(1,1) = -sx * sy * sz + cx * cz;
  mat(1,2) = -sx * cy;

  mat(2,0) = -cx * sy * cz + sx * sz;
  mat(2,1) = cx * sy * sz + sx * cz;
  mat(2,2) = cx * cy;
  }

//*============================================================*
//*==========           gr_getxformMatrix            ==========*
//*============================================================*

void
gr_getxformMatrix (GrXform& xform, PmMatrix3x3& mat) {
  /*
  if (!xform.matrix_set) {
    gr_RotationFromAngles (xform.angles[0], xform.angles[1], xform.angles[2], mat);
    xform.setMatrix(mat);
    }
  */
  }

}



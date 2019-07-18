
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
//* common:          c o m m o n   t y p e s                   *
//*============================================================*

#include "pm/common.h"
#include "pm/mth.h"

using namespace std;

namespace ProteinMechanica {

////////////////////////////////////////////////////////////////
//                       x f o r m                           //
//////////////////////////////////////////////////////////////

//*==============================================================*
//*==========              setRotation                ===========*
//*==============================================================*
// rotation from three angles.

void 
PmXform::setRotation(float xr, float yr, float zr)
  {
  set = true;
  pm_MathRotationFromAngles (xr, yr, zr, matrix);
  this->angles[0] = xr;
  this->angles[1] = yr;
  this->angles[2] = zr;
  }

//*==============================================================*
//*==========              setRotation                ===========*
//*==============================================================*
// rotate around an axis.

void 
PmXform::setRotation(float angle, PmVector3& axis)
  {
  set = true;
  pm_MathRotationAroundAxis (angle, axis, matrix);
  }

};




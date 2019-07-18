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
//* light_opengl:       o p e n g l    l i g h t               *
//*============================================================*

#include "light_opengl.h"

namespace PmGraphics {


//*============================================================*
//*==========              api_setLight              ==========*
//*============================================================*
// create a renderer.

void
GrLight::api_setLight()
  {

  GLfloat params[4];

  int id = ((LightPrvData *)prv_data)->hw_id; 
  hw_id = get_hw_light_id (id); 

  glDisable (GL_CULL_FACE);

  this->ambient.get(params);
  glLightfv (hw_id, GL_AMBIENT, params); 

  this->diffuse.get(params);
  glLightfv (hw_id, GL_DIFFUSE, params);

  this->specular.get(params);
  glLightfv (hw_id, GL_SPECULAR, params);

  glLightfv (hw_id, GL_POSITION, (GLfloat*)&this->geometry[0]);

  glEnable (hw_id);

  //glLightfv (hw_id, GL_POSITION, (GLfloat*)&this->geometry[0]);
  }


}



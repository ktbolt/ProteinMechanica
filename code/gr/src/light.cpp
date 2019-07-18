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
//* light:               l i g h t                             *
//*============================================================*

#include "pm/gr/light.h"
#include "light_prv.h"

namespace PmGraphics {


//*============================================================*
//*==========      constructors / detructor          ==========*
//*============================================================*

GrLight::GrLight() : name(NULL)
  {
  init();
  }

GrLight::GrLight(const char *name, GrLightType type)
  {
  this->name = name;
  this->type = type; 
  //fprintf (stderr, "\n>>>>>>> GrLight::GrLight:  name [%s] \n", this->name.c_str());
  init();
  }


// initialize light data.

void
GrLight::init() 
  {

  static float diffuse_vals[4] = {0.4, 0.4, 0.4, 1.0};
  static float ambient_vals[4] = {0.2, 0.2, 0.2, 1.0};

  diffuse = diffuse_vals;
  ambient = ambient_vals;
  color.set(1.0, 1.0, 1.0);

  if (type == GR_LIGHT_POSITIONAL) {
    geometry[0] = 0.0;
    geometry[1] = 0.0;
    geometry[2] = 1.0;
    geometry[3] = 0.0;
    }

  else if (GR_LIGHT_DIRECTIONAL) { 
    geometry[0] = 0.0;
    geometry[1] = 0.0;
    geometry[2] = 1.0;
    geometry[3] = 0.0;
    }
  
  LightPrvData *light_data = new LightPrvData;
  light_data->hw_id = num_lights++; 
  light_data->hw_on = false;
  light_data->hw_def = false;
  prv_data = light_data; 
  } 

GrLight::~GrLight() { 
  }


//////////////////////////////////////////////////////////////
//                    p u b l i c                          //
////////////////////////////////////////////////////////////

void
GrLight::enable() {
  api_setLight();
  }


}


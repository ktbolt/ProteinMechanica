
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
//* model:                    m o d e l                        *
//*============================================================*

#ifndef _MODEL_PM_H_ 
#define _MODEL_PM_H_

#include "pm/pm.h"
#include "body.h"
#include "joint.h"

namespace ProteinMechanica {

// PmModel 
// ---------------------
// Model  

 class PM_EXPORT PmModel {

   public:
      PmModel (const string name);

      void displayBodies(bool show);
      void setBodyMsize(float msize) { body_msize = msize; }

      void addJoint(PmJoint *joint);
      void displayJoints(bool show);
      void getJoint(const string name, PmJoint **joint);
      void setJointMsize(float msize) { joint_msize = msize; }
      void getName (string& name) { name = this->name;  }
      bool hasName(const string& name) { return name == this->name;  }

      bool addBody(PmBody *body);
      void getBody(const string name, PmBody **bdy);

   private:
      string name;
      vector <PmBody*> bodies;
      PmGraphicsModelBody *body_geom;
      float body_msize;

      vector <PmJoint*> joints;
      PmGraphicsModelJoint *joint_geom;
      float joint_msize;

      vector <PmPhysicalObj*> physical_objects;

   };


}

#endif




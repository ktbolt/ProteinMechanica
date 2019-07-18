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
//* model:             m o d e l                               *
//*============================================================*

#include "model.h"

#define ndbg_PmModel

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmModel::PmModel(const string name) 
  { 
  //fprintf (stderr, ">>>>>> PmModel:: ctor \n");
  this->name = name;
  joint_msize = 0.2;
  body_msize = 0.2;
  }

//*============================================================*
//*==========              addBody                   ==========*
//*============================================================*
// add a body to the model.

bool
PmModel::addBody(PmBody *body)
  {
  string name;
  PmBody *bdy;
  body->getName (name);
  getBody (name, &bdy);

  if (!bdy) {
    bodies.push_back(body);
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getBody                  ==========*
//*============================================================*

void
PmModel::getBody (const string name, PmBody **bdy)
  {
  string str;
  *bdy = NULL;

  for (unsigned int i = 0; i < bodies.size(); i++) {
    PmBody *body = bodies[i];
    body->getName (str);

    if (str == name) {
      *bdy = body;
      }
    }
  }

//*============================================================*
//*==========              displayBodies             ==========*
//*============================================================*
// display bodies.

void
PmModel::displayBodies (bool show)
  {
  #ifdef dbg_PmModel
  fprintf (stderr, ">>>>>> PmModel::displayBodies  \n");
  #endif
  if (!pmSystem.useGraphics()) return;
  int num_bodies = bodies.size();

  if (!num_bodies) {
    return;
    }

  PmMassProperties mass_props;
  PmVector3 *coords = new PmVector3[num_bodies];

  for (int i = 0; i < num_bodies; i++) {
    PmBody *body = bodies[i];
    body->getMassProps (mass_props);
    coords[i] = mass_props.com;
    }

  string geom_name;
  PmGraphicsAttributes atts;
  geom_name = "model[" + name + "]bodies";
  body_geom = new PmGraphicsModelBody (geom_name, num_bodies, coords);
  atts.setScale(body_msize);
  body_geom->setAttributes(atts);
  body_geom->display();
  }

//*============================================================*
//*==========              addJoint                  ==========*
//*============================================================*
// add a joint to the simulation.

void
PmModel::addJoint (PmJoint *joint)
  {
  joints.push_back (joint);
  }

//*============================================================*
//*==========              displayJoints             ==========*
//*============================================================*
// display joints.                       

void 
PmModel::displayJoints(bool show)
  {
  if (!pmSystem.useGraphics()) return;
  int num_joints = joints.size(); 
 
  if (!num_joints) {
    return;
    }

  PmVector3 *coords = new PmVector3[num_joints];

  for (int i = 0; i < num_joints; i++) {
    PmJoint *joint = joints[i];
    joint->getPosition (coords[i]);
    }

  PmGraphicsAttributes atts;
  string geom_name;
  geom_name = "model[" + name + "]joints";
  joint_geom = new PmGraphicsModelJoint (geom_name, num_joints, coords);
  atts.setScale(joint_msize);
  joint_geom->setAttributes(atts);
  joint_geom->display();
  }

//*============================================================*
//*==========              getJoint                  ==========*
//*============================================================*

void
PmModel::getJoint (const string name, PmJoint **jnt)
  {
  string str;
  *jnt = NULL;

  for (unsigned int i = 0; i < joints.size(); i++) {
    PmJoint *joint = joints[i];
    joint->getName (str);

    if (str == name) {
      *jnt = joint;
      }
    }
  }

}

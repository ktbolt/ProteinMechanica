
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
// * joint:                  j o i n t                         *
//*============================================================*

#include "joint.h"

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmJoint::PmJoint(const string name, const PmJointType type) 
  { 
  //fprintf (stderr, ">>>>>> PmJoint: ctor   name [%s] \n", name.c_str());
  this->name = name;
  this->type = type;
  body1 = NULL;
  body2 = NULL;
  position.set(0,0,0);
  current_position.set(0,0,0);
  color.set(1,1,1);
  msize = 0.1;
  position_set = false;
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  show = false;
  show_axes = false;

  file_ptr = NULL;
  binary_format = false;
  write_init = true;
  last_written = 0;
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// creates a joint of a particular type.

PmJoint*
PmJoint::create (const string name, const PmJointType type)
  {

  PmJoint *joint = NULL; 

  switch (type) {
    case PM_JOINT_HINGE:
      joint = new PmHingeJoint (name, type);
    break;

    case PM_JOINT_SLIDER:
    break;

    case PM_JOINT_GIMBAL:
    break;

    case PM_JOINT_BALL:
      joint = new PmBallJoint (name, type);
    break;

    case PM_JOINT_FREE:
      joint = new PmFreeJoint (name, type);
    break;

    case PM_JOINT_WELD:
      joint = new PmWeldJoint (name, type);
    break;

    case PM_JOINT_UNIVERSAL:
      joint = new PmUniversalJoint (name, type);
    break;

    default:
    break;
    }

  return joint;
  }

//*============================================================*
//*==========              getType                   ==========*
//*============================================================*
// get joint type string.

void 
PmJoint::getType(string& str)
  {
  switch (type) {
    case PM_JOINT_HINGE:
      str = "hinge"; 
    break;

    case PM_JOINT_SLIDER:
      str = "slider"; 
    break;

    case PM_JOINT_GIMBAL:
      str = "gimbal"; 
    break;

    case PM_JOINT_BALL:
      str = "ball"; 
    break;

    case PM_JOINT_FREE:
      str = "free"; 
    break;

    case PM_JOINT_WELD:
      str = "weld"; 
    break;

    case PM_JOINT_UNIVERSAL:
      str = "universal"; 
    break;

    default:
      str = "unknown"; 
    }
  }

//*============================================================*
//*==========              convJointType             ==========*
//*============================================================*
// convert string name into a joint type enum.

void
PmJoint::convJointType(const string str, PmJointType& type) 
  {

  if ((str == "pin") || (str == "hinge")) {
    type = PM_JOINT_HINGE;
    }
  else if (str == "slider") {
    type = PM_JOINT_SLIDER;
    }
  else if (str == "slider") {
    type = PM_JOINT_SLIDER;
    }
  else if (str == "gimbal") {
    type = PM_JOINT_GIMBAL;
    }
  else if (str == "ball") {
    type = PM_JOINT_BALL;
    }
  else if (str == "free") {
    type = PM_JOINT_FREE;
    }
  else if (str == "weld") {
    type = PM_JOINT_WELD;
    }
  else if (str == "universal") {
    type = PM_JOINT_UNIVERSAL;
    }
  else {
    type = PM_JOINT_UNKNOWN;
    }
  }

//*============================================================*
//*==========              motors                    ==========*
//*============================================================*
// get/add motors.                               

void 
PmJoint::addMotor(PmMotor *motor) {
  motors.push_back(motor);
  }

void PmJoint::getMotors(vector<PmMotor*>& motors) { motors = this->motors; }

bool PmJoint::hasMotors() { return motors.size() != 0; }

//*============================================================*
//*==========          get/setPosition               ==========*
//*============================================================*

void 
PmJoint::getCurrentPosition (PmVector3& pos) { 
  pos = current_position; 
  }

void 
PmJoint::getPosition (PmVector3& pos) { 
  pos = position; 
  }

void 
PmJoint::setPosition (const PmVector3& pos) 
  { 
  position = pos; 
  current_position = pos; 
  position_set = true; 
  }

//*============================================================*
//*==========          setWriteResultsParams         ==========*
//*============================================================*
// set writting results.

void
PmJoint::setWriteResultsParams(const string file_name, const bool binary)
  {
  this->write_file_name = file_name;
  this->write_init = true;
  this->last_written = 0;
  }

//*============================================================*
//*==========               addEnergy                ==========*
//*============================================================*
// add energy.                      

void 
PmJoint::addEnergy(const bool upd, const float val)
  {
  if (upd) {
    energy.push_back(val);
    }
  }

void
PmJoint::getEnergy(const int step, float& pe)
  {
  pe = 0.0;

  if (!energy.size()) {
    return;
    }

  if (step >= (int)energy.size()) {
    pm_ErrorWarnReport (PM, "step %d too large for getting pe for joint \"%s\" .",
                        "*", step, name.c_str());
    return;
    }

  pe = energy[step];
  }

///////////////////////////////////////////////////////////////
//                b a l l   j o i n t                       //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize data particular to a ball joint. 

void
PmBallJoint::init()
  {
  //fprintf (stderr, ">>>>>> PmBallJoint: init name [%s] \n", name.c_str());
  geometry = NULL;
  axes[0].set (1, 0, 0);
  axes[1].set (0, 1, 0);
  axes[2].set (0, 0, 1);

  for (int i = 0; i < 3; i++) {
    force_const[i] = 0.0;
    range[i] = 0.0;
    range_set[i] = false;
    }
  }

//*============================================================*
//*==========              axis                      ==========*
//*============================================================*
// set and get axes.

void 
PmBallJoint::getAxis(const int id, PmVector3& a) 
  { 
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  a = this->axes[i-1]; 
  }

void
PmBallJoint::setAxes(const PmVector3& a1, const PmVector3& a2)
  {
  float dp;
  PmVector3 a3;

  // check for orthogonal axes //

  dp = a1 * a2;

  if (fabs(dp) > 1e-6) {
    pm_ErrorWarnReport (PM, "ball joint \"%s\" axes are not orthogonal.", "*", 
                        name.c_str());
    }

  a3 = a1.cross(a2);
  this->axes[0] = a1;
  this->axes[1] = a2;
  this->axes[2] = a3;
  }

//*============================================================*
//*==========              ForceConst                ==========*
//*============================================================*
// set and get force constants.

void
PmBallJoint::getForceConst(const int id, float& val)
  {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  val = this->force_const[i-1];
  }

void
PmBallJoint::setForceConst(const int id, const float val)
  {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  this->force_const[i-1] = val;
  }

void
PmBallJoint::getForceConsts(vector<float>& k)
  {
  k.clear();
  for (int i = 0; i < 3; i++) {
    k.push_back(force_const[i]);
    }
  }

//*============================================================*
//*==========              Range                     ==========*
//*============================================================*
// set and get rotational range.

void
PmBallJoint::getRange(const int id, float& val)
  {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  val = this->range[i-1];
  }

bool 
PmBallJoint::hasRange(const int id)
  {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  return this->range_set[i-1];
  }

void
PmBallJoint::setRange(const int id, const float val)
  {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  this->range[i-1] = val;
  this->range_set[i-1] = true;
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display a ball joint.

void
PmBallJoint::display(const bool show)
  {
  /*
  fprintf (stderr, ">>>>>> PmBallJoint::display [%s] \n", name.c_str());
  fprintf (stderr, "   >>> color (%g %g %g) \n", color[0], color[1], color[2]); 
  fprintf (stderr, "   >>> pos   (%g %g %g) \n", position[0], position[1], position[2]); 
  fprintf (stderr, "   >>> msize = %g \n", msize);
  fprintf (stderr, "   >>> show = %d \n", show);
  */

  if (!pmSystem.useGraphics()) return;

  if (geometry == NULL) {
    string geom_name; 
    geom_name = "joint[ball:" + name + ']';
    //fprintf (stderr, "   >>> geom name [%s] \n", geom_name.c_str()); 
    PmVector3 *verts = new PmVector3[3]; 
    verts[0] = position;
    verts[1] = position;
    verts[2] = position;

    PmVector3 axes[3];
    axes[0] = this->axes[0];
    axes[1] = this->axes[1];
    axes[2] = this->axes[2];

    geometry = new PmGraphicsBallJoint (geom_name, verts, axes);
    }

  this->show = show;

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setScale(msize);
  atts.setShadingType(shading_type); 
  atts.setDisplayType(display_type); 
  atts.setVisibleAux(show_axes);
  atts.setVisible(show);
  geometry->setAttributes(atts); 
  geometry->display(); 
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update a ball joint.

void
PmBallJoint::update(const PmVector3& pos)
  {
  if (!geometry) {
    return;
    }

  //fprintf (stderr, ">>>>>> PmBallJoint::update[%s] \n", name.c_str());
  //fprintf (stderr, "   >>> pos   (%g %g %g) \n", pos[0], pos[1], pos[2]); 
  PmXform xform;
  xform.set = true; 
  xform.translation = pos - position;
  geometry->xform(xform);
  //geometry->update(pos);
  current_position = pos;
  }

void
PmBallJoint::update(PmXform& body_xform)
  {
  if (!geometry) {
    return;
    }

  PmXform xform;
  xform.set = true;
  xform.translation = body_xform.translation - position;

  xform.center = body_xform.center;
  xform.center = position;
  xform.matrix = body_xform.matrix;

  geometry->xform(xform);
  //geometry->update(pos);
  }

///////////////////////////////////////////////////////////////
//                f r e e   j o i n t                       //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize data particular to a free joint.

void
PmFreeJoint::init()
  {
  //fprintf (stderr, ">>>>>> PmFreeJoint: init name [%s] \n", name.c_str());
  geometry = NULL;
  }

void
PmFreeJoint::getForceConsts(vector<float>& k) {
  k.clear();
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display a free joint.

void
PmFreeJoint::display(const bool show)
  {
  /*
  fprintf (stderr, ">>>>>> PmFreeJoint::display [%s] \n", name.c_str());
  fprintf (stderr, "   >>> color (%g %g %g) \n", color[0], color[1], color[2]);
  fprintf (stderr, "   >>> pos   (%g %g %g) \n", position[0], position[1], position[2]);
  */
  if (!pmSystem.useGraphics()) return;

  if (geometry == NULL) {
    string geom_name;
    geom_name = "joint[free:" + name + ']';
    //fprintf (stderr, "   >>> geom name [%s] \n", geom_name.c_str());
    PmVector3 *verts = new PmVector3;
    verts[0] = position;
    geometry = new PmGraphicsFreeJoint (geom_name, verts);
    }

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setScale(msize);
  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  geometry->setAttributes(atts);
  geometry->display();
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update a free joint.

void
PmFreeJoint::update(const PmVector3& pos) {
  if (!geometry) {
    return;
    }

  geometry->update(pos);
  }

///////////////////////////////////////////////////////////////
//                h i n g e    j o i n t s                  //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize data particular to a hinge joint.

void
PmHingeJoint::init()
  {
  //fprintf (stderr, ">>>>>> PmHingeJoint: init name [%s] \n", name.c_str());
  geometry = NULL;
  axis.set (0, 0, 1);
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  force_const = 0.0;
  range = 360.0;
  range_set = false;
  position2_set = false;
  }

//*============================================================*
//*==========              axis                      ==========*
//*============================================================*
// set and get axes.             

void 
PmHingeJoint::getAxis(PmVector3& vec) { 
  vec = this->axis; 
  }

void 
PmHingeJoint::setAxis(const PmVector3& vec) 
  { 
  this->axis = vec; 
  this->axis.normalize();
  }

//*============================================================*
//*==========          get/setPosition2              ==========*
//*============================================================*

void
PmHingeJoint::getCurrentPosition2(PmVector3& pos) {
  pos = current_position2;
  }

void
PmHingeJoint::getPosition2 (PmVector3& pos) {
  pos = position2;
  }

void
PmHingeJoint::setPosition2 (const PmVector3& pos)
  {
  position2 = pos;
  current_position2 = pos;
  position2_set = true;
  }

//*============================================================*
//*==========              force const               ==========*
//*============================================================*
// set and get force constants.

void PmHingeJoint::getForceConst(float& k) { k = force_const; }
void PmHingeJoint::setForceConst(const float k) { force_const = k; }
void PmHingeJoint::getForceConsts(vector<float>& k) { 
     k.clear(); k.push_back(force_const);
     }

//*============================================================*
//*==========              set/get angle range       ==========*
//*============================================================*

void 
PmHingeJoint::setRange(const float r)
  {
  range_set = true;
  range = r;
  }

bool 
PmHingeJoint::getRange(float& r) {
  r = range;
  return range_set;
  }

//*============================================================*
//*==========              updateAngle               ==========*
//*============================================================*
// set joint rotation angle. 

void 
PmHingeJoint::updateAngle(float a) {
  angles.push_back(a);
  }

//*============================================================*
//*==========              getAngle                  ==========*
//*============================================================*
// get joint rotation angle.

void 
PmHingeJoint::getAngle(float& a) 
  {
  int n = angles.size();

  if (!n) {
    a = 0.0;
    return;
    }

  a = angles[n-1];
  }

void 
PmHingeJoint::getAngles(vector<float>& a) {
  a = angles;
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display a hinge joint.

void
PmHingeJoint::display(const bool show)
  {
  /*
  fprintf (stderr, ">>>>>> PmHingeJoint::display [%s] \n", name.c_str());
  fprintf (stderr, "   >>> color (%g %g %g) \n", color[0], color[1], color[2]);
  fprintf (stderr, "   >>> pos   (%g %g %g) \n", position[0], position[1], position[2]);
  */
  if (!pmSystem.useGraphics()) return;

  if (geometry == NULL) {
    string geom_name;
    geom_name = "joint[hinge:" + name + ']';
    //fprintf (stderr, "   >>> geom name [%s] \n", geom_name.c_str());

    if (position2_set) {
      PmVector3 *verts = new PmVector3[2];
      verts[0] = position;
      verts[1] = position2;
      geometry = new PmGraphicsHingeJoint (geom_name, 2, verts, axis);
      }
    else {
      PmVector3 *verts = new PmVector3;
      verts[0] = position;
      geometry = new PmGraphicsHingeJoint (geom_name, 1, verts, axis);
      }
    }

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setScale(msize);
  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  geometry->setAttributes(atts);
  geometry->display();
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update a hinge joint.

void
PmHingeJoint::update(const PmVector3& pos)
  {
  if (!geometry) {
    return;
    }

  /*
  fprintf (stderr, ">>>>>> PmHingeJoint::update \n");
  fprintf (stderr, ">>> position=%f %f %f\n", position[0], position[1], position[2]); 
  fprintf (stderr, ">>> pos=%f %f %f\n", pos[0], pos[1], pos[2]); 
  */
  PmXform xform;
  xform.set = true; 
  xform.translation = pos - position;
  geometry->xform(xform);
  //geometry->update(pos);
  }

void
PmHingeJoint::update(PmXform& body_xform)
  {
  if (!geometry) {
    return;
    }

  PmXform xform;
  xform.set = true;
  xform.translation = body_xform.translation - position;

  xform.center = body_xform.center;
  xform.center = position; 
  xform.matrix = body_xform.matrix;

  geometry->xform(xform);
  //geometry->update(pos);
  }

//*============================================================*
//*==========              writeResults              ==========*
//*============================================================*
// write hinge joint results.

void 
PmHingeJoint::writeResults(const vector<float>& time, const int lastw, const int num) 
  {
  char fname[100];

  if (this->write_file_name.empty()) {
    return;
    }

  // write file header //

  if (this->write_init) {
    sprintf (fname, "%s.pm", this->write_file_name.c_str());
    this->file_ptr = fopen (fname, "w");

    if (!this->file_ptr) {
      return;
      }

    fprintf (this->file_ptr, "# Protein Mechanica joint data file \n");
    fprintf (this->file_ptr, "# time angle            \n");
    fprintf (this->file_ptr, "angle \"%s\" \n", name.c_str());
    this->write_init = false;
    }

  // write data //

  for (int n = lastw; n < num; n++) {
    fprintf (this->file_ptr, "%f %f \n", time[n], angles[n]);
    }

  fflush (this->file_ptr);
  }

///////////////////////////////////////////////////////////////
//           u n i v e r s a l    j o i n t s               //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize data particular to a universal joint.

void
PmUniversalJoint::init()
  {
  //fprintf (stderr, ">>>>>> PmUniversalJoint: init name [%s] \n", name.c_str());
  geometry = NULL;
  axis[0].set (0, 0, 1);
  axis[1].set (1, 0, 0);
  force_const[0] = 0.0;
  force_const[1] = 0.0;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  }

//*============================================================*
//*==========              axis                      ==========*
//*============================================================*
// set and get axes.             

void 
PmUniversalJoint::getAxis(const int id, PmVector3& a) 
  { 
  a.set(0,0,0);

  if ((id < 1) || (id > 2)) { 
    return;
    }

  a = this->axis[id-1]; 
  }

void 
PmUniversalJoint::setAxis(const int id, const PmVector3& a) 
  { 
  if ((id < 1) || (id > 2)) { 
    return;
    }

  this->axis[id-1] = a;
  }

//*============================================================*
//*==========              ForceConst                ==========*
//*============================================================*
// set and get force constants. 

void 
PmUniversalJoint::getForceConst(const int id, float& k) 
  {
  if ((id < 1) || (id > 2)) { 
    return;
    }

  k = force_const[id-1]; 
  }

void 
PmUniversalJoint::setForceConst(const int id, const float k) 
  {
  if ((id < 1) || (id > 2)) { 
    return;
    }
  force_const[id-1] = k; 
  }

void
PmUniversalJoint::getForceConsts(vector<float>& k)
  {
  k.clear();
  for (int i = 0; i < 2; i++) {
    k.push_back(force_const[i]);
    }
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display a universal joint.

void
PmUniversalJoint::display(const bool show)
  {
  #ifdef dbg_PmUniversalJoint
  fprintf (stderr, ">>>>>> PmUniversalJoint::display [%s] \n", name.c_str());
  fprintf (stderr, "   >>> color (%g %g %g) \n", color[0], color[1], color[2]);
  fprintf (stderr, "   >>> pos   (%g %g %g) \n", position[0], position[1], position[2]);
  fprintf (stderr, "   >>> axis1 (%g %g %g) \n", axis1[0], axis1[1], axis1[2]);
  fprintf (stderr, "   >>> axis2 (%g %g %g) \n", axis2[0], axis2[1], axis2[2]);
  #endif

  if (!pmSystem.useGraphics()) return;

  if (geometry == NULL) {
    string geom_name;
    geom_name = "joint[universal:" + name + ']';
    //fprintf (stderr, "   >>> geom name [%s] \n", geom_name.c_str());
    PmVector3 *verts = new PmVector3[2];
    verts[0] = position;
    verts[1] = position;
    PmVector3 gaxes[2];
    gaxes[0] = this->axis[0];
    gaxes[1] = this->axis[1];
    geometry = new PmGraphicsUniversalJoint (geom_name, verts, gaxes);
    }

  this->show = show;
  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setScale(msize);
  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setVisible(show);

  geometry->setAttributes(atts);
  geometry->display();
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update a universal joint.

void
PmUniversalJoint::update(const PmVector3& pos)
  {
  if (!geometry) {
    return;
    }

  PmXform xform;
  xform.set = true; 
  xform.translation = pos - position;
  geometry->xform(xform);
  //geometry->update(pos);
  }

void
PmUniversalJoint::update(PmXform& body_xform)
  {
  if (!geometry) {
    return;
    }

  PmXform xform;
  xform.set = true;
  xform.translation = body_xform.translation - position;

  xform.center = body_xform.center;
  xform.center = position;
  xform.matrix = body_xform.matrix;

  geometry->xform(xform);
  //geometry->update(pos);
  }

///////////////////////////////////////////////////////////////
//                w e l d   j o i n t                       //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize data particular to a free joint.

void
PmWeldJoint::init()
  {
  //fprintf (stderr, ">>>>>> PmWeldJoint: init name [%s] \n", name.c_str());
  geometry = NULL;
  }

void
PmWeldJoint::getForceConsts(vector<float>& k) {
  k.clear();
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display a weld joint.

void
PmWeldJoint::display(const bool show)
  {
  if (!pmSystem.useGraphics()) return;

  if (geometry == NULL) {
    string geom_name;
    geom_name = "joint[weld:" + name + ']';
    PmVector3 *verts = new PmVector3[1];
    verts[0] = position;
    geometry = new PmGraphicsWeldJoint (geom_name, verts);
    }

  this->show = show;

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setScale(msize);
  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setVisible(show);
  geometry->setAttributes(atts);
  geometry->display();
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update a weld joint.

void
PmWeldJoint::update(const PmVector3& pos)
  {
  if (!geometry) {
    return;
    }

  PmXform xform;
  xform.set = true;
  xform.translation = pos - position;
  geometry->xform(xform);
  }

void
PmWeldJoint::update(PmXform& body_xform)
  {
  if (!geometry) {
    return;
    }

  PmXform xform;
  xform.set = true;
  xform.translation = body_xform.translation - position;

  xform.center = body_xform.center;
  xform.center = position;
  xform.matrix = body_xform.matrix;

  geometry->xform(xform);
  }

}

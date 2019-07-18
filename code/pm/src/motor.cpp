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
//* motor:                    m o t o r                        *
//*============================================================*

#include "motor.h"
#include "graphics.h"

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmMotor::PmMotor() {
  file_ptr = NULL;
  binary_format = false;
  write_init = true;
  last_written = 0;
  }

//*============================================================*
//*==========            convMotorType               ==========*
//*============================================================*

void
PmMotor::convMotorType(const string str, PmMotorType& type)
  {
  if (str == "angular") {
    type = PM_MOTOR_ANGULAR;
    }
  else if (str == "linear") {
    type = PM_MOTOR_LINEAR;
    }
  else {
    type = PM_MOTOR_UNKNOWN;
    }
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a simulation object of a particular type.

PmMotor*
PmMotor::create(const string name, const PmMotorType type)
  {
  switch (type) {
    case PM_MOTOR_ANGULAR:
      {
      PmAngularMotor *amotor = new PmAngularMotor(name);
      amotor->type = type;
      return amotor;
      }
    break;

    default:
      return NULL;
    }
  }

//*============================================================*
//*==========               enabled                  ==========*
//*============================================================*
// set/get enabled flag. 

bool PmMotor::isEnabled() { return enabled; }

void PmMotor::setEnabled(const bool flag) { enabled = flag; }


//*============================================================*
//*==========          setWriteResultsParams         ==========*
//*============================================================*
// set writting results.

void
PmMotor::setWriteResultsParams(const string file_name, const bool binary)
  {
  //fprintf (stderr, ">>>>>> PmTrace::setWriteResultsParams \n");
  this->write_file_name = file_name;
  this->write_init = true;
  this->last_written = 0;
  }

//*============================================================*
//*==========               writeResults             ==========*
//*============================================================*
// write results.

void
PmMotor::writeResults(const vector<float>& time, const int lastw, const int num)
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

    fprintf (this->file_ptr, "# Protein Mechanica motor data file \n");
    fprintf (this->file_ptr, "# time angle              \n");
    fprintf (this->file_ptr, "motor \"%s\" \n", name.c_str());
    this->write_init = false;
    }

  // write data //

  for (int n = lastw; n < num; n++) {
    fprintf (this->file_ptr, "%f %f \n", time[n], data[n]);
    }

  fflush (this->file_ptr);
  }

///////////////////////////////////////////////////////////////
//               a n g u l a r   m o t o r                  //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmAngularMotor::PmAngularMotor(const string name) 
  {
  this->name = name;
  line_geometry = NULL;
  active = false;
  enabled = false;
  joint = NULL;
  num_axes = 0;

  for (int i = 0; i < 3; i++) {
    axes[i] = 0;
    max_velocity[i] = 0;
    max_force[i] = 0;
    max_angle[i] = 0;
    angle[i] = 0;
    angle_inc[i] = 0;
    }
  }

//*============================================================*
//*==========               isActive                 ==========*
//*============================================================*
// check if a motor is active for the given time.

bool 
PmAngularMotor::isActive(float t)
  {
  return time_interval.hasTime(t);
  }

bool PmAngularMotor::isActive() { return active; }

void PmAngularMotor::setActive(const bool flag) { active = flag; }

//*============================================================*
//*==========                  addAngleData          ==========*
//*============================================================*
// add angle data. 

void 
PmAngularMotor::addAngleData(const float ang) {
  data.push_back(ang);
  }

//*============================================================*
//*==========                  angle                 ==========*
//*============================================================*
// get/set motor angles.

void 
PmAngularMotor::getAngle(const int id, float& val) { 
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  val = angle[i-1];
  }

void 
PmAngularMotor::setAngle(const int id, const float val) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  angle[i-1] = val;
  }

void 
PmAngularMotor::getMaxAngle(const int id, float& val) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  val = max_angle[i-1];
  }

void 
PmAngularMotor::setMaxAngle(const int id, const float val) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  max_angle[i-1] = val;
  }

void
PmAngularMotor::getAngleInc(const int id, float& val) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  val = angle_inc[i-1];
  }

void
PmAngularMotor::setAngleInc(const int id, const float val) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  angle_inc[i-1] = val;
  }

//*============================================================*
//*==========                  axes                  ==========*
//*============================================================*
// get/set motor axes.

void 
PmAngularMotor::setAxes(const int num, int axes[3])
  {
  if (num > 3) {
    num_axes = 3;
    }
  else {
    num_axes = num;
    }

  for (int i = 0; i < num_axes; i++) {
    if (axes[i] <= 3) { 
      this->axes[i] = axes[i];
      }
    }
  }

void 
PmAngularMotor::getAxes(int& num, int axes[3]) {
  num = num_axes;
  for (int i = 0; i < num_axes; i++) {
    axes[i] = this->axes[i];
    }
  }


//*============================================================*
//*==========               MaxForce                 ==========*
//*============================================================*
// set/get maximum force parameters.

void 
PmAngularMotor::setMaxForce(const int id, const float force) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  max_force[i] = force;
  }

void
PmAngularMotor::getMaxForce(const int id, float& force) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  force = max_force[i];
  }

//*============================================================*
//*==========               MaxVelocity              ==========*
//*============================================================*
// set/get maximum velocity parameters.

void 
PmAngularMotor::setMaxVelocity(const int id, const float vel) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  max_velocity[i] = vel;
  }

void
PmAngularMotor::getMaxVelocity(const int id, float& vel) {
  int i = id;
  if (i < 1) i = 1;
  if (i > 3) i = 3;
  vel = max_velocity[i];
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a motor.

void 
PmAngularMotor::display(PmVector3& loc)
  {
  #ifdef dbg_PmAngularMotor_display
  fprintf (stderr, "\n>>>>>>> PmAngularMotor::display \n");
  #endif

  /*
  PmGraphicsLine *line;
  PmGraphicsPoint *point;
  PmGraphicsAttributes patts;

  if (!line_geometry) {
    string geom_name;
    geom_name = "motor[" + name + "]";
    PmVector3 *verts = new PmVector3[2];
    verts[0] = loc;
    verts[1] = loc + scale*direction;
    line = new PmGraphicsLine(geom_name, 2, verts);
    line->display();
    line_geometry = dynamic_cast<PmGraphicsGeometry*>(line);
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(line_geometry);
    }

  PmVector3 verts[2]; 
  verts[0] = loc;
  verts[1] = loc + scale*direction;
  line->update(2, verts);
  */
  }

//*============================================================*
//*==========               joint                    ==========*
//*============================================================*
// set/get joint for motor.

void PmAngularMotor::setJoint(PmJoint *joint) { this->joint = joint; }
void PmAngularMotor::getJoint(PmJoint **joint) { *joint = this->joint; }

}





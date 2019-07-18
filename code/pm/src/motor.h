
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

#ifndef _MOTOR_PM_H_
#define _MOTOR_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

class PmGraphicsGeometry;
class PmJoint;

enum PmMotorType {
  PM_MOTOR_UNKNOWN,
  PM_MOTOR_LINEAR,
  PM_MOTOR_ANGULAR
  };

// PmMotor 
// -------

class PM_EXPORT PmMotor {
   public:
      PmMotor();
      ~PmMotor(){};
      void getName(string& name) { name = this->name; }
      void setTime(PmTimeInterval& time) { this->time_interval = time; }
      PmMotorType getType() { return type; }
      static void convMotorType(const string str, PmMotorType& type);
      static PmMotor* create(const string name, const PmMotorType type);
      virtual bool isActive(float time) = 0;
      virtual bool isActive() = 0;
      virtual void setActive(const bool flag) = 0;
      bool isEnabled();
      void setEnabled(const bool flag);
      virtual void display(PmVector3& loc) = 0;

      void setWriteResultsParams(const string file_name, const bool binary);
      void writeResults(const vector<float>& time, const int last_written, const int num);

   protected:
      string name;
      PmMotorType type;
      PmTimeInterval time_interval;
      bool active, enabled;

      // data //
      vector<float> data;

      // results i/o //
      bool binary_format;
      string write_file_name;
      FILE *file_ptr;
      bool write_init;
      int last_written;
   };


// PmAngularMotor 
// ------------------------
// angular motor.

class PM_EXPORT PmAngularMotor : public PmMotor {

   public:
      PmAngularMotor(const string name);
      ~PmAngularMotor(){};
      void getAngle(const int id, float& angle);
      void setAngle(const int id, const float angle);
      void getAngleInc(const int id, float& angle);
      void setAngleInc(const int id, const float angle);
      void addAngleData(const float ang);
      void getMaxAngle(const int id, float& angle);
      void setMaxAngle(const int id, const float angle);
      void setAxes(const int num, int axes[3]);
      void getAxes(int& num, int axes[3]);
      bool isActive(float time);
      bool isActive();
      void setActive(const bool flag);
      void display(PmVector3& loc);
      void getMaxForce(const int id, float& force);
      void setMaxForce(const int id, float force);
      void setJoint(PmJoint *joint);
      void getJoint(PmJoint **joint);
      void getMaxVelocity(const int id, float& vel);
      void setMaxVelocity(const int id, const float vel);

   private:
      int num_axes, axes[3];
      float max_velocity[3], max_force[3]; 
      float max_angle[3], angle[3], angle_inc[3];
      PmJoint *joint;
      PmGraphicsGeometry *line_geometry;
   };

}

#endif


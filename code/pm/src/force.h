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
//* force:                    f o r c e                        *
//*============================================================*

#ifndef _FORCE_PM_H_
#define _FORCE_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

class PmGraphicsGeometry;

enum PmForceType {
  PM_FORCE_UNKNOWN,
  PM_FORCE_EXPLICIT,
  PM_FORCE_RANDOM
  };

typedef struct PmForceVector {
  PmVector3 point;
  PmVector3 direction;
  int index;
  } PmForceVector;


// PmForce
// -------------

class PM_EXPORT PmForce {
   public:
      PmForce();
      ~PmForce(){};
      void getName(string& name) { name = this->name; }
      void setScale(const float scale);
      void setTime(PmTimeInterval& time) { this->time_interval = time; }
      bool getRampParams(const float time, float& dramp);
      void setRampParams(const float param);
      PmForceType getType() { return type; }
      void setColor(PmVector3& color);
      void setShow(bool show);
      static void convForceType(const string str, PmForceType& type);
      static PmForce* create(const string name, const PmForceType type);
      static void procQuery(PmQuery& query);
      virtual bool isActive(float time) = 0;
      virtual void display(PmVector3& loc) = 0;
   protected:
      string name;
      PmForceType type;
      PmTimeInterval time_interval;
      bool has_ramp, ramp_init;
      float ramp_start, ramp_time, dramp;
      bool active;
      float scale;
      PmVector3 color;
      bool show;
   };


// PmExplicitForce 
// ------------------------
// an explicitly set force.

class PM_EXPORT PmExplicitForce : public PmForce {

   public:
      PmExplicitForce(const string name);
      ~PmExplicitForce(){};
      void getDirection(PmVector3& dir) { dir = this->direction; }
      void getDirection(const float time, PmVector3& dir);
      void setDirection(const PmVector3& dir) { this->direction = dir; }
      void setGlobalFrame(const bool flag) { this->global_frame = flag; }
      void getGlobalFrame(bool& flag) { flag = this->global_frame; }
      void getObj(void **obj) { *obj = this->sobj; }
      void setObj(void *obj) { this->sobj = obj; }
      void getPoint(PmVector3& point) { point = this->point; }
      void setPoint(const PmVector3& point) { this->point = point; }
      bool isActive(float time);
      void setTorque(const bool flag);
      bool isTorque();
      void display(PmVector3& loc);
      void setInteractive(const bool flag);
      void getInteractive(bool& flag);
      void setAxes(PmVector3& axis1, PmVector3& axis2, PmVector3& axis3);
      void setXform(PmXform& xform); 
      void setStrength(const float val);
      void procAxesSelect(PmQuery& query);
      void setUsePickingSphere(const bool use_sphere);
      void setPickingSphereRadius(const float radius);

   private:
      PmVector3 point, current_point;
      PmVector3 direction;
      float strength, picking_sphere_radius;
      PmXform xform;
      bool global_frame;
      bool torque;
      bool interactive, use_picking_sphere;
      vector<PmVector3> interactive_axes;
      int index;
      void *sobj;
      PmGraphicsGeometry *line_geometry, *point_geometry, *axes_line_geometry;
      PmGraphicsGeometry *axes_sphere_geometry, *picking_sphere_geometry;
   };


// PmRandomForce
// -------------
// random force.

class PM_EXPORT PmRandomForce : public PmForce {

   public:
      PmRandomForce(const string name);
      ~PmRandomForce(){};
      void getObj(void **obj) { *obj = this->sobj; }
      void setObj(void *obj) { this->sobj = obj; }
      bool isActive(float time);
      void display(PmVector3& loc);
      void getData(float& sd, float& mean, long int& seed);
      void genForce(PmVector3& rvec);
      void setMean(const float mean);
      void setSeed(const long int seed);
      void setStandardDev(const float sd);

   private:
      PmVector3 point;
      PmVector3 direction;
      bool global_frame;
      int index;
      long int seed, current_seed;
      float standard_deviation, mean;
      void *sobj;
      PmGraphicsGeometry *line_geometry, *point_geometry;
   };

}

#endif


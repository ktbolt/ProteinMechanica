
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

#ifndef _JOINT_PM_H_ 
#define _JOINT_PM_H_

#include "pm/pm.h"
#include "body.h"
#include "graphics.h"

namespace ProteinMechanica {

class PmMotor;

// PmJoint
// -------

 class PM_EXPORT PmJoint {

   public:
     PmJoint() {};
     PmJoint(const string name, const PmJointType type);
     static void convJointType(const string str, PmJointType& type);
     void getBodies (PmBody **body1, PmBody **body2) { 
                     *body1 = this->body1; 
                     *body2 = this->body2; }
     void setBodies (PmBody *body1, PmBody *body2) { 
                     this->body1 = body1; 
                     this->body2 = body2; }
     static PmJoint* create(const string name, const PmJointType type);
     void setColor (PmVector3& color) { this->color = color; }
     void setDisplayType(PmGeometryDisplayType type) { this->display_type = type; }
     void addEnergy(const bool upd, const float energy);
     void getEnergy(const int step, float& val);
     void addMotor(PmMotor *motor);
     void getMotors(vector<PmMotor*>& motors);
     bool hasMotors();
     void setMsize (const float msize) { this->msize = msize; }
     void getName (string& name) { name = this->name;  }
     bool hasName(const string& name) { return name == this->name;  }
     void getCurrentPosition (PmVector3& pos);
     void getPosition (PmVector3& pos);
     void setPosition (const PmVector3& pos);
     void setShadingType(PmGeometryShadingType type) { this->shading_type = type; }
     void setShowAxes(const bool val) { this->show_axes = val; }
     void setShow(const bool val) { this->show = val; }
     PmJointType getType() { return type; }
     void getType(string& str);

     virtual void getForceConsts(vector<float>& k) = 0;
     virtual void init() = 0;
     virtual void display(const bool show) = 0;
     virtual void update(const PmVector3& pos) = 0;
     virtual void update(PmXform& body_xform) = 0;

     void setWriteResultsParams(const string file_name, const bool binary);
     virtual void writeResults(const vector<float>& time, const int last_written, 
                               const int num)=0;

   protected:
     string name;
     PmJointType type;
     bool position_set;
     PmVector3 position, current_position;
     PmBody *body1, *body2;
     vector<PmMotor*> motors;
     vector<float> energy;

     // graphics //
     PmVector3 color;
     float msize;
     PmGeometryShadingType shading_type;
     PmGeometryDisplayType display_type;
     bool show, show_axes;

     // results i/o //
     string write_file_name;
     FILE *file_ptr;
     bool binary_format;
     bool write_init;
     int last_written;
   };


// PmBallJoint
// -----------

 class PM_EXPORT PmBallJoint : public PmJoint {

   public:
     PmBallJoint(const string name, const PmJointType type) : 
       PmJoint(name, type) { init(); };
     void setAxes(const PmVector3& a1, const PmVector3& a2);
     void getAxis(const int id, PmVector3& a);
     void getForceConst(const int id, float& k);
     void setForceConst(const int id, const float k);
     void getForceConsts(vector<float>& k);
     void getRange(const int id, float& r);
     bool hasRange(const int id);
     void setRange(const int id, const float r);
     void update(const PmVector3& pos);
     void update(PmXform& body_xform);
     void display(const bool show);
     void writeResults(const vector<float>& time, const int last_written, 
                       const int num){};

   private:
     PmVector3 axes[3];
     float force_const[3];
     float range[3];
     bool range_set[3];
     PmGraphicsBallJoint *geometry;
     void init();
   };


// PmFreeJoint
// -----------

 class PM_EXPORT PmFreeJoint : public PmJoint {

   public:
     PmFreeJoint(const string name, const PmJointType type) :
       PmJoint(name, type) { init(); };
     void display(const bool show);
     void getForceConsts(vector<float>& k);
     void update(const PmVector3& pos);
     void update(PmXform& body_xform) {};
     void writeResults(const vector<float>& time, const int last_written, 
                       const int num){};

   private:
     PmGraphicsFreeJoint *geometry;
     void init();
   };


// PmHingeJoint
// ----------

 class PM_EXPORT PmHingeJoint : public PmJoint {

   public:
     PmHingeJoint(const string name, const PmJointType type) :
       PmJoint(name, type) { init(); }

     void display(const bool show);
     void getAxis(PmVector3& vec);
     void setAxis(const PmVector3& vec);
     void getForceConst(float& k);
     void setForceConst(const float k);
     void getForceConsts(vector<float>& k);
     void setRange(const float r);
     bool getRange(float& r);
     void update(const PmVector3& pos);
     void update(PmXform& body_xform);
     void updateAngle(float a);
     void getAngle(float& a);
     void getAngles(vector<float>& a);
     void writeResults(const vector<float>& time, const int last_written, 
                       const int num);
     void setPosition2 (const PmVector3& pos);
     void getPosition2 (PmVector3& pos);
     void getCurrentPosition2 (PmVector3& pos);

   private:
     PmVector3 axis, position2, current_position2;
     bool position2_set;
     float force_const; 
     vector<float> angles; 
     float range;
     bool range_set;
     PmGraphicsHingeJoint *geometry;
     void init();
   };


// PmUniversalJoint
// ----------------

 class PM_EXPORT PmUniversalJoint : public PmJoint {

   public:
     PmUniversalJoint(const string name, const PmJointType type) :
       PmJoint(name, type) { init(); }

     void display(const bool show);
     void getAxis(const int id, PmVector3& a);
     void setAxis(const int id, const PmVector3& a);
     void getForceConst(const int id, float& k);
     void getForceConsts(vector<float>& k);
     void setForceConst(const int id, const float k);
     void update(const PmVector3& pos);
     void update(PmXform& body_xform);
     void writeResults(const vector<float>& time, const int last_written, 
                       const int num){};

   private:
     PmVector3 axis[2];
     float force_const[2];
     PmGraphicsUniversalJoint *geometry;
     void init();
   };

// PmWeldJoint
// -----------

 class PM_EXPORT PmWeldJoint : public PmJoint {

   public:
     PmWeldJoint(const string name, const PmJointType type) :
       PmJoint(name, type) { init(); };
     void display(const bool show);
     void getForceConsts(vector<float>& k);
     void update(const PmVector3& pos);
     void update(PmXform& body_xform);
     void writeResults(const vector<float>& time, const int last_written, 
                       const int num){};

   private:
     PmGraphicsWeldJoint *geometry;
     void init();
   };


}

#endif




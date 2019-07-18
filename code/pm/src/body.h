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
//* body:                   b o d y                            *
//*============================================================*
// rigid body class.                      

#ifndef _BODY_PM_H_ 
#define _BODY_PM_H_

#include "pm/pm.h"
#include "pobj.h"
#include "sobj.h"
#include "state.h"

namespace ProteinMechanica {

typedef enum PmBodyType {
  PM_BODY_UNKNOWN,
  PM_BODY_GROUND,
  PM_BODY_STATIC,
  PM_BODY_RIGID,
  PM_BODY_DEFORM,
  PM_BODY_SIZE
  } PmBodyType;

class PmJoint;

typedef struct PmRigidSimResults {
  PmVector3 disp, rot;
  PmVector3 vel, avel;
  PmMatrix3x3 rot_mat;
  } PmRigidSimResults;

class PM_EXPORT PmRigidBodyState {
  public:
     PmRigidBodyState(){};
     ~PmRigidBodyState(){};
     PmVector3 displacement;
     PmVector3 rotation;
     PmVector3 velocity;
     PmVector3 angular_velocity;
     PmMatrix3x3 rotation_matrix;
  };


// PmBody
// ------
// rigid body class.

 class PM_EXPORT PmBody : public PmSimulationObj {
    public:
       PmBody(const string name, PmBodyType type);
       ~PmBody();

       PmBodyType getType() { return type; }
       void getType(string& str);
       bool isGround();
       bool isStatic();
       static void checkContact(PmVector3& color, bool show);
       void checkContact(PmBody *bdy, float tol, PmVector3& color1, PmVector3& color2, 
                         bool show);
       void getCurrentCoordinats(vector<PmVector3>& cur_coords);
       void getDomainCoords(const string seq, vector<PmVector3>& coords);
       bool getDampingFactor(float& factor, float& length);
       void setDampingFactor(const float val);
       void addKineticEnergy(bool upd, float ke);
       bool getKineticEnergy(const int step, float& ke);
       void getMassProps(PmMassProperties& props);
       void getCenterOfMass(PmVector3& com);
       void getPosition(PmVector3& pos);
       static void convBodyType(const string str, PmBodyType& type);
       void getSimResult(const int step, PmRigidSimResults& result);
       void getSimResult(PmRigidSimResults& result);
       void setSimResults(const int step);
       void addSimResults(bool upd, float time, PmVector3& disp, PmVector3& rot, 
                          PmVector3& vel, PmVector3& avel, PmMatrix3x3& mat, 
                          PmQuaternion& q);
       void getState(PmRigidBodyState& state);
       void setState(const PmRigidBodyState& state); 
       void getVelocity(PmVector3& vel);
       void getVelocity(const PmVector3& point, PmVector3& vel);
       void getXform(PmXform& xform);

    private:
       PmBodyType type;
       vector <PmJoint*> joints;
       PmMassProperties mass_properties;
       vector<PmRigidSimResults> simulation_results;
       vector<float> kinetic_energy;
       PmXform xform;
       bool has_damping;
       float damping_factor, damping_length;
       PmRigidBodyState current_state, last_state;
    };
}

#endif




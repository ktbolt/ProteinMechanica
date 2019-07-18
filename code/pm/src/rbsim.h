
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
//* rbsim:    r i g i d   b o d y   s i m u l a t i o n        *
//*============================================================*

#ifndef _RIGID_SIMULATION_PM_H_ 
#define _RIGID_SIMULATION_PM_H_

#include "pm/pm.h"
#include "sim.h"
#include "body.h"
#include "joint.h"
#include "force.h"
#include "geom.h"
#include "motor.h"
#include "msr.h"
#include "rbsolv.h"
#include "rest.h"
#include "trace.h"

namespace ProteinMechanica {

 class PM_EXPORT PmRigidSimObj {
    public:
       PmRigidSimObj() { body = NULL; };
       string name;
       PmBody *body;
    };

// PmRigidSimulation 
// ---------------------
// rigid body simulation 

 class PM_EXPORT PmRigidSimulation : public PmSimulation {

   public:
      //PmRigidSimulation (const string name) : PmSimulation(name) { };
      PmRigidSimulation (const string name);

      bool addBody(PmBody *body);
      bool getBody(const string name, PmBody **body);
      void getBodies(vector<PmBody*>& blist);

      void setDamping(const bool flag);
      bool getDamping();

      void initialize();

      bool addForce(const string name, PmForce *force, PmTimeInterval& time); 
      void displayForce(PmExplicitForce *force, PmVector3& loc);
      void getForces(vector<PmExplicitForce*>& forces);
      void getForces(vector<PmRandomForce*>& forces);

      bool addGeometry(PmGeometry *geom);
      void getGeometries(vector<PmGeometry*>& geoms);

      void halt();

      bool addJoint(PmJoint *joint);
      bool getJoint(const string name, PmJoint **joint);
      void getJoints(vector<PmJoint*>& jlist);

      bool addMeasurement(PmMeasurement *msr);
      void getMeasurement(const string name, PmMeasurement **msr);
      void getMeasurements(vector<PmMeasurement*>& msr);

      void printModelInformation(const bool pbodies, const bool pjoints); 

      void setMomentumOff(const bool flag);
      bool getMomentumOff();

      bool addMotor(PmMotor *motor, PmTimeInterval& time); 
      void getMotor(const string name, PmMotor **motor);

      void replay();

      void displayRestraint(const string name, PmVector3& pt1, PmVector3& pt2, 
                            const bool show);

      void setSolver(PmRigidBodySolver *solver);
      bool setSolverParam(const string name, const string val);
      void step(int nsteps);
      void readState(const string file_name);

      void getTime(double& time);
      void setTimeStep(const double step);
      void getTimeStep(double& step);

      bool addTrace(PmTrace *trace);
      void getTraces(vector<PmTrace*>& traces);
      void getTrace(const string name, PmTrace **trace);

      void writeResults();

   private:
     bool momentum_off, damping;
     double current_time, time_step;
     vector <PmBody*> bodies;
     vector <PmJoint*> joints;
     vector <PmPhysicalObj*> physical_objects;
  
     PmRigidBodySolver *solver;

     vector <PmRigidSimObj> simulation_objects;

     bool getBodyId (const string str, int& id);

     vector <PmExplicitForce*> explicit_forces;
     vector <PmRandomForce*> random_forces;
     vector <PmMeasurement*> measurements;
     vector <PmMotor*> motors;
     vector <PmTrace*> traces;
     vector<PmGeometry*> geometries;

     PmRigidSimResults simulation_results;

     vector<PmGraphicsGeometry*> graphics_geometry;
     void getGraphicsGeometry(const string name, PmGraphicsGeometry **geom);
     void getRestraintGeometry (const string name, PmGraphicsRestraint **rgeom);
     void updateEnergy();
     bool incrementTimeStep();
   };

}

#endif





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
//* rbsolv:      r i g i d    b o d y   s o l v e r            *
//*============================================================*

#ifndef _RIGID_BODY_SOLVER_PM_H_ 
#define _RIGID_BODY_SOLVER_PM_H_

#include "pm/pm.h"
//#include "rbsim.h"
#include "body.h"

namespace ProteinMechanica {

class PmRigidSimulation;

typedef enum {
  PM_RIGID_BODY_SOLVER_UNKNOWN,
  PM_RIGID_BODY_SOLVER_ODE,
  PM_RIGID_BODY_SOLVER_SDFAST,
  PM_RIGID_BODY_SOLVER_SIMBODY,
  PM_RIGID_BODY_SOLVER_SIZE
  } PmRigidBodySolverType;


// PmRigidBodySolver
// -----------------

 class PM_EXPORT PmRigidBodySolver {
    public:
      PmRigidBodySolver(){};
      PmRigidBodySolver(const string name, PmRigidBodySolverType type);
      static void convSolverType(const string str, PmRigidBodySolverType& type);
      static PmRigidBodySolver* create(const string name, 
                                       const PmRigidBodySolverType type);
      virtual void initialize (PmRigidSimulation *sim, double time_step) = 0;
      virtual bool setParam(const string name, const string val) = 0;
      virtual void step() = 0;
      virtual void setState(vector<string>& body_names, 
                            vector<PmRigidBodyState>& states) = 0;
      virtual void setTimeStep(const double step) = 0;

    protected:
      string name;
      bool initialzed;
      PmRigidBodySolverType type;
      PmRigidSimulation *simulation;
   };


// PmOdeSolver      
// ---------------------------
// open dynamics engine solver

class PM_EXPORT PmOdeSolver : public PmRigidBodySolver {
   public:
      PmOdeSolver(const string name);
      void initialize (PmRigidSimulation *sim, double time_step);
      bool setParam(const string name, const string val);
      void step();
      void setTimeStep(const double step);
      void setState(vector<string>& body_names, vector<PmRigidBodyState>& states);

   private:
      void *prv_data;
      void addBodies(vector<PmBody*> bodies);
      bool checkBodies(vector<PmJoint*> joints, vector<PmBody*> bodies);

      void stepMotors();

      void addForces();
      void addDampingForces();
      void addExplicitForces();
      void addRandomForces();
      void addInteractionForces();
      void addJointForces();
      void addJointMotors();
      void addRestraintForces();

      bool addJoints(vector<PmJoint*> joints);
      void update();
   };
   
}

#endif




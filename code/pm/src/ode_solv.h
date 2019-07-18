
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
//* ode_solv:             O D E  s o l v e r                   *
//*============================================================*
// interface to the Open Dynamics Engine.

#ifndef _ODE_SOLVER_PM_H_ 
#define _ODE_SOLVER_PM_H_

#include "pm/pm.h"
#include "rbsim.h"
#include "rbsolv.h"
#include "rest.h"

#define PM_ODE_SOLVER

#ifdef PM_ODE_SOLVER
#include "ode/ode.h"
#endif

namespace ProteinMechanica {


#ifdef PM_ODE_SOLVER

typedef struct PmOdeBody {
  dBodyID ode;
  PmBody *pm;
  } PmOdeBody;

typedef struct PmOdeJoint {
  dJointID motors[3];
  dJointID ode;
  PmJoint *pm;
  } PmOdeJoint;

typedef struct PmOdeSolverData {
  int step;
  double time, time_step;
  int num_joints;
  vector<PmOdeJoint> joints;
  int num_bodies;
  vector<PmOdeBody> bodies;
  dWorldID world;
  bool quick_step;
  } PmOdeSolverData;

static void
ode_getBody(PmOdeSolverData *prvd, PmBody *body, dBodyID *obody);

static void
ode_getJoint(PmOdeSolverData *prvd, PmJoint *joint, dJointID *ojoint);

static void
ode_addJointStop(PmOdeSolverData *prv_data, PmOdeJoint& ojoint);

static void
ode_activateJointMotor(PmOdeSolverData *prv_data, PmOdeJoint& ojoint, bool enable);

static void
ode_xformPoint(PmBody *body, dBodyID obody, PmVector3& pt, PmVector3& xpt);

#endif
   

}

#endif




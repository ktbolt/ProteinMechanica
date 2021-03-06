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
//* ode_solv:             O D E   s o l v e r                  *
//*============================================================*
// interface to the OpenDynamicsEngine.

#include "ode_solv.h"
#include "pm/mth.h"

// debug symbols //

//#define dbg_PmOdeSolver

namespace ProteinMechanica {

static int vel_param[] = { dParamVel, dParamVel2, dParamVel3 };
static int force_param[] = { dParamFMax, dParamFMax2, dParamFMax3 };


//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmOdeSolver::PmOdeSolver(const string name)
  {
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, ">>>>>> PmOdeSolver: ctor  [%s] \n", name.c_str());
  #endif
  this->name = name;
  initialzed = false;
  }

//*============================================================*
//*==========             initialize                 ==========*
//*============================================================*

void
PmOdeSolver::initialize(PmRigidSimulation *sim, double time_step)
  {
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "\n-------------------------------\n");
  fprintf (stderr, ">>>>>> ODE Solver: initialize \n");
  fprintf (stderr, "   >>> size of dReal[%d]\n", sizeof(dReal));
  #endif
  simulation = sim;

#ifdef PM_ODE_SOLVER

  if (initialzed) {
    return;
    }

  dInitODE();

  // get bodies from simulation //
  vector<PmBody*> bodies;
  sim->getBodies(bodies);
  int num_bodies = bodies.size();
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "   >>> num bodies [%d] \n", num_bodies);
  #endif

  if (num_bodies == 0) {
    return;
    }

  PmOdeSolverData *prvd = new PmOdeSolverData;
  prv_data = prvd;

  // get joints //

  vector<PmJoint*> joints;
  sim->getJoints(joints);
  int num_joints = joints.size();
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "   >>> num joints [%d] \n", num_joints);
  #endif

  // check to make sure the bodies defined for the joints //
  // have been added to the simualtion.                   //

  if (!checkBodies(joints, bodies)) {
    return;
    }

  // create ode world //

  prvd->world = dWorldCreate();
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "   >>> size of ODE real [%d] \n", sizeof(dReal));
  #endif

  // set error control //

  //dReal erp, cfm;
  //dWorldSetERP (prvd->world, erp);
  //dWorldSetCFM (prvd->world, cfm);

  // set quick step option //

  prvd->quick_step = false; 

  // add bodies //

  addBodies(bodies);

  // add joints //

  if (!addJoints(joints)) {
    return;
    }

  prvd->num_bodies = num_bodies;
  prvd->num_joints = num_joints;
  prvd->time = 0.0;
  prvd->time_step = time_step;
  prvd->step = 0;
  initialzed = true;

  // add initial state //

  update();

#else

#endif
  }

//*============================================================*
//*==========             setParam                   ==========*
//*============================================================*
// set an ode solver parameter.

bool 
PmOdeSolver::setParam(const string name, const string val)
  {
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "\n-------------------------------\n");
  fprintf (stderr, ">>>>>> ODE Solver: set parameter \n");
  #endif
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;

  if (name == "finite_rotation") {
    return true;
    }

  else if (name == "momentum") {
    return true;
    }

  else if (name == "quick_step") {
    prvd->quick_step = (val == "true");
    return true;
    }

  return false;
  }

//*============================================================*
//*==========                 step                   ==========*
//*============================================================*
// step the ode solver.

void 
PmOdeSolver::step()
  {
#ifdef PM_ODE_SOLVER

  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "\n-------------------------------\n");
  fprintf (stderr, ">>>>>> ODE Solver: step \n");
  #endif

  if (!initialzed) {
    return;
    }

  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "       time step [%g] \n", prvd->time_step);
  #endif

  // step motors //
  stepMotors();

  // add forces to bodies //
  addForces();

  // step dynamics // 

  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "       step dynamics \n" );
  #endif

  if (prvd->quick_step) { 
    dWorldQuickStep (prvd->world, prvd->time_step);
    }
  else {
    dWorldStep (prvd->world, prvd->time_step);
    }

  prvd->step += 1;
  prvd->time += prvd->time_step;

  //  if momentum is off then zero velocities  //

  if (simulation->getMomentumOff()) { 
    //fprintf (stderr, ">>> momentum off \n");

    for (unsigned int i = 0; i < prvd->bodies.size(); i++) {
      dBodyID obody = prvd->bodies[i].ode;

      if (obody) { 
        dBodySetLinearVel (obody, 0.0, 0.0, 0.0);
        dBodySetAngularVel (obody, 0.0, 0.0, 0.0);
        }
      }
    }

  // update objects state //

  update();

#endif
  }

//*============================================================*
//*==========                 setState               ==========*
//*============================================================*
// set the state for the bodies in the simulation. 

void 
PmOdeSolver::setState(vector<string>& body_names, vector<PmRigidBodyState>& states)
  {
  //fprintf (stderr, ">>>>>> PmOdeSolver::setState \n");
  PmBody *body;
  dBodyID obody;
  const dReal *pos;
  dMatrix3 R;
  PmVector3 disp, rot, new_pos, vel, avel;
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;

  for (int i = 0; i < body_names.size(); i++) {
    //fprintf (stderr, "--- body %s --- \n", body_names[i].c_str());
    pmSystem.getBody(body_names[i], &body);

    if (!body) {
      continue;
      }

    ode_getBody(prvd, body, &obody);
    pos = dBodyGetPosition (obody);
    //fprintf (stderr, ">>> pos = %f %f %f \n", pos[0], pos[1], pos[2]);

    disp = states[i].displacement;

    // set position //

    for (int j = 0; j < 3; j++) {
      new_pos[j] = pos[j] + disp[j];
      }

    dBodySetPosition (obody, new_pos[0], new_pos[1], new_pos[2]);

    // set rotation //

    R[0] = states[i].rotation_matrix(0,0);
    R[1] = states[i].rotation_matrix(0,1);
    R[2] = states[i].rotation_matrix(0,2);
    R[3] = 0.0; 

    R[4] = states[i].rotation_matrix(1,0);
    R[5] = states[i].rotation_matrix(1,1);
    R[6] = states[i].rotation_matrix(1,2);
    R[7] = 0.0; 

    R[8] = states[i].rotation_matrix(2,0);
    R[9] = states[i].rotation_matrix(2,1);
    R[10] = states[i].rotation_matrix(2,2);
    R[11] = 0.0;
 
    dBodySetRotation (obody, R);

    // set velocity //

    /*
    fprintf (stderr, ">> vel = %f %f %f \n", states[i].velocity[0], 
             states[i].velocity[1],  states[i].velocity[2]);
    */

    dBodySetLinearVel (obody, states[i].velocity[0], states[i].velocity[1],  
                       states[i].velocity[2]);

    dBodySetAngularVel (obody, states[i].angular_velocity[0],  
                        states[i].angular_velocity[1], states[i].angular_velocity[2]);
    }

  update();
  }

//*============================================================*
//*==========            setTimeStep                 ==========*
//*============================================================*
// set time step.

void
PmOdeSolver::setTimeStep(const double step)
  {
#ifdef PM_ODE_SOLVER
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  prvd->time_step = step;
#endif
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             addBodies                  ==========*
//*============================================================*
// add bodies to ode solver.

void
PmOdeSolver::addBodies(vector<PmBody*> bodies)
  {
#ifdef PM_ODE_SOLVER
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  int num_bodies = bodies.size();
  PmOdeBody obody;
  PmBody *body;
  PmMassProperties props;
  dMass omass_props;

  PmBodyType type;
  string name;
  dReal xp, yp, zp;
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "   --- add bodies --- \n");
  #endif

  for (int i = 0; i < num_bodies; i++) {
    body = bodies[i];
    body->getName(name);
    //body->setRotOrderZYX(true);
    type = body->getType();
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "   >>> body[%s]  type[%d] \n", name.c_str(), type);
    #endif

    if (type != PM_BODY_GROUND) {

      // create an ode body //

      obody.ode = dBodyCreate (prvd->world);

      // set mass properties //

      body->getMassProps(props);
      xp = props.com[0];
      yp = props.com[1];
      zp = props.com[2];
      dBodySetPosition (obody.ode, xp, yp, zp);


      #ifdef dbg_PmOdeSolver
      fprintf (stderr, "       mass[%g] \n", props.mass);
      fprintf (stderr, "       com (%g %g %g) \n", xp, yp, zp); 
      fprintf (stderr, "       I   (%g %g %g) \n", props.inertia(0,0),
                                                   props.inertia(1,1), 
                                                   props.inertia(2,2));
      fprintf (stderr, "           (%g %g %g) \n", props.inertia(0,1),
                                                   props.inertia(0,2), 
                                                   props.inertia(1,2));
      #endif


      // center of mass in body frame //

      props.com[0] = 0.0;
      props.com[1] = 0.0;
      props.com[2] = 0.0;

      dMassSetParameters (&omass_props, props.mass,
                          props.com[0], props.com[1], props.com[2],
                          props.inertia(0,0), props.inertia(1,1),
                          props.inertia(2,2), props.inertia(0,1),
                          props.inertia(0,2), props.inertia(1,2));

      dBodySetMass (obody.ode, &omass_props);
      dBodySetFiniteRotationMode (obody.ode, 1);
      dBodySetFiniteRotationAxis (obody.ode, 0.0, 0.0, 0.0);
      }
    else {
      obody.ode = 0;
      }
 
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "       map pmbody[%x]  to ode body[%x] \n", body, obody.ode); 
    #endif
    obody.pm = body;
    prvd->bodies.push_back(obody);
    }
#endif
  }

//*============================================================*
//*==========             checkBodies                ==========*
//*============================================================*
// check bodies.                     

bool
PmOdeSolver::checkBodies(vector<PmJoint*> joints, vector<PmBody*> bodies)
  {
#ifdef PM_ODE_SOLVER
  //PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  int num_bodies = bodies.size();
  int num_joints = joints.size();
  PmBody *body1, *body2;

  // check to make sure the bodies defined for the joints   //
  // have been added to the simualtion.                     //

  for (int i = 0; i < num_joints; i++) {
    PmJoint *joint = joints[i];
    joint->getBodies(&body1, &body2);

    if (!body1 || !body2) {
      string name;
      joint->getName(name);
      pm_ErrorReport ("pm> ", "joint \"%s\" bodies not defined.", "*", name.c_str());
      return false;
      }

    bool body_found = false;
    for (int j = 0; j < num_bodies; j++) {
      if (body1 == bodies[j]) {
        body_found = true;
        break;
        }
      }

    if (!body_found) {
      string name;
      body1->getName(name);
      pm_ErrorReport ("PM", "body \"%s\" not defined.", "*", name.c_str());
      return false;
      }

    body_found = false;
    for (int j = 0; j < num_bodies; j++) {
      if (body2 == bodies[j]) {
        body_found = true;
        break;
        }
      }

    if (!body_found) {
      string name;
      body2->getName(name);
      pm_ErrorReport ("PM", "body \"%s\" not defined.", "*", name.c_str());
      return false;
      }
    }

  return true;
#endif
  }

//*============================================================*
//*==========                 addForces              ==========*
//*============================================================*
// add forces to the bodies. 

//#define dbg_PmOdeSolver

void
PmOdeSolver::addForces()
  {
#ifdef PM_ODE_SOLVER
  #ifdef dbg_PmOdeSolver
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  fprintf (stderr, ">>> ODE Solver: addForces  time [%g] \n", prvd->time);
  #endif

  if (!initialzed) {
    return;
    }

  // reset forces //

  vector<PmBody*> bodies;
  simulation->getBodies(bodies);

  for (unsigned int i = 0; i < bodies.size(); i++) {
    bodies[i]->resetForces();
    }

  // add joint forces //

  addJointForces();

  // add restraint forces //

  addRestraintForces();

  // add explicit forces //

  addExplicitForces();

  // add random forces //

  addRandomForces();

  // add damping forces //

  addDampingForces();

  // add interaction forces //

  addInteractionForces();
  }

//*============================================================*
//*==========           addExplcitForces             ==========*
//*============================================================*
// add explicit forces. 

void
PmOdeSolver::addExplicitForces()
  {
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, ">>>>>> PmOdeSolver::addExplcitForces \n");
  #endif
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  vector<PmExplicitForce*> forces;
  simulation->getForces(forces);
  int num = forces.size();
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "    num forces[%d] \n", num);
  #endif

  if (!num) {
    return;
    }

  dBodyID obody;
  PmVector3 com, force, loc, vec, floc;
  bool global_frame = false;

  // ODE types //

  const dReal *pos, *r;
  dVector3 gvec;
  double global_pt[3]; 
  PmBody *body;
  float time = prvd->time;

  // apply each force to the appropriate body //

  for (int i = 0; i < num; i++) {
    if (!forces[i]->isActive(time)) {
      continue;
      }

    forces[i]->getGlobalFrame(global_frame);
    forces[i]->getObj((void**)&body);
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "--- %d  body[%x] ---  \n", i, body);
    #endif

    ode_getBody(prvd, body, &obody);
    body->getCenterOfMass(com);
    forces[i]->getDirection(time, force);

    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "com (%g %g %g) \n", com[0], com[1], com[2]);
    fprintf (stderr, "force (%g %g %g) \n", force[0], force[1], force[2]);
    #endif

    // apply a torque //

    if (forces[i]->isTorque()) {
      dBodyAddTorque(obody, force[0], force[1], force[2]);
      dBodyGetRelPointPos (obody, 0, 0, 0, gvec);

      for (int j = 0; j < 3; j++) {
        floc[j] = gvec[j];
        }

      forces[i]->display(floc);
      }

    // apply force in global frame //

    else if (global_frame) {
      forces[i]->getPoint(loc);
      pos = dBodyGetPosition (obody);
      r = dBodyGetRotation (obody);
      vec = loc - com;
      dBodyVectorToWorld (obody, vec[0], vec[1], vec[2], gvec);
      #ifdef dbg_PmOdeSolver
      fprintf (stderr, "pos (%g %g %g) \n", pos[0], pos[1], pos[2]);
      #endif

      for (int j = 0; j < 3; j++) {
        global_pt[j] = pos[j] + gvec[j];
        floc[j] = pos[j] + gvec[j];
        }

      dBodyAddForceAtPos (obody, force[0], force[1], force[2], global_pt[0], 
                          global_pt[1], global_pt[2]);

      /*
      fprintf (stderr, "    global_pt (%g %g %g) \n", global_pt[0], global_pt[1], 
               global_pt[2]);
      */

      //simulation->displayForce(forces[i], floc);

      PmMatrix3x3 m;
      PmXform body_xform;
      m(0,0) = r[0]; m(1,0) = r[4]; m(2,0) = r[8];
      m(0,1) = r[1]; m(1,1) = r[5]; m(2,1) = r[9];
      m(0,2) = r[2]; m(1,2) = r[6]; m(2,2) = r[10];

      body_xform.translation[0] = pos[0] - com[0];
      body_xform.translation[1] = pos[1] - com[1];
      body_xform.translation[2] = pos[2] - com[2];

      body->getCenterOfMass (com);
      body_xform.center = com;
      body_xform.setRotation(m);

      forces[i]->setXform(body_xform);
      forces[i]->display(floc);
      }

    // apply force in local frame //

    else {
      forces[i]->getPoint(loc);
      dBodyAddForceAtRelPos (obody, force[0], force[1], force[2], loc[0], loc[1], loc[2]);
      dBodyGetRelPointPos (obody, loc[0], loc[1], loc[2], gvec);

      for (int j = 0; j < 3; j++) {
        floc[j] = gvec[j];
        }

      forces[i]->display(floc);
      }
    }
#endif
  }

//*============================================================*
//*==========           addRandomForces              ==========*
//*============================================================*
// add random forces.

void
PmOdeSolver::addRandomForces()
  {
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, ">>>>>> PmOdeSolver::addRandomForces \n");
  #endif
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  vector<PmRandomForce*> forces;
  simulation->getForces(forces);
  int num = forces.size();
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "    num forces[%d] \n", num);
  #endif

  if (!num) {
    return;
    }

  dBodyID obody;
  PmVector3 com, force, loc, vec, floc;

  // ODE types //

  dVector3 gvec;
  PmBody *body;
  float time = prvd->time;

  // apply each random force to the appropriate body //

  for (int i = 0; i < num; i++) {
    if (!forces[i]->isActive(time)) {
      continue;
      }

    forces[i]->getObj((void**)&body);
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "--- %d  body[%x] ---  \n", i, body);
    #endif
    ode_getBody(prvd, body, &obody);
    body->getCenterOfMass(com);

    // generate and apply a random force //

    forces[i]->genForce(force);
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "com (%g %g %g) \n", com[0], com[1], com[2]);
    fprintf (stderr, "force (%g %g %g) \n", force[0], force[1], force[2]);
    #endif

    loc.set(0,0,0); 
    dBodyAddForceAtRelPos (obody, force[0], force[1], force[2], loc[0], loc[1], loc[2]);
    dBodyGetRelPointPos (obody, loc[0], loc[1], loc[2], gvec);

    for (int j = 0; j < 3; j++) {
      floc[j] = gvec[j];
      }

    forces[i]->display(floc);

    // generate and apply a random torque //

    forces[i]->genForce(force);
    dBodyAddTorque(obody, force[0], force[1], force[2]);
    }
  }

//*============================================================*
//*==========           addDampingForces             ==========*
//*============================================================*
// add damping forces.

void
PmOdeSolver::addDampingForces()
  {
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, ">>>>>> PmOdeSolver::addDampingForces\n");
  #endif
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;

  if (!simulation->getDamping()) { 
    return;
    }

  const dReal *vel, *avel;
  PmBody *body;
  dBodyID obody;
  PmVector3 com, force; 
  float c, factor, length;

  vector<PmBody*> bodies;
  simulation->getBodies(bodies);

  for (unsigned int i = 0; i < bodies.size(); i++) {
    body = bodies[i];

    if (body->isGround()) {
      continue;
      }

    if (!body->getDampingFactor(factor, length)) {
      continue;
      }

    // generate and apply a damping force //

    ode_getBody(prvd, body, &obody);
    vel = dBodyGetLinearVel(obody);
    c = factor;
    force[0] = -c*vel[0];
    force[1] = -c*vel[1];
    force[2] = -c*vel[2];
    dBodyAddForce (obody, force[0], force[1], force[2]);
    //fprintf (stderr, "   >>> force = %g %g %g \n", force[0], force[1], force[2]);

    // generate and apply a damping torque //

    avel = dBodyGetAngularVel(obody);
    length /= 2.0;

    if (length > 0.1) {
      c = factor*length*length;
      force[0] = -c*avel[0];
      force[1] = -c*avel[1];
      force[2] = -c*avel[2];
      dBodyAddTorque(obody, force[0], force[1], force[2]);
      }
    }
  }

//*============================================================*
//*==========           addInteractionForces         ==========*
//*============================================================*
// add interaction forces.

void
PmOdeSolver::addInteractionForces()
  {
  #ifdef dbg_PmOdeSolver_addInteractionForces
  fprintf (stderr, "\n------ addInteractionForces ------\n");
  #endif
  vector<PmInteraction*> interactions;
  string intr_name;
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  float time = prvd->time;

  //===== get the interactions defined for a simulation =====//

  simulation->getInteractions(interactions);
  bool upd = simulation->checkUpdateState();

  //===== compute interaction forces =====//

  for (unsigned int i = 0; i < interactions.size(); i++) {
    //fprintf (stderr, "--- interaction[%s] ---\n", intr_name.c_str());
    //if (interactions[i]->isActive(time)) { 
    interactions[i]->compForces(upd, time);
    }

  //===== apply forces generated from interactions =====//

  vector<PmBody*> bodies;
  vector<PmForceVector> forces;
  string bname;
  PmVector3 dir, loc;
  dBodyID obody;
  simulation->getBodies(bodies);

  for (unsigned int i = 0; i < bodies.size(); i++) {
    bodies[i]->getForces(forces);
    ode_getBody(prvd, bodies[i], &obody);
    //bodies[i]->getName(bname);
    //fprintf (stderr, ">>> body[%s]  num forces[%d]\n", bname.c_str(), forces.size());

    for (unsigned int j = 0; j < forces.size(); j++) {
      loc = forces[j].point;
      dir = forces[j].direction;
      dBodyAddForceAtPos (obody, dir[0], dir[1], dir[2], loc[0], loc[1], loc[2]);
      #ifdef dbg_PmOdeSolver_addInteractionForces
      fprintf (stderr, "-------------------\n");
      fprintf (stderr, ">>> loc (%g %g %g) \n", loc[0], loc[1], loc[2]); 
      fprintf (stderr, ">>> dir (%g %g %g) \n", dir[0], dir[1], dir[2]); 
      #endif
      }
    }

  // now apply any torques generated from the interactions //

  for (unsigned int i = 0; i < bodies.size(); i++) {
    bodies[i]->getTorques(forces);
    ode_getBody(prvd, bodies[i], &obody);

    for (unsigned int j = 0; j < forces.size(); j++) {
      dir = forces[j].direction;
      dBodyAddTorque(obody, dir[0], dir[1], dir[2]); 
      //fprintf (stderr, ">>> add torque (%g %g %g) \n", dir[0], dir[1], dir[2]); 
      }
    }
  }

//*============================================================*
//*==========           addJointForces               ==========*
//*============================================================*
// add joint forces.

void
PmOdeSolver::addJointForces()
  {
#ifdef PM_ODE_SOLVER
  #ifdef dbg_PmOdeSolver_addInteractionForces
  fprintf (stderr, "\n------ addJointForces ------\n");
  #endif

  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  PmBody *body1, *body2;
  PmJoint *joint;
  PmJointType jtype;
  dJointID ojoint;
  float energy;
  bool upd = simulation->checkUpdateState();

  // for each joint type add torque if it has a non-zero force constant //

  int num_joints = prvd->joints.size();

  for (int i = 0; i < num_joints; i++) {
    joint = prvd->joints[i].pm;

    if (joint->hasMotors()) {
      continue;
      }

    ojoint = prvd->joints[i].ode;
    jtype = joint->getType();
    joint->getBodies(&body1, &body2);

    // universal joint //

    if (jtype == PM_JOINT_UNIVERSAL) {
      PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);
      float k1, k2, ang1, ang2, dang, torque1, torque2;
      ujoint->getForceConst(1,k1);
      ujoint->getForceConst(2,k2);
      //fprintf (stderr, "   >>> ujoint  k1[%f]  k2[%f] \n", k1, k2);

      if (k1 + k2 == 0.0) {
        continue;
        }

      ang1 = dJointGetUniversalAngle1(ojoint);
      ang2 = dJointGetUniversalAngle2(ojoint);
      //fprintf (stderr, "   >>> ang1[%f]  ang2[%f] \n", ang1, ang2); 

      dang = ang1;
      if (dang > M_PI)  dang -= 2*M_PI;
      if (dang < -M_PI) dang += 2*M_PI;
      torque1 = dang*k1;

      dang = ang2;
      if (dang > M_PI)  dang -= 2*M_PI;
      if (dang < -M_PI) dang += 2*M_PI;
      torque2 = dang*k2;
      dJointAddUniversalTorques (ojoint, -torque1, -torque2);
      }

    // ball joint //

    else if (jtype == PM_JOINT_BALL) {
      PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);
      float k1, k2, k3, ang1, ang2, ang3, dang, torque1, torque2, torque3;
      PmVector3 axis1, axis2, axis3, u1, v1, w1, u2, v2, w2, dir, u, v, w;
      dBodyID obody1, obody2, obody;
      const dReal *r1, *r2;  
      PmMatrix3x3 m1, m2, rmat; 
      float dcm[3][3], rang[3]; 

      bjoint->getForceConst(1,k1);
      bjoint->getAxis(1,axis1);
      bjoint->getForceConst(2,k2);
      bjoint->getAxis(2,axis2);
      bjoint->getForceConst(3,k3);
      bjoint->getAxis(3,axis3);

      /*
      fprintf (stderr, "\n");
      fprintf (stderr, ">>> axis1 = %f %f %f \n", axis1[0], axis1[1], axis1[2]); 
      fprintf (stderr, ">>> axis2 = %f %f %f \n", axis2[0], axis2[1], axis2[2]); 
      fprintf (stderr, ">>> axis3 = %f %f %f \n", axis3[0], axis3[1], axis3[2]); 
      fprintf (stderr, ">>> k1 = %f   k2 = %f  k3 = %f \n", k1, k2, k3); 
      */

      if (k1 + k2 + k3 == 0.0) {
        continue;
        }

      ode_getBody(prvd, body1, &obody1);
      ode_getBody(prvd, body2, &obody2);

      // obody1 is 0 then it is the ground //

      if (obody1) { 
        r1 = dBodyGetRotation(obody1);
        m1(0,0) = r1[0]; m1(1,0) = r1[4]; m1(2,0) = r1[8];
        m1(0,1) = r1[1]; m1(1,1) = r1[5]; m1(2,1) = r1[9];
        m1(0,2) = r1[2]; m1(1,2) = r1[6]; m1(2,2) = r1[10];
        u1 = m1*axis1;
        v1 = m1*axis2;
        w1 = m1*axis3;
        obody = obody1;
        }
      else {
        m1.diag(1.0);
        u1 = axis1;
        v1 = axis2;
        w1 = axis3;
        }

      if (obody2) { 
        r2 = dBodyGetRotation(obody2);
        m2(0,0) = r2[0]; m2(1,0) = r2[4]; m2(2,0) = r2[8];
        m2(0,1) = r2[1]; m2(1,1) = r2[5]; m2(2,1) = r2[9];
        m2(0,2) = r2[2]; m2(1,2) = r2[6]; m2(2,2) = r2[10];
        u2 = m2*axis1;
        v2 = m2*axis2;
        w2 = m2*axis3;
        obody = obody2;
        }
      else {
        m2.diag(1.0);
        u2 = axis1;
        v2 = axis2;
        w2 = axis3;
        }

      /*
      fprintf (stderr, ">>> u1 = %f %f %f \n", u1[0], u1[1], u1[2]); 
      fprintf (stderr, "    v1 = %f %f %f \n", v1[0], v1[1], v1[2]); 
      fprintf (stderr, "    w1 = %f %f %f \n", w1[0], w1[1], w1[2]); 
      fprintf (stderr, "\n");
      fprintf (stderr, ">>> u2 = %f %f %f \n", u2[0], u2[1], u2[2]); 
      fprintf (stderr, "    v2 = %f %f %f \n", v2[0], v2[1], v2[2]); 
      fprintf (stderr, "    w2 = %f %f %f \n", w2[0], w2[1], w2[2]); 
      */

      // compute rotation from one axes to the other //

      pm_MathFrameRotation(u1, v1, w1, u2, v2, w2, rmat);

      /*
      u = rmat*u1 - u2;
      v = rmat*v1 - v2;
      w = rmat*w1 - w2;
      fprintf (stderr, ">>> diff vec = %f %f %f\n", u.length(), v.length(), w.length());  
      */

      // extract angles from rotation matrix //
      
      dcm[0][0] = rmat(0,0);  dcm[0][1] = rmat(0,1);  dcm[0][2] = rmat(0,2);
      dcm[1][0] = rmat(1,0);  dcm[1][1] = rmat(1,1);  dcm[1][2] = rmat(1,2);
      dcm[2][0] = rmat(2,0);  dcm[2][1] = rmat(2,1);  dcm[2][2] = rmat(2,2);
   
      pm_MathExtractRotations(dcm, rang);
      //fprintf (stderr, ">>> rang = %f %f %f \n", rang[0], rang[1], rang[2]); 

      // axis 1 //

      if (k1 != 0.0) {
        ang1 = rang[0];
        dang = ang1;
        if (dang > M_PI)  dang -= 2*M_PI;
        if (dang < -M_PI) dang += 2*M_PI;
        torque1 = dang*k1;
        dir = axis1;
        dir.set(1,0,0); 

        if (obody1) {
          //dir = u1;
          //dBodyAddRelTorque(obody1, torque1*dir[0], torque1*dir[1], torque1*dir[2]);
          dBodyAddTorque(obody1, torque1*dir[0], torque1*dir[1], torque1*dir[2]);
          }

        if (obody2) {
          //dir = u2;
          torque1 = -torque1;
          //dBodyAddRelTorque(obody2, torque1*dir[0], torque1*dir[1], torque1*dir[2]);
          dBodyAddTorque(obody2, torque1*dir[0], torque1*dir[1], torque1*dir[2]);
          }
        }

      // axis 2 //

      if (k2 != 0.0) {
        ang2 = rang[1];
        dang = ang2;
        if (dang > M_PI)  dang -= 2*M_PI;
        if (dang < -M_PI) dang += 2*M_PI;
        torque2 = dang*k2;
        dir = axis2;
        dir.set(0,1,0); 

        if (obody1) {
          //dir = v1;
          //dBodyAddRelTorque(obody1, torque2*dir[0], torque2*dir[1], torque2*dir[2]);
          dBodyAddTorque(obody1, torque2*dir[0], torque2*dir[1], torque2*dir[2]);
          }

        if (obody2) {
          //dir = v2;
          torque2 = -torque2;
          //dBodyAddRelTorque(obody2, torque2*dir[0], torque2*dir[1], torque2*dir[2]);
          dBodyAddTorque(obody2, torque2*dir[0], torque2*dir[1], torque2*dir[2]);
          }
        }

      // axis 3 //

      if (k3 != 0.0) {
        ang3 = rang[2];
        dang = ang3;
        if (dang > M_PI)  dang -= 2*M_PI;
        if (dang < -M_PI) dang += 2*M_PI;
        torque3 = dang*k3;
        dir = axis3;
        dir.set(0,0,1); 

        if (obody1) {
          //dir = w1;
          //dBodyAddRelTorque(obody1, torque3*dir[0], torque3*dir[1], torque3*dir[2]);
          dBodyAddTorque(obody1, torque3*dir[0], torque3*dir[1], torque3*dir[2]);
          }

        if (obody2) {
          //dir = w2;
          torque3 = -torque3;
          //dBodyAddRelTorque(obody2, torque3*dir[0], torque3*dir[1], torque3*dir[2]);
          dBodyAddTorque(obody2, torque3*dir[0], torque3*dir[1], torque3*dir[2]);
          }
        }
      }

    // hinge joint //

    else if (jtype == PM_JOINT_HINGE) {
      PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
      float k, ang, dang, torque;
      hjoint->getForceConst(k);

      if (k == 0.0) {
        continue;
        }

      ang = dJointGetHingeAngle(ojoint);
      dang = ang;
      if (dang > M_PI)  dang -= 2*M_PI;
      if (dang < -M_PI) dang += 2*M_PI;
      torque = -dang*k;
      dJointAddHingeTorque (ojoint, torque);

      energy = 0.5*k*dang*dang;
      joint->addEnergy(upd, energy);

      /*
      fprintf (stderr, ">>> joint %d   ang = %f   dang = %f  torq = %f \n", i, ang,
               dang, torque);
      */
      }
    }

#endif
  }

//*============================================================*
//*==========               stepMotors               ==========*
//*============================================================*
// advance any motors attached to joints. 

void
PmOdeSolver::stepMotors()
 {
#ifdef PM_ODE_SOLVER
  //fprintf (stderr, ">>> PmOdeSolver::stepMotors \n");
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  PmJoint *joint;
  vector<PmMotor*> motors;
  dJointID mjoint;
  float avel, ang, ainc, max_angle;
  int num_axes, axes[3], axis;
  float time = prvd->time;
  float dt = prvd->time_step;
  
  // check each joint for motors // 

  int num_joints = prvd->joints.size();

  for (int i = 0; i < num_joints; i++) {
    joint = prvd->joints[i].pm;

    if (!joint->hasMotors()) {
      continue;
      }

    joint->getMotors(motors);
    //fprintf (stderr, "---- joint[%d]  %d motors ---- \n", i, motors.size());

    for (unsigned int j = 0; j < motors.size(); j++) {
      PmMotor *motor = motors[j];
      PmAngularMotor *amotor = dynamic_cast<PmAngularMotor*>(motor);
      amotor->getAxes(num_axes, axes);
      mjoint = prvd->joints[i].motors[j];

      if (!motor->isEnabled()) {
        continue;
        }

      if (!motor->isActive(time)) { 
        if (motor->isActive()) {
          ode_activateJointMotor(prvd, prvd->joints[i], false);
          //fprintf (stderr, "---- joint[%d] time[%f] de-activate \n", i, time);
          motor->setActive(false);
          }

        continue;
        }
      else if (!motor->isActive()) {
        ode_activateJointMotor(prvd, prvd->joints[i], true);
        motor->setActive(true);
        //fprintf (stderr, "---- joint[%d] time[%f] activate \n", i, time);
        }


      // the angle increment is determined by the //
      // velocity (degs/radians) and time step    //

      for (int k = 0; k < num_axes; k++) {
        axis = axes[k] - 1;
        amotor->getMaxVelocity(k+1, avel);
        ainc = avel * dt;
        //amotor->getAngleInc(k+1, ainc);
        ang = dJointGetAMotorAngle (mjoint, axis);
        ang += ainc;
        amotor->addAngleData(ang);
        amotor->getMaxAngle(k+1, max_angle);
        max_angle = max_angle * (M_PI / 180.0);
        //fprintf (stderr, "\n---- joint[%d] ang[%f] max_angle[%f]\n", i, ang, max_angle); 
        if (fabs(ang) <= max_angle) {
          dJointSetAMotorAngle (mjoint, axis, ang);
          }
        else {
          dJointSetAMotorParam (mjoint, vel_param[k], 0.0);
          dJointSetAMotorParam (mjoint, force_param[k], 0.0);
          motor->setEnabled(false);
          }
        }
      }
    }

#endif
 }

//*============================================================*
//*==========           addRestraintForces           ==========*
//*============================================================*
// add restraint forces.

void
PmOdeSolver::addRestraintForces()
  {
  #ifdef dbg_addRestraintForces
  fprintf (stderr, ">>> addRestraintForces \n");
  #endif

  //===== get the restraints defined for a simulation =====//

  vector<PmRestraint*> restraints;

  simulation->getRestraints(restraints);
  int num = restraints.size();

  if (num == 0) {
    return;
    }

  bool upd = simulation->checkUpdateState();
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  float time = prvd->time;

  PmVector3 center1, center2;
  float power_params[2];

  #ifdef dbg_addRestraintForces
  fprintf (stderr, ">>> num restraints [%d] \n", num);
  #endif

  //===== compute restraint forces =====//

  for (unsigned int i = 0; i < restraints.size(); i++) {
    restraints[i]->compForces(upd, time);
    }
  }

//*============================================================*
//*==========             ode_attach_joints          ==========*
//*============================================================*
// attach bodies to joints.

void
ode_AttachJoints (PmOdeJoint& ojoint, PmBody *body1, dBodyID obody1, PmBody *body2, 
                  dBodyID obody2)
  {
  //fprintf (stderr, "\n###### ode_AttachJoints body1[%x]  body2[%x] \n", body1, body2);
  //fprintf (stderr, "\n###### ode_AttachJoints obody1[%x]  obody2[%x] \n", obody1, obody2);

#ifdef PM_ODE_SOLVER

  if (body1->isGround()) {
    dJointAttach (ojoint.ode, obody2, 0);
    }
  else if (body2->isGround()) {
    dJointAttach (ojoint.ode, obody1, 0);
    }
  else {
    dJointAttach (ojoint.ode, obody1, obody2);
    }

#endif
  }

//*============================================================*
//*============================================================*
// add joints to ode solver.

bool
PmOdeSolver::addJoints(vector<PmJoint*> joints)
  {

#ifdef PM_ODE_SOLVER
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  int num_joints = joints.size();
  PmBody *body1, *body2; 
  #ifdef dbg_PmOdeSolver
  fprintf (stderr, "\n   --- add joints --- \n");
  #endif

  for (int i = 0; i < num_joints; i++) {
    string jname;
    PmOdeJoint ojoint;
    PmVector3 pos;
    dBodyID obody1, obody2;
    vector<PmMotor*> motors;

    PmJoint *joint = joints[i];
    ojoint.pm = joint;
    joint->getName(jname);
    PmJointType jtype = joint->getType();

    // get joint data //

    joint->getPosition (pos);
    joint->getBodies(&body1, &body2);
    joint->getMotors(motors);

    // get ode body associated with pm bodies //

    ode_getBody(prvd, body1, &obody1);
    ode_getBody(prvd, body2, &obody2);
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "   >>> joint [%s]\n", jname.c_str());
    fprintf (stderr, "   >>> num motors [%d]\n", motors.size()); 
    #endif

    // create a two-axis universal joint //

    if (jtype == PM_JOINT_UNIVERSAL) {
      #ifdef dbg_PmOdeSolver
      fprintf (stderr, "       type universal \n");
      #endif
      ojoint.ode = dJointCreateUniversal (prvd->world, 0);
      ode_AttachJoints (ojoint, body1, obody1, body2, obody2);
      dJointSetUniversalAnchor (ojoint.ode, pos[0], pos[1], pos[2]);
      PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);

      if (!ujoint) {
        pm_ErrorReport ("pm> ", "joint \"%s\" failed conversion.", "*", name.c_str());
        return false;
        }

      PmVector3 axis1, axis2;
      ujoint->getAxis(1,axis1);
      dJointSetUniversalAxis1 (ojoint.ode, axis1[0], axis1[1], axis1[2]);
      ujoint->getAxis(2,axis2);
      dJointSetUniversalAxis2 (ojoint.ode, axis2[0], axis2[1], axis2[2]);

      #ifdef dbg_PmOdeSolver
      fprintf (stderr, "       axes:  1 (%g %g %g)\n",axis1[0],axis1[1],axis1[2]); 
      fprintf (stderr, "       axes:  2 (%g %g %g)\n",axis2[0],axis2[1],axis2[2]); 
      #endif
      ode_addJointStop(prvd, ojoint);
      }

    // create a single-axis hinge joint //

    else if (jtype == PM_JOINT_HINGE) {
      #ifdef dbg_PmOdeSolver
      fprintf (stderr, "       type hinge \n");
      #endif
      ojoint.ode = dJointCreateHinge(prvd->world, 0);
      ode_AttachJoints (ojoint, body1, obody1, body2, obody2);
      dJointSetHingeAnchor(ojoint.ode, pos[0], pos[1], pos[2]);
      PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);

      if (!hjoint) {
        pm_ErrorReport ("pm> ", "joint \"%s\" failed conversion.", "*", name.c_str());
        return false;
        }

      PmVector3 axis;
      hjoint->getAxis(axis);
      dJointSetHingeAxis (ojoint.ode, axis[0], axis[1], axis[2]);
      ode_addJointStop(prvd, ojoint);
      }

    // three-axis ball joint //

    else if (jtype == PM_JOINT_BALL) {
      //PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);
      #ifdef dbg_PmOdeSolver
      fprintf (stderr, "       type ball \n");
      fprintf (stderr, "       anchor (%g %g %g) \n", pos[0], pos[1], pos[2]);
      #endif
      ojoint.ode = dJointCreateBall(prvd->world, 0);
      ode_AttachJoints (ojoint, body1, obody1, body2, obody2);
      dJointSetBallAnchor (ojoint.ode, pos[0], pos[1], pos[2]);
      ode_addJointStop(prvd, ojoint);
      }

    else if (jtype == PM_JOINT_WELD) {
      ojoint.ode = dJointCreateFixed(prvd->world, 0);
      ode_AttachJoints (ojoint, body1, obody1, body2, obody2);
      dJointSetFixed (ojoint.ode);
      }

    else {
      string name;
      joint->getName(name);
      pm_ErrorReport ("pm> ", "joint \"%s\" type not supported.", "*", name.c_str());
      return false;
      }

    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "       bodies: [%x] [%x] \n", obody1, obody2);
    #endif
    prvd->joints.push_back(ojoint);
    }

  // create motors for the joints //

  addJointMotors();

  return true;
#endif
  }

//*============================================================*
//*==========             addJointMotors             ==========*
//*============================================================*
// add a motor to a joint.

void 
PmOdeSolver::addJointMotors()
  {
#ifdef PM_ODE_SOLVER
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  PmBody *body1, *body2; 
  PmJoint *joint;
  dBodyID obody1, obody2;
  vector<PmMotor*> motors;

  for (unsigned int i = 0; i < prvd->joints.size(); i++) {
    joint = prvd->joints[i].pm;

    if (!joint->hasMotors()) {
      continue;
      }

    joint->getMotors(motors);
    PmJointType jtype = joint->getType();
    joint->getBodies(&body1, &body2);
    joint->getMotors(motors);
    ode_getBody(prvd, body1, &obody1);
    ode_getBody(prvd, body2, &obody2);

    // create motors for the joint //

    for (unsigned int j = 0; j < motors.size(); j++) {
      int num_axes, axes[3];
      PmOdeJoint mjoint;
      PmVector3 jaxis;
      float force, vel;
      PmMotor *motor = motors[j];
      PmAngularMotor *amotor = dynamic_cast<PmAngularMotor*>(motor);
      amotor->getAxes(num_axes, axes);

      // create ode motor //

      mjoint.ode = dJointCreateAMotor(prvd->world, 0);

      #ifdef PmOdeSolver_addJointMotors_use_euler

      dJointSetAMotorMode (mjoint.ode, dAMotorEuler);

      #else

      dJointSetAMotorMode (mjoint.ode, dAMotorUser);

      ode_AttachJoints (mjoint, body1, obody1, body2, obody2);
      dJointSetAMotorNumAxes (mjoint.ode, num_axes);
      prvd->joints[i].motors[j] = mjoint.ode;

      for (int k = 0; k < num_axes; k++) {
        if (jtype == PM_JOINT_BALL) {
          PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);
          bjoint->getAxis(axes[k], jaxis);
          }
  
        dJointSetAMotorAxis (mjoint.ode, k, 2, jaxis[0], jaxis[1], jaxis[2]);
        //dJointSetAMotorAxis (mjoint.ode, k, axes[k]-1, jaxis[0], jaxis[1], jaxis[2]);

        // set max vel and max force //

        if (amotor->isEnabled()) { 
          amotor->getMaxVelocity(k+1, vel);
          amotor->getMaxForce(k+1, force);
          dJointSetAMotorParam (mjoint.ode, vel_param[k], vel);
          dJointSetAMotorParam (mjoint.ode, force_param[k], force);
          }
        }

      if (jtype == PM_JOINT_HINGE) {
        PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
        float range, hi_ang, low_ang;

        if (hjoint->getRange(range)) {
          hi_ang  =  range*M_PI / 180.0;
          low_ang = -range*M_PI / 180.0;
          dJointSetAMotorParam(mjoint.ode, dParamLoStop, low_ang); 
          dJointSetAMotorParam(mjoint.ode, dParamHiStop, hi_ang); 
          //dJointSetHingeParam (mjoint.ode, dParamLoStop, low_ang); 
          //dJointSetHingeParam (mjoint.ode, dParamHiStop, hi_ang); 
          //fprintf (stderr, ">>> add hinge stop. range = %f \n", range);
          }
        }

      #endif
      }
    }
#endif
  }

//*============================================================*
//*==========             addJointStop               ==========*
//*============================================================*
// add a stop to a joint.

static void 
ode_addJointStop(PmOdeSolverData *prvd, PmOdeJoint& ojoint)
  {
#ifdef PM_ODE_SOLVER
  PmJoint *joint = ojoint.pm;
  PmJointType jtype = joint->getType();
  PmOdeJoint mjoint;

  // joint stops for a ball joint are implemented //
  // using an angular motor.                      //

  if (jtype == PM_JOINT_BALL) {
    PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);
    int num_axes = 0;

    for (int i = 0; i < 3; i++) {
      if (bjoint->hasRange(i+1)) {
        num_axes += 1;
        break; 
        }
      }

    if (!num_axes) {
      return;
      }

    //fprintf (stderr, ">>> add stop \n");
    PmBody *body1, *body2; 
    dBodyID obody1, obody2;
    PmVector3 jaxis;
    float hi_ang, fmax, low_ang, val;

    // get joint information //
    joint->getBodies(&body1, &body2);
    ode_getBody(prvd, body1, &obody1);
    ode_getBody(prvd, body2, &obody2);

    mjoint.ode = dJointCreateAMotor(prvd->world, 0);
    dJointSetAMotorNumAxes(mjoint.ode, 3);
    ode_AttachJoints(mjoint, body1, obody1, body2, obody2);
    fmax = 1.0;

    // set axis 1 //

    bjoint->getAxis(1, jaxis);

    if (bjoint->hasRange(1)) {
      bjoint->getRange(1, val);
      hi_ang  = val*M_PI / 180.0;
      low_ang = -val*M_PI / 180.0;
      }
    else {
      hi_ang  = 0.0; 
      low_ang = 0.0;
      }

    dJointSetAMotorAxis(mjoint.ode, 0, 1, jaxis[0], jaxis[1], jaxis[2]);
    dJointSetAMotorParam(mjoint.ode, dParamFMax, fmax);
    dJointSetAMotorParam(mjoint.ode, dParamLoStop, low_ang); 
    dJointSetAMotorParam(mjoint.ode, dParamHiStop, hi_ang); 
    dJointSetAMotorParam(mjoint.ode, dParamVel, 0.0);
    //dJointSetAMotorParam(mjoint.ode, dParamStopERP, 0.8);
    //dJointSetAMotorParam(mjoint.ode, dParamStopCFM, 0.0);

    // axis 2 //

    if (bjoint->hasRange(2)) {
      bjoint->getRange(2, val);
      hi_ang  = val*M_PI / 180.0;
      low_ang = -val*M_PI / 180.0;
      }
    else {
      hi_ang  = 0.0;
      low_ang = 0.0;
      }

    dJointSetAMotorParam(mjoint.ode, dParamFMax2, fmax);
    dJointSetAMotorParam(mjoint.ode, dParamLoStop2, low_ang); 
    dJointSetAMotorParam(mjoint.ode, dParamHiStop2, hi_ang); 
    dJointSetAMotorParam(mjoint.ode, dParamVel2, 0.0);
    //dJointSetAMotorParam(mjoint.ode, dParamStopERP2, 0.8);
    //dJointSetAMotorParam(mjoint.ode, dParamStopCFM2, 0.0);

    // axis 3 //

    if (bjoint->hasRange(3)) {
      bjoint->getRange(3, val);
      hi_ang  = val*M_PI / 180.0;
      low_ang = -val*M_PI / 180.0;
      }
    else {
      hi_ang  = 0.0;
      low_ang = 0.0;
      }

    bjoint->getAxis(2, jaxis);
    dJointSetAMotorAxis(mjoint.ode, 2, 2, jaxis[0], jaxis[1], jaxis[2]);
    dJointSetAMotorParam(mjoint.ode, dParamFMax3, fmax);
    dJointSetAMotorParam(mjoint.ode, dParamLoStop3, low_ang); 
    dJointSetAMotorParam(mjoint.ode, dParamHiStop3, hi_ang); 
    dJointSetAMotorParam(mjoint.ode, dParamVel3, 0.0);
    //dJointSetAMotorParam(mjoint.ode, dParamStopERP3, 0.8);
    //dJointSetAMotorParam(mjoint.ode, dParamStopCFM3, 0.0);

    dJointSetAMotorMode(mjoint.ode, dAMotorEuler);

    /*
    for (int i = 0; i < 3; i++) {
      if (bjoint->hasRange(i+1)) {
        bjoint->getRange(i+1, val);
        val = val*M_PI / 180.0;
        bjoint->getAxis(i+1, jaxis);
        dJointSetAMotorAxis(mjoint.ode, num_axes, 2, jaxis[0], jaxis[1], jaxis[2]);
        fprintf (stderr, ">>> add stop %f \n", val);
        fprintf (stderr, ">>> axis %f %f %f \n", jaxis[0], jaxis[1], jaxis[2]);

        if (num_axes == 0) { 
          dJointSetAMotorParam(mjoint.ode, dParamHiStop, val); 
          dJointSetAMotorParam(mjoint.ode, dParamFMax, fmax);
          dJointSetAMotorParam(mjoint.ode, dParamVel, 0.0);
          }
        else if (num_axes == 1) { 
          dJointSetAMotorParam(mjoint.ode, dParamHiStop2, val); 
          dJointSetAMotorParam(mjoint.ode, dParamFMax2, fmax);
          dJointSetAMotorParam(mjoint.ode, dParamVel2, 0.0);
          }
        else if (num_axes == 2) { 
          dJointSetAMotorParam(mjoint.ode, dParamHiStop3, val); 
          dJointSetAMotorParam(mjoint.ode, dParamFMax3, fmax);
          dJointSetAMotorParam(mjoint.ode, dParamVel3, 0.0);
          }
       
        num_axes += 1;
        }
      }
    */
    }

  else if (jtype == PM_JOINT_HINGE) {
    PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
    float range, hi_ang, low_ang;

    if (hjoint->getRange(range)) {
      hi_ang  =  range*M_PI / 180.0;
      low_ang = -range*M_PI / 180.0;
      dJointSetHingeParam (ojoint.ode, dParamLoStop, low_ang); 
      dJointSetHingeParam (ojoint.ode, dParamHiStop, hi_ang); 
      }
    }

  else if (jtype == PM_JOINT_UNIVERSAL) {
    //PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);
    }

#endif
  }

//*============================================================*
//*==========             activateJointMotor         ==========*
//*============================================================*
// activate/deactivate a joint motor.

static void 
ode_activateJointMotor(PmOdeSolverData *prv_data, PmOdeJoint& ojoint, bool activate)
  {
#ifdef PM_ODE_SOLVER
  float force, vel;
  int num_axes, axes[3];
  vector<PmMotor*> motors;
  PmJoint *joint = ojoint.pm;
  joint->getMotors(motors);

  for (unsigned int j = 0; j < motors.size(); j++) {
    dJointID mjoint = ojoint.motors[j];
    PmMotor *motor = motors[j];

    if (!motor->isEnabled()) {
      continue;
      }

    PmAngularMotor *amotor = dynamic_cast<PmAngularMotor*>(motor);
    amotor->getAxes(num_axes, axes);

    for (int k = 0; k < num_axes; k++) {
      if (activate) { 
        amotor->getMaxVelocity(k+1, vel);
        amotor->getMaxForce(k+1, force);
        dJointSetAMotorParam (mjoint, vel_param[k], vel);
        dJointSetAMotorParam (mjoint, force_param[k], force);
        }
      else {
        dJointSetAMotorParam (mjoint, vel_param[k], 0.0);
        dJointSetAMotorParam (mjoint, force_param[k], 0.0);
        }
      }
    }
#endif
  }

//*============================================================*
//*==========             getBody                    ==========*
//*============================================================*
// get an ode body id for the given PmBody.

static void
ode_getBody(PmOdeSolverData *prvd, PmBody *body, dBodyID *obody)
  {
#ifdef PM_ODE_SOLVER
  int num_bodies = prvd->bodies.size();

  for (int i = 0; i < num_bodies; i++) {
    if (body == prvd->bodies[i].pm) {
      *obody = prvd->bodies[i].ode;
      return;
      }
    }

  *obody = NULL;
#endif
  }

//*============================================================*
//*==========             getJoint                   ==========*
//*============================================================*
// get an ode joint id for the given PmJoint.

static void
ode_getJoint(PmOdeSolverData *prvd, PmJoint *joint, dJointID *ojoint)
  {
#ifdef PM_ODE_SOLVER
  int num_joints = prvd->joints.size();

  for (int i = 0; i < num_joints; i++) {
    if (joint == prvd->joints[i].pm) {
      *ojoint = prvd->joints[i].ode;
      return;
      }
    }

  *ojoint = NULL;
#endif
  }

//*============================================================*
//*==========                 update                 ==========*
//*============================================================*
// update body, joints and trace states. note that we still need
// to update the transformation for a body to update the geoemtry 
// of any potentials defined for the body.

void
PmOdeSolver::update() 
  {
#ifdef PM_ODE_SOLVER
  #define ndbg_PmOdeSolver
  #ifdef dbg_PmOdeSolver
  PmOdeSolverData *prvddbg = (PmOdeSolverData*)prv_data;
  fprintf (stderr, "\n-------------------------------\n");
  fprintf (stderr, ">>>>>> ODE Solver: update\n");
  fprintf (stderr, "   >>> num bodies[%d]\n", prvddbg->bodies.size());
  fprintf (stderr, "   >>> num joints[%d]\n", prvddbg->joints.size());
  #endif

  if (!initialzed) {
    return;
    }

  const dReal *pos, *r, *vel, *avel;
  PmVector3 com, rvec, disp, vec, linvel, angvel;
  bool upd;
  PmMatrix3x3 m; 
  PmQuaternion q;
  PmMassProperties props;
  float ke, total_ke;
  PmOdeSolverData *prvd = (PmOdeSolverData*)prv_data;
  float time = prvd->time;

  // check to see if we should update state //
  // for storing results and display.       //

  upd = simulation->checkUpdateState();

  // update bodies //

  total_ke = 0.0;

  for (unsigned int i = 0; i < prvd->bodies.size(); i++) {
    #ifdef dbg_PmOdeSolver
    fprintf (stderr, "------ body [%d] ------ \n", i);
    #endif
    PmBody *body = prvd->bodies[i].pm;
    dBodyID obody = prvd->bodies[i].ode;

    if ((body->getType() == PM_BODY_GROUND) || (body->getType() == PM_BODY_STATIC)) {
      #ifdef dbg_PmOdeSolver
      fprintf (stderr, ">>> ground \n");
      #endif
      continue;
      }

    body->getCenterOfMass(com);
    pos = dBodyGetPosition (obody);
    r = dBodyGetRotation (obody);

    disp[0] = pos[0] - com[0];
    disp[1] = pos[1] - com[1];
    disp[2] = pos[2] - com[2];

    #ifdef PmOdeSolver_update_comp_other_rot 
    rot[0][0] = r[0]; rot[0][1] = r[4]; rot[0][2] = r[8];
    rot[1][0] = r[1]; rot[1][1] = r[5]; rot[1][2] = r[9];
    rot[2][0] = r[2]; rot[2][1] = r[6]; rot[2][2] = r[10];

    pm_MathExtractRotations(rot, ang);
    rvec[0] = 180.0*ang[0] / M_PI;
    rvec[1] = 180.0*ang[1] / M_PI;
    rvec[2] = 180.0*ang[2] / M_PI;

    #ifdef dbg_PmOdeSolver
    fprintf (stderr, ">>> pos  (%g %g %g) \n", pos[0], pos[1], pos[2]);
    fprintf (stderr, ">>> disp (%g %g %g) \n", disp[0], disp[1], disp[2]);
    fprintf (stderr, ">>> ang  (%g %g %g) \n", ang[0], ang[1], ang[2]);
    #endif
    #endif

    m(0,0) = r[0]; m(1,0) = r[4]; m(2,0) = r[8];
    m(0,1) = r[1]; m(1,1) = r[5]; m(2,1) = r[9];
    m(0,2) = r[2]; m(1,2) = r[6]; m(2,2) = r[10];

    // convert the rotation matrix into a quaternion //
    //pm_MathMatrixToQuaternion(m, q);


    // set kinetic energy //

    body->getMassProps(props);
    vel = dBodyGetLinearVel(obody);
    linvel[0] = vel[0];
    linvel[1] = vel[1];
    linvel[2] = vel[2];
    ke = 0.5*props.mass*linvel*linvel;

    avel = dBodyGetAngularVel(obody);
    angvel[0] = avel[0];
    angvel[1] = avel[1];
    angvel[2] = avel[2];
    rvec = props.inertia * angvel;
    ke += 0.5*rvec*rvec;
    total_ke += ke;
    body->addKineticEnergy(upd, ke);

    // set the simulation results //

    body->addSimResults(upd, time, disp, rvec, linvel, angvel, m, q);

    // set body state //

    PmRigidBodyState state;

    for (int i = 0; i < 3; i++) {
      state.velocity[i] = vel[i];
      state.angular_velocity[i] = avel[i];
      state.displacement[i] = disp[i];
      }
 
    state.rotation_matrix = m;
    body->setState(state);
    }

  //fprintf (stderr, ">>>> total ke = %g \n", total_ke);

  // update joints //

  PmBody *body1, *body2, *body;
  dBodyID obody;
  PmVector3 jpos, xjpos;
  PmXform body_xform;
  dJointID ojoint;

  for (unsigned int i = 0; i < prvd->joints.size(); i++) {
    //fprintf (stderr, "------ joint[%d] ------ \n", i);
    PmJoint *joint = prvd->joints[i].pm;
    PmJointType jtype = joint->getType();
    joint->getPosition (jpos);
    joint->getBodies(&body1, &body2);
    ojoint = prvd->joints[i].ode;

    if (body1->isGround()) {
      body = body2;
      }
    else {
      body = body1;
      }

    ode_getBody(prvd, body, &obody);
    ode_xformPoint(body, obody, jpos, xjpos);
    /*
    fprintf (stderr, ">>> jpos=%g %g %g -> xjpos=%g %g %g\n", jpos[0], jpos[1], jpos[2],
                                                             xjpos[0], xjpos[1], xjpos[2]);
    */
    r = dBodyGetRotation (obody);
    PmMatrix3x3 m;
    m(0,0) = r[0]; m(1,0) = r[4]; m(2,0) = r[8];
    m(0,1) = r[1]; m(1,1) = r[5]; m(2,1) = r[9];
    m(0,2) = r[2]; m(1,2) = r[6]; m(2,2) = r[10];

    body_xform.translation = xjpos;
    body->getCenterOfMass (com);
    body_xform.center = com;
    body_xform.setRotation(m);

    joint->update(xjpos);
    joint->update(body_xform);

    //dReal a = dJointGetHingeAngle(mjoint);
    //fprintf (stderr, "---- joint[%d] ang[%f] \n", i, ang);

    if (jtype == PM_JOINT_HINGE) {
      PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
      dReal a = dJointGetHingeAngle(ojoint);
      hjoint->updateAngle(a);
      }
    }

  // update objects not needed for the    //
  // simulation every saved state upate.  // 

  if (!upd) {
    return;
    }

  // update geometies //

  vector<PmGeometry*> geoms;
  simulation->getGeometries(geoms);
  int ngeoms = geoms.size();

  for (int i = 0; i < ngeoms; i++) {
    geoms[i]->getBody(&body);
    ode_getBody(prvd, body, &obody);
    body->getCenterOfMass(com);
    pos = dBodyGetPosition (obody);

    if (!obody) {
      fprintf (stderr, "\n**** WARNING: geometry does not have an obody. \n");
      continue;
      }

    r = dBodyGetRotation (obody);
    //ode_xformPoint(body, obody, com, xjpos);

    PmMatrix3x3 m;
    m(0,0) = r[0]; m(1,0) = r[4]; m(2,0) = r[8];
    m(0,1) = r[1]; m(1,1) = r[5]; m(2,1) = r[9];
    m(0,2) = r[2]; m(1,2) = r[6]; m(2,2) = r[10];

    body_xform.translation[0] = pos[0] - com[0];
    body_xform.translation[1] = pos[1] - com[1];
    body_xform.translation[2] = pos[2] - com[2];

    body->getCenterOfMass (com);
    body_xform.center = com;
    body_xform.setRotation(m);
    geoms[i]->update(body_xform);
    geoms[i]->display(true);
    }

  // update traces //

  vector<PmTrace*> traces;
  PmVector3 loc, tloc;
  bool global_frame = false;
  dVector3 gvec;

  simulation->getTraces(traces);
  int num = traces.size();
  //fprintf (stderr, ">>>>>> PmOdeSolver::update   ntraces %d \n", num);

  for (int i = 0; i < num; i++) {
    traces[i]->getGlobalFrame(global_frame);
    traces[i]->getObj((void**)&body);
    ode_getBody(prvd, body, &obody);
    body->getCenterOfMass(com);
    traces[i]->getPoint(loc);
    //fprintf (stderr, "   >>> loc %f %f %f \n", loc[0], loc[1], loc[2]);
    //fprintf (stderr, "   >>> com %f %f %f \n", com[0], com[1], com[2]);

    if (global_frame) {
      pos = dBodyGetPosition (obody);
      //r = dBodyGetRotation (obody);
      vec = loc - com;
      dBodyVectorToWorld (obody, vec[0], vec[1], vec[2], gvec);

      for (int j = 0; j < 3; j++) {
        tloc[j] = pos[j] + gvec[j];
        }
      }
    else {
      dBodyGetRelPointPos (obody, loc[0], loc[1], loc[2], gvec);

      for (int j = 0; j < 3; j++) {
        tloc[j] = gvec[j];
        }
      }

    traces[i]->addPoint(tloc);
    traces[i]->display(true);
    }

  // update measurements //

  vector<PmMeasurement*> msr;
  vector<PmVector3> points;
  vector<PmVector3> upd_points;
  vector<void*> pobjs;
  PmVector3 ploc;
  global_frame = true;
  simulation->getMeasurements(msr);
  num = msr.size();
  //fprintf (stderr, ">>>>>> PmOdeSolver::update   nmsr %d \n", num);

  for (int i = 0; i < num; i++) {
    //msr[i]->getGlobalFrame(global_frame);
    msr[i]->getObjs(pobjs);
    msr[i]->getPoints(points);
    upd_points.clear();

    for (unsigned int j = 0; j < points.size(); j++) {
      loc = points[j];
      body = (PmBody*)pobjs[j];
      ode_getBody(prvd, body, &obody);
      body->getCenterOfMass(com);

      if (global_frame) {
        pos = dBodyGetPosition (obody);
        r = dBodyGetRotation (obody);
        vec = loc - com;
        dBodyVectorToWorld (obody, vec[0], vec[1], vec[2], gvec);

        for (int k = 0; k < 3; k++) {
          ploc[k] = pos[k] + gvec[k];
          }
        }
      else {
        dBodyGetRelPointPos (obody, loc[0], loc[1], loc[2], gvec);

        for (int k = 0; k < 3; k++) {
          ploc[k] = gvec[k];
          }
        }

      //fprintf (stderr, "ploc %f %f %f \n", ploc[0], ploc[1], ploc[2]);
      upd_points.push_back(ploc);
      }

    msr[i]->updatePoints(upd_points);
    msr[i]->display();
    }
  }

//*============================================================*
//*==========                 XformPoint             ==========*
//*============================================================*
// transform a point on a body.

static void
ode_xformPoint(PmBody *body, dBodyID obody, PmVector3& pt, PmVector3& xpt) 
  {
#ifdef PM_ODE_SOLVER
  const dReal *pos; 
  dReal gvec[3]; 
  PmVector3 com, vec;

  if (body->getType() == PM_BODY_GROUND) {
    xpt = pt;
    return; 
    }

  //fprintf (stderr, ">>>>>> ode_xformPoint \n");
  body->getCenterOfMass(com);
  pos = dBodyGetPosition (obody);
  vec = pt - com;
  dBodyVectorToWorld (obody, vec[0], vec[1], vec[2], gvec);
  //fprintf (stderr, "   >>> gvec=%g %g %g\n", gvec[0], gvec[1], gvec[2]);
  //fprintf (stderr, "   >>> pos=%g %g %g\n", pos[0], pos[1], pos[2]);

  for (int i = 0; i < 3; i++) {
    xpt[i] = pos[i] + gvec[i];
    }

#endif
  }

#endif

}

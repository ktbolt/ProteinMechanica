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

#include "db.h"
#include "rbsim.h"

#define ndbg_PmRigidSimulation

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmRigidSimulation::PmRigidSimulation (const string name) 
  { 
  #ifdef dbg_PmRigidSimulation
  fprintf (stderr, ">>> PmRigidSimulation:: ctor  [%s] \n", name.c_str());
  #endif
  this->name = name;
  this->type = PM_SIMULATION_RIGID;
  solver = NULL;
  current_time = 0.0;
  time_step = 0.01;
  simulation_active = true;
  num_inc = 0;
  inc_count = 0;

  replay_active = false;
  replay_start = 0;
  replay_end = 0;
  replay_step = 0;
  replay_inc = 1;

  num_steps = 0;
  state_save_freq = 1;
  momentum_off = false;
  damping = false;

  PmSimulation::simulations.push_back(this);
  }

//*============================================================*
//*==========              addBody                   ==========*
//*============================================================*
// add a body to the simulation.

bool
PmRigidSimulation::addBody(PmBody *body)
  {
  #ifdef dbg_PmRigidSimulation
  fprintf (stderr, ">>> PmRigidSimulation::addBody  \n");
  #endif
  bodies.push_back(body);

  PmRigidSimObj rsobj;
  body->getName(rsobj.name); 
  rsobj.body = body;
  simulation_objects.push_back(rsobj);
  return true; 
  }

//*============================================================*
//*==========              getBody                   ==========*
//*============================================================*
// get a body from the simulation.

bool
PmRigidSimulation::getBody(const string str, PmBody **body)
  {
  *body = NULL;
  string bname;
  int n = bodies.size();

  for (int i = 0; i < n; i++) {
    PmBody *bdy = bodies[i];
    bdy->getName(bname);

    if (bname == str) {
      *body = bdy;
      return true;
      }
    }

  return false;
  }

void 
PmRigidSimulation::getBodies(vector<PmBody*>& blist) { 
  blist = bodies; 
  }

void 
PmRigidSimulation::setDamping(const bool flag) { 
  damping = flag; 
  }

bool 
PmRigidSimulation::getDamping() { 
  return damping; 
  }

void 
PmRigidSimulation::getJoints(vector<PmJoint*>& jlist) { 
  jlist = joints; 
  }

void 
PmRigidSimulation::setMomentumOff(const bool flag) { 
  momentum_off = flag; 
  }

bool 
PmRigidSimulation::getMomentumOff() { 
  return momentum_off; 
  }

void 
PmRigidSimulation::setSolver(PmRigidBodySolver *solver) { 
  this->solver = solver; 
  }

void 
PmRigidSimulation::getTime(double& time) { 
  time = this->current_time; 
  }

void 
PmRigidSimulation::getTimeStep(double& step) { 
  step = this->time_step; 
  }


//*============================================================*
//*==========              initalize                 ==========*
//*============================================================*
// initialize simulation.

void
PmRigidSimulation::initialize() 
  { 
  #ifdef dbg_PmRigidSimulation
  fprintf (stderr, ">>> PmRigidSimulation:: initialize [%s] \n", name.c_str());
  fprintf (stderr, "    mass scale[%g] \n", pmSystem.getUnitsMassScale());
  #endif

  if (!solver) {
    pm_ErrorReport (PM, "no solver set for simulation \"%s\".", "*", name.c_str());
    return;
    }

  // initialize solver //

  solver->initialize(this, time_step);
  }

//*============================================================*
//*==========              addForce                  ==========*
//*============================================================*
// add a force to a body of the simulation.

bool 
PmRigidSimulation::addForce (const string str, PmForce *force, PmTimeInterval& time)
  {
  PmBody *body;

  if (!getBody(str, &body)) {
    return false;
    }

  PmForceType type = force->getType();
  force->setTime(time);

  if (type == PM_FORCE_EXPLICIT) {
    PmExplicitForce *eforce = dynamic_cast<PmExplicitForce*>(force);
    eforce->setObj(body);
    explicit_forces.push_back(eforce);
    }

  else if (type == PM_FORCE_RANDOM) {
    PmRandomForce *rforce = dynamic_cast<PmRandomForce*>(force);
    rforce->setObj(body);
    random_forces.push_back(rforce);
    }

  return true;
  }

//*============================================================*
//*==========              getForces                 ==========*
//*============================================================*
// get the forces defined for the simulation.

void
PmRigidSimulation::getForces(vector<PmExplicitForce*>& forces) {
  forces = this->explicit_forces; 
  }

//*============================================================*
//*==========              getForces                 ==========*
//*============================================================*
// get the forces defined for the simulation.

void
PmRigidSimulation::getForces(vector<PmRandomForce*>& forces) {
  forces = this->random_forces;
  }

//*============================================================*
//*==========              addGeometry               ==========*
//*============================================================*
// add a geometry to the simulation.

bool 
PmRigidSimulation::addGeometry(PmGeometry* geom) 
  {
  geometries.push_back(geom);
  return true;
  }

//*============================================================*
//*==========              getGeometries             ==========*
//*============================================================*
// get the geometries defined for the simulation.

void
PmRigidSimulation::getGeometries(vector<PmGeometry*>& geoms) {
  geoms = this->geometries;
  }

//*============================================================*
//*==========              halt                      ==========*
//*============================================================*
// halt the simulation.

void
PmRigidSimulation::halt()
  {
  num_inc = 0;
  inc_count = 0;
  }

//*============================================================*
//*==========              addJoint                  ==========*
//*============================================================*
// add a joint to the simulation.

bool
PmRigidSimulation::addJoint (PmJoint *joint) {
  joints.push_back (joint);
  return true;
  }

//*============================================================*
//*==========              getJoint                  ==========*
//*============================================================*
// get a joint from the simulation.

bool
PmRigidSimulation::getJoint(const string str, PmJoint **joint)
  {
  *joint = NULL;
  string jname;
  int n = joints.size();

  for (int i = 0; i < n; i++) {
    PmJoint *jnt = joints[i];
    jnt->getName(jname);

    if (jname == str) {
      *joint = jnt;
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========              addMeasurement            ==========*
//*============================================================*
// add a measurement to a simulation.

bool 
PmRigidSimulation::addMeasurement(PmMeasurement *msr) {
  measurements.push_back(msr);
  return true;
  }

//*============================================================*
//*==========              getMeasurement            ==========*
//*============================================================*
// get a measurement by name for a simulation.

void
PmRigidSimulation::getMeasurement(const string name, PmMeasurement **msr)
  {
  string mname;
  *msr = NULL;

  for (unsigned int i = 0; i < this->measurements.size(); i++) {
    this->measurements[i]->getName(mname);

    if (name == mname) {
      *msr = this->measurements[i];
      return;
      }
    }
  }

//*============================================================*
//*==========              getMeasurements           ==========*
//*============================================================*
// get the measurements for a simulation.

void
PmRigidSimulation::getMeasurements(vector<PmMeasurement*>& msr) {
  msr = this->measurements;
  }

//*============================================================*
//*==========          printModelInformation         ==========*
//*============================================================*
// print model information.                     

void 
PmRigidSimulation::printModelInformation(const bool pbodies, const bool pjoints)
  {
  int num;

  if (pbodies) {
    string bname, ptypes, btype, pname;
    PmPhysicalObj *pobj;
    PmObjectType ptype;

    num = bodies.size();
    fprintf (stderr, "\n---------- bodies ----------\n");
    fprintf (stderr, "number of bodies = %d \n", num);

    for (int i = 0; i < num; i++) {
      PmBody *body = bodies[i];
      body->getName(bname);
      body->getType(btype);
      body->getPhysicalObj(&pobj);

      fprintf (stderr, "%2d:  name = %s  type = %s ", i+1, bname.c_str(), 
               btype.c_str());
  
      if (pobj) {
        pobj->getObjectType(ptypes);
        pobj->getName(pname);
        ptype = pobj->getObjectType();
        fprintf (stderr, "  derived from: name = %s  type = %s ",pname.c_str(), 
                 ptypes.c_str());
        }

      fprintf (stderr, "\n");
      }

    fprintf (stderr, "\n");
    }


  // print joint information //

  if (pjoints) {
    PmBody *body1, *body2;
    string jname, jtype, bname1, bname2;
    PmVector3 pos;
    vector<float> k;
    int nk;

    num = joints.size();
    fprintf (stderr, "\n---------- joints ----------\n");
    fprintf (stderr, "number of joints = %d \n", num);

    for (int i = 0; i < num; i++) {
      PmJoint *joint = joints[i];
      joint->getName(jname);
      joint->getType(jtype);
      joint->getBodies(&body1, &body2);
      joint->getPosition(pos);

      fprintf (stderr, "%2d:  name = %s  type = %s ", i+1, jname.c_str(), 
               jtype.c_str());

      if (body1 && body2) {
        body1->getName(bname1);
        body2->getName(bname2);
        fprintf (stderr, " bodies = { %s %s } ", bname1.c_str(), bname2.c_str());
        }

      fprintf (stderr, "\n     ");
      fprintf (stderr, "position = [ %g, %g %g ] ", pos[0], pos[1], pos[2]); 

      fprintf (stderr, "\n     ");
      joint->getForceConsts(k);
      nk = k.size(); 

      if (nk) {
        if (nk == 1) {
          fprintf (stderr, "force_const = "); 
          }
        else { 
          fprintf (stderr, "force_const = { "); 
          }

        for (int j = 0; j < nk; j++) {
          fprintf (stderr, "%g ", k[j]); 
          }

        if (nk != 1) {
          fprintf (stderr, "} "); 
          }
        }

      fprintf (stderr, "\n");
      }

    fprintf (stderr, "\n");
    }
  }

//*============================================================*
//*==========              addMotor                  ==========*
//*============================================================*
// add a motor to a simulation.

bool
PmRigidSimulation::addMotor(PmMotor *motor, PmTimeInterval& time)
  {
  motor->setTime(time);
  motor->setEnabled(true);
  motors.push_back(motor);
  return true;
  }

//*============================================================*
//*==========              getMotor                  ==========*
//*============================================================*
// get a motor from a simulation.

void 
PmRigidSimulation::getMotor(const string name, PmMotor **motor)
  {
  string mname;
  *motor = NULL;

  for (unsigned int i = 0; i < this->motors.size(); i++) {
    this->motors[i]->getName(mname);

    if (name == mname) {
      *motor = this->motors[i];
      return;
      }
    }
  }

//*============================================================*
//*==========              setSolverParam            ==========*
//*============================================================*
// set a solver parameter.             

bool 
PmRigidSimulation::setSolverParam(const string name, const string val)
  {
  if (!solver) {
    return false;
    }

  return solver->setParam (name, val);
  }

//*============================================================*
//*==========                 replay                 ==========*
//*============================================================*
// replay simulation.   

void
PmRigidSimulation::replay()
  {
  //fprintf (stderr, ">>> PmRigidSimulation[%s]_replay \n", name.c_str());
  //fprintf (stderr, "    num steps [%d] \n", num_steps);
  //fprintf (stderr, "    replay step [%d] \n", replay_step);
  //fprintf (stderr, "    replay active [%d] \n", replay_active);
  //fprintf (stderr, "    state freq [%d] \n", state_save_freq); 

  if (!replay_active) {
    return;
    }

  if (state_save_freq*replay_step >= num_steps-1) {
    replay_active = false;
    replay_step = 0;
    return;
    }

  for (unsigned int j = 0; j < bodies.size(); j++) {
    PmBody *body = bodies[j];
    body->setSimResults(replay_step);
    }

  replay_step += replay_inc;
  pmSystem.updateGraphics(true);
  }

//*============================================================*
//*==========            incrementTimeStep           ==========*
//*============================================================*
// increment the time step.

bool
PmRigidSimulation::incrementTimeStep()
  {
  bool res = false;
  write_results = false;
  bool verbose = pmSystem.getCmdVerbose();
  //fprintf (stderr, ">>>>>> PmRigidSimulation::incrementTimeStep name=%s\n",name.c_str());

  //===== evaluate a boolean while expression =====//

  if (while_exp) {
    res = while_exp->evaluate();

    if (!res) {
      write_results = true;
      while_exp = NULL;
      }
    }

  //===== check is inc_count is > 0 =====//

  else {

    if (this->inc_count == 0) { 
      res = false;
      }

    else {
      this->inc_count--;
      res = true;

      if (this->inc_count == 0) { 
        write_results = true;
        }
      }
    }

  if (write_results) {
    if (verbose) fprintf (stderr, "\n");
    this->writeResults();
    write_results = false;
    pmSystem.setCmdWait(false);
    }

  return res;
  }

//*============================================================*
//*==========                 step                   ==========*
//*============================================================*
// perform a time step.   

void 
PmRigidSimulation::step(int nsteps)
  {
  if (!solver) {
    return;
    }

  if (nsteps != 1) {
    for (int i = 0; i < nsteps; i++) {
      this->step(1);
      }
    return;
    }

  #define ndbg_PmRigidSimulation
  #ifdef dbg_PmRigidSimulation
  fprintf (stderr, ">>> PmRigidSimulation[%s]::step \n", name.c_str());
  fprintf (stderr, "    num inc    [%d] \n", num_inc); 
  fprintf (stderr, "    state freq [%d] \n", state_save_freq); 
  #endif

  //===== save initial state =====//

  if (inc_count && (current_time == 0.0)) {
    time.push_back(current_time);
    //this->writeResults();
    num_state_upd += 1;
    updateEnergy();
    }

  // check if the solver should take a time step //

  //if (this->inc_count) {

  if (incrementTimeStep()) {
    simulation_active = true;
    solver->step();
    this->num_steps++;
    this->current_time += this->time_step;
    //fprintf (stderr, ">>> PmRigidSimulation[%s]::step \n", name.c_str());

    if (this->checkUpdateState() || (num_steps == 1)) { 
      this->updateDisplay();
      pmSystem.updateGraphics(true);

      if (!silent_step_) {
        fprintf (stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        fprintf (stderr, ">>> step [%d]  t [%g] ", num_steps, current_time);
        }

      if (num_steps != 1) {
        num_state_upd += 1;
        time.push_back(current_time);
        updateEnergy();
        }
       
      /*
      else if ((num_steps == 1) && (state_save_freq == 0)) {
        time.push_back(current_time);
        updateEnergy();
        }
      */
      }

    /*
    if (write_results) {
      fprintf (stderr, "\n");
      this->writeResults();
      write_results = false;
      pmSystem.setCmdWait(false);
      }
    */
    }
  else {
    simulation_active = false;
    }
  }

//*============================================================*
//*==========            setTimeStep                 ==========*
//*============================================================*
// set time step.

void 
PmRigidSimulation::setTimeStep(const double step)
  {
  this->time_step = step;
  solver->setTimeStep(time_step);
  }

//*============================================================*
//*==========            updateEnergy                ==========*
//*============================================================*
// update energy. 

void
PmRigidSimulation::updateEnergy()
  {
  /*
  fprintf (stderr, "\n>>> PmRigidSimulation:: update energy \n"); 
  fprintf (stderr, ">>> write_energy = %d \n", write_energy); 
  */

  if (!write_energy) {
    return;
    }

  PmPotentialType ptype;
  vector<PmInteraction*> interactions;
  this->getInteractions(interactions);

  for (unsigned int i = 0; i < interactions.size(); i++) {
    interactions[i]->getPotentialType(ptype); 

    if ((ptype == PM_POTENTIAL_CONTACT) && write_contact_energy) {
      interactions[i]->addEnergy(current_time);
      }
    else {
      interactions[i]->addEnergy(current_time);
      }
    }
  }

////////////////////////////////////////////////////////////////
//                    g r a p h i c s                        //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              displayForces             ==========*
//*============================================================*
// display forces.           

void
PmRigidSimulation::displayForce(PmExplicitForce *force, PmVector3& loc)
  {
/*
  fprintf (stderr, ">>> PmRigidSimulation[%s]::displayForce \n", name.c_str());
  fprintf (stderr, "    loc (%g %g %g) \n", loc[0], loc[1], loc[2]); 
  fprintf (stderr, "    scale[%g] \n", force.scale); 

 
    string geom_name;
    geom_name = "force[" + name + "]";
    PmVector3 *verts = new PmVector3[2];
    verts[0] = loc;
    verts[1] = loc + force.scale*force.direction;
    PmGraphicsLine *geometry = new PmGraphicsLine(geom_name, 2, verts);

    geometry->display();
*/
  }

//*============================================================*
//*==========          getGraphicsGeometry           ==========*
//*============================================================*
// get a graphics geometry object.

void
PmRigidSimulation::getGraphicsGeometry(const string name, PmGraphicsGeometry **geom)
  {
  *geom = NULL;

  for (unsigned int i = 0; i < graphics_geometry.size(); i++) {
    if (name == graphics_geometry[i]->name) {
      *geom = graphics_geometry[i];
      }
    }
  }

//*============================================================*
//*==========              displayRestraint          ==========*
//*============================================================*
// display a restraint.

void
PmRigidSimulation::displayRestraint(const string name, PmVector3& pt1, PmVector3& pt2,
                                    const bool show)
  {
  PmGraphicsRestraint *geom;
  getRestraintGeometry(name, &geom);

  if (!geom) {
    geom = new PmGraphicsRestraint(name, pt1, pt2);
    geom->display();
    graphics_geometry.push_back(geom);
    }
  else {
    geom->update(pt1, pt2);
    }
  }

//*============================================================*
//*==========          getRestraintGeometry          ==========*
//*============================================================*
// get a graphics restraint geometry object.

void
PmRigidSimulation::getRestraintGeometry (const string name, PmGraphicsRestraint **rgeom)
  {
  *rgeom = NULL;
  PmGraphicsGeometry *geom = NULL;
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    fprintf (stderr, ">>> getRestraintGeometry:: name[%s] not found. \n", name.c_str());
    return;
    }

  *rgeom = dynamic_cast<PmGraphicsRestraint*>(geom);
  }

//*============================================================*
//*==========              addTrace                  ==========*
//*============================================================*
// add trace.                                     

bool 
PmRigidSimulation::addTrace(PmTrace *trace)
  {
  traces.push_back(trace);
  return true;
  }

//*============================================================*
//*==========              getTrace                  ==========*
//*============================================================*
// get a trace by name for a simulation.

void 
PmRigidSimulation::getTrace(const string name, PmTrace **trace)
  {
  string tname;
  *trace = NULL;

  for (unsigned int i = 0; i < this->traces.size(); i++) {
    this->traces[i]->getName(tname);

    if (name == tname) {
      *trace = this->traces[i];
      return;
      }
    }
  }

//*============================================================*
//*==========              getTraces                 ==========*
//*============================================================*
// get the traces for a simulation.

void 
PmRigidSimulation::getTraces(vector<PmTrace*>& traces) {
  traces = this->traces;
  }

////////////////////////////////////////////////////////////////
//              r e a d   s t a t e                          //
//////////////////////////////////////////////////////////////

void 
PmRigidSimulation::readState(const string file_name)
  {
  char fname[1000], line[10000], *s;
  FILE *fp;
  int n, len;
  int num_bodies;
  vector<string> body_names;
  vector<PmRigidBodyState> states;
  PmRigidBodyState state;
  float time; 
  PmVector3 disp, rot, vel, avel;
  PmMatrix3x3 mat;
  string str;
  //fprintf (stderr, "------ read state ------ \n"); 

  sprintf (fname, "%s_state.pm", file_name.c_str());
  //fprintf (stderr, ">>> file_name = %s  \n", fname);
  len = 10000;

  if ((fp = fopen(fname, "r")) == NULL) {
    pm_ErrorReport (PM, "can't open file \"%s\".", "*", fname);
    return;
    }

  while (fgets(line, len, fp) != NULL) {
    //fprintf (stderr, ">>> line=\"%s\" \n", line);

    if (line[0] == '#') {
      continue;
      }

    else if (sscanf(line, "number of bodies %d", &num_bodies) ==1) { 
      //fprintf (stderr, ">>> number of bodies=%d \n", num_bodies);
      }

    else if (!strncmp(line, "body names", 10)) {
      s = strtok (line, " ");
      s = strtok (NULL, " ");
      s = strtok (NULL, " ");
      s = strtok (NULL, " ");

      // first read body information //

      while (s != NULL) { 
        s = strtok (NULL, " ");

        if ((s == NULL) || (*s == '}')) {
          break;
          }

        str = s;
        body_names.push_back(s);
        }

      /*
      fprintf (stderr, ">>> bodies = ");
      for (int i = 0; i < body_names.size(); i++) {
        fprintf (stderr, "%s ", body_names[i].c_str());
        }
      fprintf (stderr, "\n");
      */
      break;
      }
    }

  // now read state information //

  states.reserve(num_bodies);

  while (fscanf(fp, "%f\n", &time) == 1) { 
    //fprintf (stderr, ">>> time = %f \n", time);
    for (int i = 0; i < num_bodies; i++) {
      fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f  %f %f %f \n", 
             &disp[0], &disp[1], &disp[2], 
             &mat(0,0), &mat(0,1), &mat(0,2),
             &mat(1,0), &mat(1,1), &mat(1,2),
             &mat(2,0), &mat(2,1), &mat(2,2),
             &vel[0], &vel[1], &vel[2], 
             &avel[0], &avel[1], &avel[2]);

      state.displacement = disp;
      state.rotation_matrix = mat;
      state.velocity = vel;
      state.angular_velocity = avel;
      states[i] = state;
      }
    }

  // set state //

  if (!solver) {
    return;
    }

  solver->setState(body_names, states);
  }

////////////////////////////////////////////////////////////////
//              r e s u l t s   i / o                        //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              writeResults              ==========*
//*============================================================*
// write simulation results to a file. 

void
PmRigidSimulation::writeResults()
  {
  /*
  fprintf (stderr, ">>>>>> PmRigidSimulation::writeResults \n"); 
  fprintf (stderr, ">>> writeResults  last_written %d  num_state_upd %d \n", 
           last_written, num_state_upd);
  fprintf (stderr, ">>>> time size = %d \n", time.size()); 

  for (int i = 0; i < time.size(); i++) {
    fprintf (stderr, " %f \n", time[i]); 
    }
  */

  int state_end_num, state_start_num;
  int num_res_bodies;
  bool verbose = pmSystem.getCmdVerbose();

  if (write_current) {
    //fprintf (stderr, ">>> write current \n");
    state_start_num = time.size()-1;
    state_end_num   = state_start_num+1;
    }
  else {
    state_start_num = last_written;
    state_end_num   = num_state_upd;
    }

  /*
  fprintf (stderr, ">>> state_start_num = %d  \n", state_start_num);
  fprintf (stderr, ">>> state_end_num = %d  \n", state_end_num);
  */

  // write state data //

  if (this->state_fp) {
    vector<PmInteraction*> interactions;
    vector<PmRestraint*> restraints;
    string intr_name, bname;
    PmRigidSimResults result;
    int nb = bodies.size();
    num_res_bodies = 0;

    // write header //

    if (this->write_init) {
      for (int i = 0; i < nb; i++) {
        PmBody *body = bodies[i];
        if (body->hasPhysicalObj()) {
          num_res_bodies += 1;
          }
        }

      fprintf (this->state_fp, "# Protein Mechanica simulation state file \n"); 
      fprintf (this->state_fp, "# time displacement(3) rotation(9) vel(3) avel(3) \n"); 
      fprintf (this->state_fp, "number of bodies %d \n", num_res_bodies); 
      fprintf (this->state_fp, "body names = { "); 

      for (int i = 0; i < nb; i++) {
        PmBody *body = bodies[i];

        if (body->hasPhysicalObj()) {
          body->getName(bname);
          fprintf (this->state_fp, " %s ", bname.c_str()); 
          }
        }

      fprintf (this->state_fp, "}\n"); 
      }

    // write state data //

    for (int n = state_start_num; n < state_end_num; n++) {
      fprintf (this->state_fp, "%f \n", time[n]); 

      for (int i = 0; i < nb; i++) {
        PmBody *body = bodies[i];

        if (body->hasPhysicalObj()) {
          body->getSimResult(n, result);
          fprintf (this->state_fp, "   %f %f %f ", 
                   result.disp[0], result.disp[1], result.disp[2]); 
          fprintf (this->state_fp, "   %f %f %f %f %f %f %f %f %f", 
                   result.rot_mat(0,0), result.rot_mat(0,1), result.rot_mat(0,2), 
                   result.rot_mat(1,0), result.rot_mat(1,1), result.rot_mat(1,2), 
                   result.rot_mat(2,0), result.rot_mat(2,1), result.rot_mat(2,2));
          fprintf (this->state_fp, "   %f %f %f %f %f %f \n", 
                   result.vel[0],  result.vel[1],  result.vel[2],
                   result.avel[0], result.avel[1], result.avel[2]); 
          }
        }

      //fprintf (this->state_fp, "\n"); 
      }

    fflush (this->state_fp);
    }

  // write joint data //

  if (this->joints_fp) {
    string jname;
    PmVector3 pos;
    int num_joints = joints.size();

    // write header //

    if (this->write_init) {
      fprintf (this->joints_fp, "# Protein Mechanica simulation joints file \n");
      fprintf (this->joints_fp, "# time position(3)                         \n");
      fprintf (this->joints_fp, "number of joints %d \n", num_joints);
      fprintf (this->joints_fp, "joint names = { ");

      for (int i = 0; i < num_joints; i++) {
        PmJoint *joint = joints[i];
        joint->getName(jname);
        fprintf (this->joints_fp, " %s ", jname.c_str());
        }

      fprintf (this->joints_fp, "}\n");
      }

    fprintf (this->joints_fp, "%f  ", time[state_end_num-1]); 

    for (int i = 0; i < num_joints; i++) {
      PmJoint *joint = joints[i];
      joint->getCurrentPosition (pos);
      fprintf (this->joints_fp, "%f %f %f  ", pos[0], pos[1], pos[2]); 
      }

    fprintf (this->joints_fp, "\n");
    }

  // write potential energy data //

  if (this->energy_fp) { 
    string intr_name, rest_name;
    float t, energy, total_energy, total_pe, total_ke;
    PmPotentialType ptype;
    bool rest_energy;
    this->getInteractions(interactions);
    this->getRestraints(restraints);

    if (this->write_init) {
      fprintf (this->energy_fp, "# Protein Mechanica simulation energy file \n"); 
      fprintf (this->energy_fp, "# time potential_energy                    \n"); 
      fprintf (this->energy_fp, "interaction names = { "); 
   
      for (unsigned int i = 0; i < interactions.size(); i++) {
        interactions[i]->getPotentialType(ptype); 

        if ((ptype != PM_POTENTIAL_CONTACT) || 
            ((ptype == PM_POTENTIAL_CONTACT) && write_contact_energy)) {
          interactions[i]->getName(intr_name);
          fprintf (this->energy_fp, " %s ", intr_name.c_str());
          }
        }

      fprintf (this->energy_fp, "}\n");

      //===== write restraint energy =====//

      rest_energy = false;

      for (unsigned int i = 0; i < restraints.size(); i++) {
        if (restraints[i]->hasEnergy()) { 
          rest_energy = true;
          break; 
          }
        }

      if (rest_energy) { 
        fprintf (this->energy_fp, "restraint names = { "); 

        for (unsigned int i = 0; i < restraints.size(); i++) {
          restraints[i]->getName(rest_name); 
          fprintf (this->energy_fp, " %s ", rest_name.c_str());
          }

        fprintf (this->energy_fp, "}\n");
        }

      fflush (this->energy_fp);
      }

    total_energy = 0.0;
    bool verbose = pmSystem.getCmdVerbose();

    for (int n = state_start_num; n < state_end_num; n++) {
      total_pe = 0.0;
      fprintf (this->energy_fp, "%f ", time[n]); 

      for (unsigned int i = 0; i < interactions.size(); i++) {
        interactions[i]->getPotentialType(ptype);

        if ((ptype != PM_POTENTIAL_CONTACT) ||
            ((ptype == PM_POTENTIAL_CONTACT) && write_contact_energy)) {

          if (interactions[i]->getEnergy(n, t, energy)) {
            fprintf (this->energy_fp, "%f ", energy); 
            total_pe += energy;
            }
          }
        }

      if (rest_energy) { 
        for (unsigned int i = 0; i < restraints.size(); i++) {
          restraints[i]->getEnergy(energy); 
          fprintf (this->energy_fp, " %f ", energy); 
          total_pe += energy;
          }
        }

      for (unsigned int i = 0; i < joints.size(); i++) {
        joints[i]->getEnergy(n, energy);
        total_pe += energy;
        }

      fprintf (this->energy_fp, "%f ", total_pe); 
      fprintf (this->energy_fp, "\n");

      total_ke = 0.0;

      for (unsigned int i = 0; i < bodies.size(); i++) {
        if (bodies[i]->getKineticEnergy(n, energy)) {
          total_ke += energy;
          }
        }

      if (verbose) {
        fprintf (stderr, " time = %f  pe = %g  ke = %g  total = %g \n", time[n], 
                 total_pe, total_ke, total_pe + total_ke); 
        }
      }

    fflush (this->energy_fp);
    }

  // write kinetic energy data //

  if (this->kenergy_fp) { 
    string bname;
    float energy, total_ke;
    int nb = bodies.size();

    if (this->write_init) {
      fprintf (this->kenergy_fp, "# Protein Mechanica simulation energy file \n"); 
      fprintf (this->kenergy_fp, "# time kinetic_energy                    \n"); 
      fprintf (this->kenergy_fp, "body names = { ");

      for (int i = 0; i < nb; i++) {
        PmBody *body = bodies[i];

        if (body->hasPhysicalObj()) {
          body->getName(bname);
          fprintf (this->kenergy_fp, " %s ", bname.c_str());
          }
        }

      fprintf (this->kenergy_fp, "}\n");
      }

    for (int n = state_start_num; n < state_end_num; n++) {
      fprintf (this->kenergy_fp, "%f ", time[n]); 
      total_ke = 0.0;

      for (unsigned int i = 0; i < bodies.size(); i++) {
        if (bodies[i]->hasPhysicalObj()) {
          bodies[i]->getKineticEnergy(n, energy);
          fprintf (this->kenergy_fp, "%f ", energy); 
          total_ke += energy;
          }
        }

      fprintf (this->kenergy_fp, "%f ", total_ke); 
      fprintf (this->kenergy_fp, "\n");

      if (verbose) {
        fprintf (stderr, " total kinetic energy %f \n", total_ke); 
        }
      }

    fflush (this->kenergy_fp);
    }

  // write strain data //

  if (this->strain_fp) {
    string intr_name;
    float t, strain;
    PmPotentialType ptype;

    if (this->write_init) {
      fprintf (this->strain_fp, "# Protein Mechanica simulation strain file \n");
      fprintf (this->strain_fp, "# time strain                              \n");
      this->getInteractions(interactions);
      fprintf (this->strain_fp, "interaction names = { ");

      for (unsigned int i = 0; i < interactions.size(); i++) {
        interactions[i]->getPotentialType(ptype);

        if (ptype == PM_POTENTIAL_SPRING) {
          interactions[i]->getName(intr_name);
          fprintf (this->strain_fp, " %s ", intr_name.c_str());
          }
        }

      fprintf (this->strain_fp, "}\n");
      }

    for (int n = state_start_num; n < state_end_num; n++) {
      fprintf (this->strain_fp, "%f ", time[n]);
      this->getInteractions(interactions);

      for (unsigned int i = 0; i < interactions.size(); i++) {
        interactions[i]->getPotentialType(ptype);

        if (ptype == PM_POTENTIAL_SPRING) {
          interactions[i]->getStrain(n, t, strain);
          fprintf (this->strain_fp, "%f ", strain);
          }
        }

      fprintf (this->strain_fp, "\n");
      }

    fflush (this->strain_fp);
    }

  // write joint data //

  for (unsigned int i = 0; i < joints.size(); i++) {
    joints[i]->writeResults(time, state_start_num, state_end_num);
    }

  // write trace data //

  for (unsigned int i = 0; i < traces.size(); i++) {
    traces[i]->writeResults(time, state_start_num, state_end_num);
    }

  // write motor data //

  for (unsigned int i = 0; i < motors.size(); i++) {
    motors[i]->writeResults(time, state_start_num, state_end_num);
    }

  // write measurement data //

  for (unsigned int i = 0; i < measurements.size(); i++) {
    measurements[i]->writeResults(time, state_start_num, state_end_num);
    }

  // write domain data //

  if (this->pdb_db) { 
    PmPhysicalObj *pobj;
    vector<PmMolecule*> domains; 
    vector<PmBody*> domain_bodies; 
    vector<PmXform> xforms;
    PmXform xform;
    PmVector3 com;
    PmRigidSimResults result;
    PmDbPdbInterface *pdb = dynamic_cast<PmDbPdbInterface*>(this->pdb_db);

    if (this->write_init) {
      pdb->writeHeader(); 
      }

    for (unsigned int i = 0; i < bodies.size(); i++) {
      if (bodies[i]->isStatic()) {
        continue;
        }

      bodies[i]->getPhysicalObj(&pobj);

      if (pobj && ((pobj->getObjectType() == PM_OBJECT_DOMAIN) || 
                   (pobj->getObjectType() == PM_OBJECT_MOLECULE))) {
        domains.push_back(dynamic_cast<PmMolecule*>(pobj));
        domain_bodies.push_back(bodies[i]);
        }
      }

    for (int n = state_start_num; n < state_end_num; n++) {
      //fprintf (stderr, "   >>> write domains \n");
      xforms.clear();

      for (unsigned int i = 0; i < domain_bodies.size(); i++) {
        domain_bodies[i]->getSimResult(n, result);
        domain_bodies[i]->getCenterOfMass(com);
        xform.center = com;
        xform.matrix = result.rot_mat;
        xform.translation = result.disp;
        xforms.push_back(xform);
        }

      pdb->writeDynDomainAtoms(n, this->write_init, domains, xforms);
      }
    }

  last_written = num_state_upd;
  this->write_init = false;
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              getBodyId                 ==========*
//*============================================================*

bool
PmRigidSimulation::getBodyId (const string str, int& id)
  {
  string bname;
  int n = bodies.size();

  for (int i = 0; i < n; i++) {
    PmBody *bdy = bodies[i];
    bdy->getName(bname);

    if (bname == str) {
      id = i;
      return true;
      }
    }

  id = -1;
  return false;
  }




}

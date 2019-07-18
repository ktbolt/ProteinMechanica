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
//* pm:                     s y s t e m                        *
//*============================================================*

#include "pm/pm.h"
#include "bsim.h"
#include "db.h"
#include "force.h"
#include "geom.h"
#include "graphics.h"
#include "joint.h"
#include "model.h"
#include "mol.h"
#include "motor.h"
#include "msr.h"
#include "sim.h"
#include "solid.h"
#include "trace.h"

namespace ProteinMechanica {

//void pm_CmdScriptProc (const char *script);
void pm_CmdContinueProc();
void pm_CmdAddVariable(const string name, const string value);

PmSystem pmSystem;

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSystem::PmSystem() 
    {
    //fprintf (stderr, ">>>>>> PmSystem: ctor \n");
    cmd_echo = true;
    cmd_verbose = true;
    num_cmd_fp = 0;
    cmd_wait = false;
    cmd_error = false;

    use_graphics = false;
    graphics = NULL; 
    update_graphics = false; 
    current_molecule = NULL;
    current_domain = NULL;
    current_surface = NULL;

    // initialize units
    units.mass.scale = 1.0;

    replay_simulations = false;
    extent_set = false;
    silent = false;
    }

//*============================================================*
//*==========              init                      ==========*
//*============================================================*

void 
PmSystem::init(const bool silent)
  {
  //fprintf (stderr, ">>>>>> PmSystem: init \n");

#ifdef USE_GRAPHICS 
  use_graphics = true;
#else
  use_graphics = false;
#endif

  graphics = new PmGraphicsInterface (use_graphics);
  graphics->init (); 

  if (!silent) {
    fprintf (stderr, "====== Protein Mechanica %s %s ======\n", __DATE__, __TIME__);

    if (!use_graphics) {
      fprintf (stderr, ">>> no graphics \n");
      }
    }
  else {
    cmd_echo = false;
    cmd_verbose = false;
    }

  // set default units // 

  units.base = PM_UNIT_BASE_LMT;

  units.mass.type = PM_UNIT_MASS;
  units.mass.name = "amu";
  units.mass.scale = 1.0;
  units.mass.unit.value = PM_UNITS_AMU;
  units.mass.unit.name = "kg";

  units.length.type = PM_UNIT_LENGTH;
  units.length.name = "nanometer";
  units.length.scale = 1.0;
  units.length.unit.value = PM_UNITS_NANO;
  units.length.unit.name = "m";

  units.temperature.type = PM_UNIT_TEMPERATURE;
  units.temperature.name = "Kelvin";
  units.temperature.scale = 1.0;
  units.temperature.unit.value = 298.15;
  units.temperature.unit.name = "K";

  units.energy.type = PM_UNIT_ENERGY;
  units.energy.name = "Joule";
  units.energy.scale = 1.0;
  units.energy.unit.value = PM_CONS_KB * units.temperature.unit.value;
  units.energy.unit.name = "J";

  units.time.type = PM_UNIT_TIME;
  units.time.name = "picosecond";
  units.time.scale = 1.0;
  units.time.unit.value = PM_UNITS_PICO;
  units.time.unit.name = "s";

  deriveUnits();
  }

//*============================================================*
//*==========              setExtent                 ==========*
//*============================================================*

void 
PmSystem::setExtent (PmExtent& extent)
  {
  if (use_graphics) {
    graphics->setExtent (extent); 
    }

  extent_set = true;
  this->extent = extent;
  }

//*============================================================*
//*==========              updateExtent              ==========*
//*============================================================*

void
PmSystem::updateExtent(PmExtent& extent)
  {
  /*
  fprintf (stderr, "\n>>>>>> PmSystem::updateExtent \n");
  fprintf (stderr, "       extent_set = %d \n", extent_set);
  fprintf (stderr, "       x = %f %f \n", extent.min[0], extent.max[0]);
  fprintf (stderr, "       y = %f %f \n", extent.min[1], extent.max[1]);
  fprintf (stderr, "       z = %f %f \n", extent.min[2], extent.max[2]);
  */

  if (!extent_set) {
    setExtent(extent);
    }

  else {
    this->extent.update(extent);

    if (use_graphics) {
      graphics->setExtent (this->extent);
      }
    }
  }

//*============================================================*
//*==========              addBody                   ==========*
//*============================================================*

void
PmSystem::addBody (PmBody *body) {
  bodies.push_back (body);
  //current_body = body;
  }

//*============================================================*
//*==========              getBody                   ==========*
//*============================================================*

void
PmSystem::getBody (const string name, PmBody **body)
  {
  string str;
  *body = NULL;

  for (unsigned int i = 0; i < bodies.size(); i++) {
    PmBody *bdy = bodies[i];
    bdy->getName (str);

    if (str == name) {
      *body = bdy;
      return;
      }
    }
  }

//*============================================================*
//*==========              getBodies                 ==========*
//*============================================================*

void 
PmSystem::getBodies (vector<PmBody*>& bodies)
  {
  bodies = this->bodies;
  }

//*============================================================*
//*==========              potentials                ==========*
//*============================================================*

void 
PmSystem::addPotential(PmPotential *pot) 
  {
  potentials.push_back (pot);
  }

void 
PmSystem::getPotential(const string name, PmPotential **p_pot)
  {
  string str;
  *p_pot = NULL;

  for (unsigned int i = 0; i < potentials.size(); i++) {
    PmPotential *pot = potentials[i];
    pot->getName (str);

    if (str == name) {
      *p_pot = pot;
      return;
      }
    }
  }

void
PmSystem::addPotentialGeometry(PmPotentialGeom *pot)
  {
  potential_geometries.push_back (pot);
  }

void
PmSystem::getPotentialGeometry(const string name, PmPotentialGeom **p_pot)
  {
  string str;
  *p_pot = NULL;

  for (unsigned int i = 0; i < potential_geometries.size(); i++) {
    PmPotentialGeom *pot = potential_geometries[i];
    pot->getName (str);

    if (str == name) {
      *p_pot = pot;
      return;
      }
    }
  }

//*============================================================*
//*==========              get/set WorkingDir        ==========*
//*============================================================*

void 
PmSystem::getWorkingDir(string& val) {
  val = this->working_directory;
  }

void 
PmSystem::setWorkingDir(const string val) {
  this->working_directory = val;
  }

//*============================================================*
//*==========              getVariableValue          ==========*
//*============================================================*
// get the value of a system variable.

float 
PmSystem::getVariableValue(const string name)
  {
  //fprintf (stderr, "\n>>>>>> PmSystem::getVariableValue name=%s \n", name.c_str()); 
  float value;
  char c;
  vector<string> names(10);
  int n = 0;

  for (unsigned int i = 0; i < name.size(); i++) {
    c = name[i];

    if (c == ':') {
      n += 1;
      }
    else {
      names[n].push_back(c);
      }
    }

  //===== process potential energy =====//

  if (name.find("potential_energy") != string::npos) {
    float energy;

    // total system energy //

    if (n == 0) {
      current_simulation->getEnergy(energy);
      //fprintf (stderr, ">>> total potential energy=%f \n", energy); 
      value = energy;
      }
    }

  //===== process time step =====//

  else if (name.find("step") != string::npos) {
    int step;
    current_simulation->getNumSteps(step);
    //fprintf (stderr, ">>> step=%d \n", step); 
    value = step;
    }

  return value;
  }

//*============================================================*
//*==========              Database                  ==========*
//*============================================================*

void
PmSystem::addDatabase(PmDbInterface *db) {
  databases.push_back (db);
  }

void
PmSystem::getDatabase(const string name, PmDbInterface **p_db)
  {
  string str;
  *p_db = NULL;

  for (unsigned int i = 0; i < databases.size(); i++) {
    PmDbInterface *db = databases[i];
    db->getName (str);

    if (str == name) {
      *p_db = db;
      return;
      }
    }
  }

//*============================================================*
//*==========              Domain                    ==========*
//*============================================================*

void
PmSystem::addDomain (PmMolecule *domain) {
  domains.push_back (domain);
  current_domain = domain;
  }

void
PmSystem::getDomain (const string name, PmMolecule **domain)
  {
  string str;
  *domain = NULL;

  for (unsigned int i = 0; i < domains.size(); i++) {
    PmMolecule *mol = domains[i];
    mol->getName (str);

    if (str == name) {
      *domain = mol;
      return;
      }
    }
  }

void 
PmSystem::getDomains(vector<PmMolecule*>& domains) {
  domains = this->domains;
  }

//*============================================================*
//*==========              addForce                  ==========*
//*============================================================*

void 
PmSystem::addForce(PmForce *force) {
  //fprintf (stderr, "\n###### add force=%s \n", force);
  forces.push_back(force);
  }

//*============================================================*
//*==========              getForce                  ==========*
//*============================================================*

void 
PmSystem::getForce(const string name, PmForce **force)
  {
  string str;
  *force = NULL;

  for (unsigned int i = 0; i < forces.size(); i++) {
    PmForce *f = forces[i];
    f->getName (str);

    if (str == name) {
      *force = f;
      return;
      }
    }
  }

//*============================================================*
//*==========              add/getGeometry           ==========*
//*============================================================*

void 
PmSystem::addGeometry(PmGeometry *geom)
  {
  geometries.push_back(geom);
  }

void 
PmSystem::getGeometry(const string name, PmGeometry **geom)
  {
  string str;
  *geom = NULL;

  for (unsigned int i = 0; i < geometries.size(); i++) {
    PmGeometry *g = geometries[i];
    g->getName (str);

    if (str == name) {
      *geom = g;
      return;
      }
    }
  }

//*============================================================*
//*==========              getGraphicsContext        ==========*
//*============================================================*

void
PmSystem::getGraphicsContext (void **context) {
  if (!use_graphics) return;
  graphics->getContext (context);
  }

//*============================================================*
//*==========         getGraphicsInterace            ==========*
//*============================================================*

PmGraphicsInterface *
PmSystem::getGraphicsInterace() {
  return graphics;
  }

//*============================================================*
//*==========         updateGraphics                 ==========*
//*============================================================*
// set update graphics flag.

void
PmSystem::updateGraphics(const bool flag) {
  update_graphics = flag;
  }

//*============================================================*
//*==========         updateGraphics                 ==========*
//*============================================================*
// update graphics, if needed.

void
PmSystem::updateGraphics() 
  {
  if (update_graphics) {
    //fprintf (stderr, ">>>>>> PmSystem::updateGraphics \n");
    graphics->update();
    }

  update_graphics = false;
  }

//*============================================================*
//*==========         forceUpdateGraphics            ==========*
//*============================================================*
// force the update of graphics.

void
PmSystem::forceUpdateGraphics() {
  graphics->update();
  }

//*============================================================*
//*==========              useGraphics               ==========*
//*============================================================*
// check if graphics is enabled.

int
PmSystem::useGraphics()  
  {
  return (use_graphics);
  }

//*============================================================*
//*==========              addGrid                   ==========*
//*============================================================*

void
PmSystem::addGrid (PmGrid *grid) {
  grids.push_back(grid);
  }

//*============================================================*
//*==========              getGrid                   ==========*
//*============================================================*

void
PmSystem::getGrid (const string name, PmGrid **grid)
  {
  string str;
  *grid = NULL;

  for (unsigned int i = 0; i < grids.size(); i++) {
    PmGrid *g = grids[i];
    g->getName (str);

    if (str == name) {
      *grid = g;
      return;
      }
    }
  }

//*============================================================*
//*==========              getCmdFilePtr             ==========*
//*============================================================*
// get cmd file ptr.                

FILE* 
PmSystem::getCmdFilePtr() 
  { 
  if (num_cmd_fp == 0) {
    return NULL;
    }

  num_cmd_fp -= 1;
  return cmd_fp[num_cmd_fp]; 
  }

//*============================================================*
//*==========              setCmdFilePtr             ==========*
//*============================================================*
// set cmd file ptr.

void 
PmSystem::setCmdFilePtr(FILE *fp) 
  { 

  if (num_cmd_fp >= (int)cmd_fp.size()) {
    cmd_fp.push_back(fp); 
    }
  else {
    cmd_fp[num_cmd_fp] = fp; 
    }

  num_cmd_fp += 1;
  }

//*============================================================*
//*==========              addCmdVariable            ==========*
//*============================================================*
// add a command variable.

void 
PmSystem::addCmdVariable(const string name, const string value) {
  pm_CmdAddVariable(name, value);
  }

//*============================================================*
//*==========              setCmdWait                ==========*
//*============================================================*
// set to wait for a command.

void 
PmSystem::setCmdWait(const bool val) 
  { 
  cmd_wait = val; 
  int n = num_cmd_fp;

  if (!val) {
    for (int i = 0; i < n; i++) {
      pm_CmdContinueProc();
      }
    }
  }

//*============================================================*
//*==========              procInterrupt             ==========*
//*============================================================*
// process an interrupt.               

void
PmSystem::procInterrupt()
  {

  // halt all simulations //

  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->halt();
    }
  }

//*============================================================*
//*==========              addJoint                  ==========*
//*============================================================*

void
PmSystem::addJoint (PmJoint *joint) {
  joints.push_back (joint);
  }

//*============================================================*
//*==========              getJoint                  ==========*
//*============================================================*

void
PmSystem::getJoint (const string name, PmJoint **joint)
  {
  string str;
  *joint = NULL;

  for (unsigned int i = 0; i < joints.size(); i++) {
    PmJoint *jnt = joints[i];
    jnt->getName (str);

    if (str == name) {
      *joint = jnt;
      return;
      }
    }
  }

void 
PmSystem::getJoints (vector<PmJoint*>& joints)
  {
  joints = this->joints;
  }

//*============================================================*
//*==========              procEvents                ==========*
//*============================================================*
// event processing loop.                  

void 
PmSystem::procEvents (const char *script) 
  {

  // process script //

  /*
  if (script) {
    pm_CmdScriptProc (script);
    }
  */

  graphics->procEvents (script); 

  if (!use_graphics) {
    int num_asim = 0;

    while (1) {
      this->stepSimulations();
      num_asim = 0;

      for (unsigned int i = 0; i < simulations.size(); i++) {
        PmSimulation *sim = simulations[i];
        if (sim->active()) num_asim += 1;
        }

      if (num_asim == 0) {
        break;
        }
      }
    } 
  }

//*============================================================*
//*==========              addMeasurement            ==========*
//*============================================================*

void 
PmSystem::addMeasurement(PmMeasurement *msr) {
  measurements.push_back (msr);
  }

void 
PmSystem::getMeasurement(const string name, PmMeasurement **msr) {
  string str;
  *msr = NULL;

  for (unsigned int i = 0; i < measurements.size(); i++) {
    measurements[i]->getName (str);

    if (str == name) {
      *msr = measurements[i];
      return;
      }
    }
  }

//*============================================================*
//*==========              addModel                  ==========*
//*============================================================*

void
PmSystem::addModel (PmModel *model) {
  models.push_back (model);
  current_model = model;
  }

//*============================================================*
//*==========              getModel                  ==========*
//*============================================================*

void
PmSystem::getModel (const string name, PmModel **rmodel)
  {
  string str;
  *rmodel = NULL;

  for (unsigned int i = 0; i < models.size(); i++) {
    PmModel *model = models[i];
    model->getName (str);

    if (str == name) {
      *rmodel = model;
      return;
      }
    }
  }

//*============================================================*
//*==========              addMolecule               ==========*
//*============================================================*

void 
PmSystem::addMolecule (PmMolecule *mol) {
  molecules.push_back (mol);
  current_molecule = mol;
  }

//*============================================================*
//*==========              getMolecule               ==========*
//*============================================================*

void 
PmSystem::getMolecule (const string name, PmMolecule **molecule) 
  {
  //fprintf (stderr, ">>>>>> PmSystem::getMolecule  name[%s] \n", name.c_str());
  string str;
  *molecule = NULL;

  for (unsigned int i = 0; i < molecules.size(); i++) {
    PmMolecule *mol = molecules[i];
    mol->getName (str);
    //fprintf (stderr, ">>> str [%s] \n", str.c_str());

    if (str == name) {
      *molecule = mol;
      return;
      }
    }
  }

//*============================================================*
//*==========              addMotor                  ==========*
//*============================================================*

void
PmSystem::addMotor(PmMotor *motor) {
  motors.push_back(motor);
  }

//*============================================================*
//*==========              getMotor                  ==========*
//*============================================================*

void
PmSystem::getMotor(const string name, PmMotor **motor)
  {
  string str;
  *motor = NULL;

  for (unsigned int i = 0; i < motors.size(); i++) {
    PmMotor *obj = motors[i];
    obj->getName(str);

    if (str == name) {
      *motor = obj;
      return;
      }
    }
  }

//*============================================================*
//*==========              addParticle               ==========*
//*============================================================*

void
PmSystem::addParticle (PmParticle *part) {
  particles.push_back(part);
  }

//*============================================================*
//*==========              getParticle               ==========*
//*============================================================*

void
PmSystem::getParticle (const string name, PmParticle **part)
  {
  string str;
  *part = NULL;

  for (unsigned int i = 0; i < particles.size(); i++) {
    PmParticle *p = particles[i];
    p->getName (str);
    fprintf (stderr, ">>> str = \"%s\" \n", str.c_str());

    if (str == name) {
      *part = p;
      return;
      }
    }
  }

//*============================================================*
//*==========              procQuery                 ==========*
//*============================================================*

void 
PmSystem::procQuery(PmQuery& query)
  {
  //fprintf (stderr, ">>>>>> PmSystem::procQuery  \"%s\" \n", query.name.c_str());
  string obj_type, obj_name;
  unsigned int n1, n2;

  n1 = query.name.find("[");
  n2 = query.name.find("]");

  if ((n1 == string::npos) || (n2 == string::npos)) {
    return;
    }

  obj_type.assign(query.name, 0, n1);
  obj_name.assign(query.name, n1+1, n2-n1-1);
  //fprintf (stderr, "   >>> object type \"%s\" \n", obj_type.c_str());
  //fprintf (stderr, "   >>> object name \"%s\" \n", obj_name.c_str());

  if (obj_type == "domain") {
    PmMolecule::procQuery(query);
    }
  else if (obj_type == "potential") {
    PmPotential::procQuery(query);
    }
  else if (obj_type == "force") {
    PmForce::procQuery(query);
    }
  }

//*============================================================*
//*==========              procQuit                  ==========*
//*============================================================*

void
PmSystem::procQuit()
  {
  exit(0);
  }
//*============================================================*
//*==========              addSimulation             ==========*
//*============================================================*

void
PmSystem::addBindSimulation (PmBindingSimulation *sim){
  binding_simulations.push_back (sim);
  //current_simulation = sim;
  }

void
PmSystem::getBindSimulation (const string name, PmBindingSimulation **simulation)
  {
  string str;
  *simulation = NULL;

  for (unsigned int i = 0; i < binding_simulations.size(); i++) {
    PmBindingSimulation *sim = binding_simulations[i];
    sim->getName(str);

    if (str == name) {
      *simulation = sim;
      return;
      }
    }
  }
  
//*============================================================*
//*==========              addSimulation             ==========*
//*============================================================*

void
PmSystem::addSimulation (PmSimulation *sim){
  simulations.push_back (sim);
  current_simulation = sim;
  string str;
  sim->getName(str);
  //fprintf (stderr, ">>>>>> PmSystem:: addSimulation  [%s] \n", str.c_str());
  }

//*============================================================*
//*==========              getSimulation             ==========*
//*============================================================*

void
PmSystem::getSimulation (const string name, PmSimulation **simulation)
  {
  string str;
  *simulation = NULL;

  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->getName(str);

    if (str == name) {
      *simulation = sim;
      return;
      }
    }
  }

void
PmSystem::getSimulations (vector<PmSimulation*> simulations)
  {
  simulations = this->simulations;
  }

//*============================================================*
//*==========        setSimulationsNumInc            ==========*
//*============================================================*
// set the number of increments each simulation will perform.

void 
PmSystem::setSimulationsNumInc(const int num_inc)
  {
  //fprintf (stderr, ">>>>>> PmSystem::setSimulationsNumInc [%d] \n", num_inc); 
  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->setNumInc(num_inc);
    }
  }

//*============================================================*
//*==========        setSimulationsStateSaveFreq     ==========*
//*============================================================*
// set the frequency to save simulation state.

void
PmSystem::setSimulationsStateSaveFreq(const int save_freq)
  {
  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->setStateSaveFreq(save_freq);
    }
  }

//*============================================================*
//*==========              replaySimulations         ==========*
//*============================================================*

void
PmSystem::setReplaySimulations(const bool flag, const int start, const int end,
                               const int inc)
  {
  //fprintf (stderr, ">>>>>> PmSystem::setReplaySimulations \n"); 
  //fprintf (stderr, ">>> replay_simulations = %d \n", flag); 
  replay_simulations = flag;

  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->setReplay(flag, start, end, inc);
    }
  }

//*============================================================*
//*==========              replaySimulations         ==========*
//*============================================================*

void
PmSystem::replaySimulations()
  {

  if (!replay_simulations) { 
    return;
    }

  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->replay();
    }
  }

//*============================================================*
//*==========              stepSimulations           ==========*
//*============================================================*

void
PmSystem::stepSimulations()
  {
  //fprintf (stderr, ">>>>>> PmSystem::stepSimulations \n"); 
  for (unsigned int i = 0; i < simulations.size(); i++) {
    PmSimulation *sim = simulations[i];
    sim->step(1);
    }
  }

//*============================================================*
//*==========              addSolid                  ==========*
//*============================================================*

void 
PmSystem::addSolid(PmSolid *solid) {
  solids.push_back (solid);
  }

//*============================================================*
//*==========              getSolid                  ==========*
//*============================================================*

void 
PmSystem::getSolid(const string name, PmSolid **solid)
  {
  string str;
  *solid = NULL;

  for (unsigned int i = 0; i < solids.size(); i++) {
    PmSolid *obj = solids[i];
    obj->getName(str);

    if (str == name) {
      *solid = obj;
      return;
      }
    }
  }

//*============================================================*
//*==========              addSurface                ==========*
//*============================================================*

void
PmSystem::addSurface (PmSurface *surf) {
  surfaces.push_back (surf);
  current_surface = surf;
  }

//*============================================================*
//*==========              getSurface                ==========*
//*============================================================*

void
PmSystem::getSurface (const string name, PmSurface **surface)
  {
  string str;
  *surface = NULL;

  for (unsigned int i = 0; i < surfaces.size(); i++) {
    PmSurface *surf = surfaces[i];
    surf->getName(str);

    if (str == name) {
      *surface = surf;
      return;
      }
    }
  }

//*============================================================*
//*==========              addTrace                  ==========*
//*============================================================*

void
PmSystem::addTrace(PmTrace *trace) {
  traces.push_back (trace);
  }

//*============================================================*
//*==========              getTrace                  ==========*
//*============================================================*

void
PmSystem::getTrace(const string name, PmTrace **trace)
  {
  string str;
  *trace = NULL;

  for (unsigned int i = 0; i < traces.size(); i++) {
    PmTrace *trc = traces[i];
    trc->getName(str);

    if (str == name) {
      *trace = trc;
      return;
      }
    }
  }

//*============================================================*
//*==========              convUnitType              ==========*
//*============================================================*

void
PmSystem::convUnitType(const string str, PmUnitType& type) 
  {
  if (str == "energy") {
    type = PM_UNIT_ENERGY;
    }
  else if (str == "length") {
    type = PM_UNIT_LENGTH;
    }
  else if (str == "mass") {
    type = PM_UNIT_MASS;
    }
  else if (str == "temperature") {
    type = PM_UNIT_TEMPERATURE;
    }
  else if (str == "time") {
    type = PM_UNIT_TIME;
    }
  else { 
    type = PM_UNIT_UNKNOWN;
    }
  }

//*============================================================*
//*==========              deriveUnits               ==========*
//*============================================================*

void
PmSystem::deriveUnits()
  {

  float energy_u, mass_u, length_u, time_u;

  float T, tmp;

  float kb_r, kT_r, eta_r;

  mass_u = units.mass.unit.value;
  length_u = units.length.unit.value;
  time_u = units.time.unit.value;

  if (units.base == PM_UNIT_BASE_LMT) {
    tmp = length_u / time_u;
    energy_u = mass_u * tmp * tmp;
    units.energy.unit.value = energy_u;
    }

  T = units.temperature.unit.value;
  kb_r = PM_CONS_KB * time_u*time_u / (mass_u * length_u*length_u);
  kT_r = kb_r * T;
  eta_r = PM_CONS_ETA * (length_u*time_u) / mass_u;

  units.reduced.eta = eta_r;
  units.reduced.kb = kb_r;
  units.reduced.kT = kT_r;
  }

//*============================================================*
//*==========              printUnits                ==========*
//*============================================================*

void
PmSystem::printUnits()
  {
  if (!cmd_verbose) {
    return;
    }

  fprintf (stderr, "\n---------------- energy ---------------- \n");
  fprintf (stderr, ">>>>>> name  [%s] \n", units.energy.name.c_str());
  fprintf (stderr, ">>>>>> scale [%g] \n", units.energy.scale);
  fprintf (stderr, ">>>>>> value [%g %s] \n", units.energy.unit.value,
           units.energy.unit.name.c_str());

  fprintf (stderr, "\n---------------- length ---------------- \n");
  fprintf (stderr, ">>>>>> name  [%s] \n", units.length.name.c_str());
  fprintf (stderr, ">>>>>> scale [%g] \n", units.length.scale);
  fprintf (stderr, ">>>>>> value [%g %s] \n", units.length.unit.value,
           units.length.unit.name.c_str());

  fprintf (stderr, "\n---------------- mass ---------------- \n");
  fprintf (stderr, ">>>>>> name  [%s] \n", units.mass.name.c_str());
  fprintf (stderr, ">>>>>> scale [%g] \n", units.mass.scale);
  fprintf (stderr, ">>>>>> value [%g %s] \n", units.mass.unit.value,
           units.mass.unit.name.c_str());

  fprintf (stderr, "\n---------------- temperature ---------------- \n");
  fprintf (stderr, ">>>>>> name  [%s] \n", units.temperature.name.c_str());
  fprintf (stderr, ">>>>>> scale [%g] \n", units.temperature.scale);
  fprintf (stderr, ">>>>>> value [%g %s] \n", units.temperature.unit.value,
           units.temperature.unit.name.c_str());

  fprintf (stderr, "\n---------------- time ---------------- \n");
  fprintf (stderr, ">>>>>> name  [%s] \n", units.time.name.c_str());
  fprintf (stderr, ">>>>>> scale [%g] \n", units.time.scale);
  fprintf (stderr, ">>>>>> value [%g %s] \n", units.time.unit.value,
           units.time.unit.name.c_str());

  fprintf (stderr, "\n---------------- reduced units ---------------- \n");
  fprintf (stderr, ">>>>>> eta [%g] \n", units.reduced.eta);
  fprintf (stderr, ">>>>>> kb  [%g] \n", units.reduced.kb);
  fprintf (stderr, ">>>>>> kT  [%g] \n", units.reduced.kT);
  }

//*============================================================*
//*==========              setUnitsScale             ==========*
//*============================================================*

void
PmSystem::setUnitsScale (PmUnitType type, const string name, float scale)
  {
  switch (type) {
    case PM_UNIT_MASS:
      units.mass.name = name;
      units.mass.scale = scale;
      units.mass.unit.value *= scale;
    break;

    case PM_UNIT_LENGTH:
      units.length.name = name;
      units.length.scale = scale;
      units.length.unit.value *= scale;
    break;

    default:
    break;
    }
  }

//*============================================================*
//*==========              pm_SystemIdle             ==========*
//*============================================================*
// idle function callback from a graphics system. when the
// graphics system is idle this function will be called.

void
pm_SystemIdle() 
  {

  pmSystem.updateGraphics(false);

  // step simulations //

  pmSystem.stepSimulations();

  // replay simulations //

  pmSystem.replaySimulations();

  // update graphics //

  pmSystem.updateGraphics();
  }

}


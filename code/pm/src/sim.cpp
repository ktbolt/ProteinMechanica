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
// * sim:           s i m u l a t i o n                        *
//*============================================================*

#include "sim.h"
#include "rbsim.h"
#include "db.h"

#define ndbg_PmInteraction_compForces
#define ndbg_PmInteraction_PmInteraction_name_list

namespace ProteinMechanica {

vector<PmSimulation*> PmSimulation::simulations;

////////////////////////////////////////////////////////////////
//                     s i m u l a t i o n                   //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSimulation::PmSimulation () 
  {
  //fprintf (stderr, ">>> PmSimulation: ctor \n");
  write_state = false;
  state_fp = NULL;

  write_energy = false;
  write_contact_energy = false;
  energy_fp = NULL;

  write_joints = false;
  joints_fp = NULL;

  write_kenergy = false;
  kenergy_fp = NULL;

  write_strain = false;
  strain_fp = NULL;

  write_domains = false;
  pdb_db = NULL;

  binary_format = true;
  write_init = true;
  last_written = 0;
  num_state_upd = 0;
  write_current = false;
  write_results = false;
  silent_step_ = false;

  while_exp = NULL;
  simulation_active = false;
  }

PmSimulation::PmSimulation (const string name) 
  { 
  //fprintf (stderr, ">>> PmSimulation: ctor [%s] \n", name.c_str());
  this->name = name;
  PmSimulation::simulations.push_back(this);
  }




//*============================================================*
//*==========              getName                   ==========*
//*============================================================*

void 
PmSimulation::getName(string& name) { 
  name = this->name;  
  }

//*============================================================*
//*==========              getType                   ==========*
//*============================================================*

PmSimulationType 
PmSimulation::getType() { 
  return type; 
  }

//*============================================================*
//*==========              convSimulationType        ==========*
//*============================================================*
// convert a string simulation type to a symbol.

void 
PmSimulation::convSimulationType(const string str, PmSimulationType& type)
  {
  type = PM_SIMULATION_UNKNOWN;

  if (str == "rigid") {
    type = PM_SIMULATION_RIGID;
    }
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a simulation object of a particular type.

PmSimulation* 
PmSimulation::create(const string name, const PmSimulationType type)
  {
  switch (type) {
    case PM_SIMULATION_RIGID:
      {
      PmRigidSimulation *rsim = new PmRigidSimulation(name);
      rsim->type = type;
      return rsim;
      }
    break;

    default:
      return NULL;
    }
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of grapical objects related to simulations.

void 
PmSimulation::updateDisplay()
  {
  //fprintf (stderr, ">>> PmSimulation: updateDisplay \n");

  for (int i = 0; i < interactions.size(); i++) {
    interactions[i]->updateDisplay();
    }
  }

//*============================================================*
//*==========              setSilentStep             ==========*
//*============================================================*

void
PmSimulation::setSilentStep(const bool flag) {
  silent_step_ = flag; 
  }

//*============================================================*
//*==========              addInteraction            ==========*
//*============================================================*
// add an interaction.

void 
PmSimulation::addInteraction(PmInteraction *inter)
  {
  interactions.push_back(inter);
  }

//*============================================================*
//*==========              getInteractions           ==========*
//*============================================================*

void 
PmSimulation::getInteractions(vector<PmInteraction*>& intr) { 
  intr = this->interactions; 
  }

//*============================================================*
//*==========              setNumInc                 ==========*
//*============================================================*
// set number of solver steps to increment.

void 
PmSimulation::setNumInc(int n) 
  { 
  num_inc = n; 
  inc_count = n;
  simulation_active = true;
  }

//*============================================================*
//*==========              setNumSteps               ==========*
//*============================================================*
// get the number of solver steps.

void
PmSimulation::getNumSteps(int& n) {
  n = num_steps;
  }

//*============================================================*
//*==========              setReplay                 ==========*
//*============================================================*
// set a simulation to replay it's results.

void
PmSimulation::setReplay(const bool flag, const int start, const int end, const int inc)
  {
  //fprintf (stderr, "\n>>>>>> PmSimulation::setReplay name = %s \n", name.c_str()); 
  replay_active = flag;
  replay_start = start;
  replay_end = end;
  replay_inc = inc;
  }

//*============================================================*
//*==========              addRestraint              ==========*
//*============================================================*
// add a restraint.                                            

void 
PmSimulation::addRestraint(PmRestraint *res)
  {
  restraints.push_back(res);
  }

//*============================================================*
//*==========              getRestraints             ==========*
//*============================================================*

void 
PmSimulation::getRestraints(vector<PmRestraint*>& res) { 
  res = this->restraints; 
  }

//*============================================================*
//*==========              checkWhileExp             ==========*
//*============================================================*

bool
PmSimulation::checkWhileExp()
  {
  string res;

  if (while_exp) {
    return while_exp->evaluate();
    }
  }

//*============================================================*
//*==========              setWhileExp               ==========*
//*============================================================*

void 
PmSimulation::setWhileExp (PmBooleanExpression *bexp)
  {
  while_exp = bexp;
  }

//*============================================================*
//*==========               active                   ==========*
//*============================================================*

bool 
PmSimulation::active() {
  return simulation_active;
  }

//*============================================================*
//*==========         setWriteState                  ==========*
//*============================================================*
// set flag to write state data. 

void 
PmSimulation::setWriteContact(const bool flag) {
  write_contact_energy = flag;
  }

void 
PmSimulation::setWriteDomains(const bool flag) {
  write_domains = flag;
  }

void 
PmSimulation::setWriteJoints(const bool flag) {
  write_joints = flag;
  }

void 
PmSimulation::setWriteState(const bool flag) {
  write_state = flag;
  }

void 
PmSimulation::setWriteEnergy(const bool flag) {
  write_energy = flag;
  }

void 
PmSimulation::setWriteKineticEnergy(const bool flag) {
  write_kenergy = flag;
  }

void 
PmSimulation::setWriteStrain(const bool flag) {
  write_strain = flag;
  }

void
PmSimulation::setWriteCurrent(const bool flag) {
  write_current = flag;
  }

//*============================================================*
//*==========         setWriteResultsParams          ==========*
//*============================================================*
// set the paramaters for writing results.

void 
PmSimulation::setWriteResultsParams (string file_name, bool bformat)
  {
  char fname[100];
  this->write_file_name = file_name;
  this->binary_format = bformat;
  this->write_init = true;

  if (write_state && !this->state_fp) {
    sprintf (fname, "%s_state.pm", file_name.c_str());
    this->state_fp = fopen (fname, "w");
    }

  if (write_energy && !this->energy_fp) {
    sprintf (fname, "%s_energy.pm", file_name.c_str());
    this->energy_fp = fopen (fname, "w");
    }

  if (write_joints && !this->joints_fp) {
    sprintf (fname, "%s_joints.pm", file_name.c_str());
    this->joints_fp = fopen (fname, "w");
    }

  if (write_kenergy && !this->kenergy_fp) {
    sprintf (fname, "%s_kenergy.pm", file_name.c_str());
    this->kenergy_fp = fopen (fname, "w");
    }

  if (write_strain && !this->strain_fp) {
    sprintf (fname, "%s_strain.pm", file_name.c_str());
    this->strain_fp = fopen (fname, "w");
    }

  if (write_domains) {
    PmDbInterfaceSelect db_sel;
    sprintf (fname, "%s_domains.pdb", file_name.c_str());
    //this->domains_fp = fopen (fname, "w");
    pdb_db = db_sel.create(fname, PM_DB_PDB);
    pdb_db->open(fname, PM_DB_MODE_WRITE, "");
    }
  }

//*============================================================*
//*==========              getRigid                  ==========*
//*============================================================*
// get a rigid simulation. 

PmRigidSimulation * 
PmSimulation::getRigid(const string name) 
  {
  string str;
  PmRigidSimulation *rsim = NULL; 
  int n = PmSimulation::simulations.size();

  for (int i = 0; i < n; i++) {
    PmSimulation *sim = PmSimulation::simulations[i];
    sim->getName (str);

    if (name == str) {
      rsim = dynamic_cast<PmRigidSimulation*>(sim);
      return rsim;
      }
    }

  return rsim;
  }

//*============================================================*
//*==========          setStateSaveFreq              ==========*
//*============================================================*
// get/set save frequency. 

void 
PmSimulation::setStateSaveFreq(const int freq) { 
  state_save_freq = freq; 
  }

void 
PmSimulation::getStateSaveFreq(int& freq) { 
  freq = state_save_freq; 
  }

//*============================================================*
//*==========              checkUpdateState          ==========*
//*============================================================*
// check for state update. 

bool
PmSimulation::checkUpdateState() 
  {
  if (num_steps == 0) {
    return true;
    }

  return !(num_steps % state_save_freq);
  }

//*============================================================*
//*==========              getKineticEnergy          ==========*
//*============================================================*

void
PmSimulation::getKineticEnergy(vector<float>& ke) {
  ke = kinetic_energy;
  }

void
PmSimulation::getKineticEnergy(const int n, float& val) 
  {
  if (this->kinetic_energy.size() <= n) {
    val = 0.0;
    return;
    }

  val = this->kinetic_energy[n];
  }

//*============================================================*
//*==========              getEnergy                 ==========*
//*============================================================*
// get potential energy.

void
PmSimulation::getEnergy(float& val) 
  {
  float ival;
  val = 0.0;

  for (int i = 0; i < interactions.size(); i++) {
    interactions[i]->getEnergy(ival);
    val += ival;
    }
  }

////////////////////////////////////////////////////////////////
//                 i n t e r a c t i o n s                   //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmInteraction::PmInteraction(const string name, vector<PmPotential*> list)
  {
  #ifdef dbg_PmInteraction_PmInteraction_name_list
  fprintf (stderr, "\n>>>>>> PmInteraction::ctor name = %s \n", name.c_str());
  fprintf (stderr, ">>> list.size = %d \n", list.size());
  #endif
  this->name = name;
  this->active = true;
  this->activated = true;
  this->verbose = false;
  this->general = true;

  for (int i = 0; i < list.size(); i++) {
    potentials.push_back(list[i]);
    }
  }

//*============================================================*
//*==========              getBodies                 ==========*
//*============================================================*

void
PmInteraction::getBodies(vector<PmBody*>& bodies) {
  bodies = this->bodies;
  }

//*============================================================*
//*==========              getName                   ==========*
//*============================================================*

void
PmInteraction::getName(string& name) {
  name = this->name;
  }

//*============================================================*
//*==========         getPotentialType               ==========*
//*============================================================*

void 
PmInteraction::getPotentialType(PmPotentialType& ptype) 
  {
  if (potentials.size()) { 
    ptype = potentials[0]->getType();
    }
  else {
    ptype = PM_POTENTIAL_UNKNOWN;
    }
  }

//*============================================================*
//*==========              getPotentials             ==========*
//*============================================================*

void 
PmInteraction::getPotentials(vector<PmPotential*>& potentials) 
  {

  potentials.clear();

  for (int i = 0; i < this->potentials.size(); i++) {
    potentials.push_back(this->potentials[i]);
    }
  }

//*============================================================*
//*==========           getNumInteractions           ==========*
//*============================================================*

void 
PmInteraction::getNumInteractions(int& num) 
  {
  num = 0;

  for (int i = 0; i < potentials.size(); i++) {
    num += potentials[i]->getNumInteractions();
    }
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*

void
PmInteraction::updateDisplay() 
  {
  /*
  fprintf (stderr, ">>> PmInteraction::updateDisplay: %s \n", name.c_str());
  fprintf (stderr, ">>> potentials.size  = %d  \n", potentials.size());
  */

  for (int i = 0; i < potentials.size(); i++) {
    potentials[i]->updateDisplay();
    }
  }

//*============================================================*
//*==========              setVerbose                ==========*
//*============================================================*

void 
PmInteraction::setVerbose(const bool flag) {
  verbose = flag;
  }

//*============================================================*
//*==========              addEnergy                 ==========*
//*============================================================*

void
PmInteraction::addEnergy(const float t)
  {
  float val;
  this->getEnergy(t, val);
  this->energy.push_back(val);
  /*
  fprintf (stderr, "\n>>>>>> PmInteraction::addEnergy time=%g val=%g size=%d \n", t,
           val, this->energy.size());
  */
  }

//*============================================================*
//*==========              getEnergy                 ==========*
//*============================================================*

void
PmInteraction::getEnergy(vector<float>& vals) {
  vals = energy;
  }

void
PmInteraction::getEnergy(float& tval) 
  {
  float pval;
  tval = 0.0;

  for (int i = 0; i < potentials.size(); i++) {
    potentials[i]->getEnergy(pval);
    tval += pval;
    }
  }

bool
PmInteraction::getEnergy(const int n, float& t, float& val) 
  {
  int ne = this->energy.size();
  //fprintf (stderr, ">>>>>> PmInteraction::getEnergy:  ne=%d \n", ne); 

  if ((ne == 0) || (ne <= n) || (n < 0)) {
    t = 0.0;
    val = 0.0;
    return false;
    }

  getTime(n, t);
  val = this->energy[n];
  return true;
  }

void
PmInteraction::getEnergy (const float time, float& val)
  {
  bool has_time = time_interval.hasTime(time);
  float tval = 0.0, pval;

  for (int i = 0; i < potentials.size(); i++) {
    potentials[i]->getEnergy(time, has_time, pval);
    tval += pval;
    }

  val = tval;
  }

//*============================================================*
//*==========              getTime                   ==========*
//*============================================================*

bool 
PmInteraction::getTime(const int n, float& t) 
  {
  int en = this->time.size();

  if ((en == 0) || (en <= n) || (n < 0)) {
    t = 0.0;
    return false;
    }

  t = this->time[n];
  return true;
  }

void 
PmInteraction::getTime(float& t)
  {
  t = this->time.back();
  }

//*============================================================*
//*==========              setActive                 ==========*
//*============================================================*

void
PmInteraction::setActive(const bool flag) 
  {
  activated = flag;
  active = flag;
  }

//*============================================================*
//*==========               isActive                 ==========*
//*============================================================*
// check if an interaction is active for the given time.

bool
PmInteraction::isActive(float t)
  {
  //fprintf (stderr, "\n>>>>>> PmInteraction::isActive: t=%g %d \n", t);

  if (!activated) {
    active = false;
    }

  else {
    bool has_time = time_interval.hasTime(t);

    // if the interaction is active then deactivate it //

    if (active && !has_time) {
      active = false;
      }

    else if (!active && has_time) {
      active = true;
      }
    }

  return active;
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*

void
PmInteraction::compForces(const bool upd, const float time) 
  {
  #ifdef dbg_PmInteraction_compForces
  fprintf (stderr, ">>>>>> PmInteraction::compForces \n");
  fprintf (stderr, ">>> potentials.size = %d \n", potentials.size());
  #endif

  if (!activated) {
    return;
    }

  bool has_time = time_interval.hasTime(time);

  for (int i = 0; i < potentials.size(); i++) {
    potentials[i]->compForces(upd, time, has_time, verbose);
    }
  }

//*============================================================*
//*==========              addStrain                 ==========*
//*============================================================*

void
PmInteraction::addStrain(const float t)
  {
  float val;
  this->getStrain(t, val);
  this->strain.push_back(val);
  }

//*============================================================*
//*==========              getStrain                 ==========*
//*============================================================*

void
PmInteraction::getStrain(const float time, float& val)
  {
  float tval = 0.0;
  bool has_time = time_interval.hasTime(time);

  for (int i = 0; i < potentials.size(); i++) {
    potentials[i]->getStrain(time, has_time, tval);
    val += tval;
    }
  }

//*============================================================*
//*==========              getStrain                 ==========*
//*============================================================*

void
PmInteraction::getStrain(vector<float>& vals) {
  vals = strain;
  }

bool
PmInteraction::getStrain(const int n, float& t, float& val) 
  {
  int ns = this->strain.size();

  if ((ns == 0) || (ns <= n) || (n < 0)) {
    t = 0.0;
    val = 0.0;
    return false;
    }

  getTime(n, t);
  val = this->strain[n];
  return true;
  }

//*============================================================*
//*==========              setTimeInterval           ==========*
//*============================================================*
// set the time interval for an interaction.

void 
PmInteraction::setTimeInterval(PmTimeInterval& tint) {
  time_interval = tint;
  }

}

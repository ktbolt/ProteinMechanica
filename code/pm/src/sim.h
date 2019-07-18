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
//* sim:                 s i m u l a t i o n                   *
//*============================================================*
// includes object definitions for PmSimulation and PmInteraction.

#ifndef _SIMULATION_PM_H_ 
#define _SIMULATION_PM_H_

#include "pm/pm.h"
#include "pot.h"
#include "bexp.h"

namespace ProteinMechanica {

class PmMolecule;
class PmSimulationObj;
class PmRigidSimulation;
class PmRestraint;

typedef enum {
  PM_SIMULATION_UNKNOWN,
  PM_SIMULATION_DEFORMED,
  PM_SIMULATION_RIGID,
  PM_SIMULATION_SIZE
  } PmSimulationType;


//*============================================================*
//*==========              PmInteraction             ==========*
//*============================================================*
// physical interaction between simulation objects. 

class PM_EXPORT PmInteraction {

    public:
       PmInteraction(const string name, vector<PmPotential*> list); 

       void getBodies(vector<PmBody*>& bodies);

       void addEnergy(const float t);
       void getEnergy(const float time, float& val);
       void getEnergy(float& val);
       void getEnergy(vector<float>& energy);
       bool getEnergy(const int n, float& t, float& energy);

       void compForces(const bool upd, const float time);

       void getName(string& name);
       void getNumInteractions(int& num);

       void getPotentials(vector<PmPotential*>& potentials);
       void getPotentialType(PmPotentialType& ptype);

       void addStrain(const float t);
       void getStrain(const float time, float& val);
       void getStrain(vector<float>& strain);
       bool getStrain(const int n, float& t, float& strain);

       void getTime(float& time);
       bool getTime(const int n, float& time);
       void setTimeInterval(PmTimeInterval& tint);
       bool isActive(float time);

       void updateDisplay();
       void setVerbose(const bool flag);
       void setActive(const bool flag);

    private:
       string name;
       vector<PmBody*> bodies, bodies1, bodies2; 
       vector<PmPotential*> potentials; 
       bool active, activated;
       bool verbose;
       bool general;
       PmTimeInterval time_interval;
       vector<float> time;
       vector<float> energy;
       vector<float> strain;
  };

//*============================================================*
//*==========              PmSimulation              ==========*
//*============================================================*

 class PM_EXPORT PmSimulation {
    public:
      PmSimulation();
      PmSimulation(const string name);
      void getName(string& name);
      PmSimulationType getType();
      static void convSimulationType(const string str, PmSimulationType& type);
      static PmSimulation* create(const string name, const PmSimulationType type);
      static vector<PmSimulation*> simulations; 
      static PmRigidSimulation* getRigid(const string name);

      void setReplay(const bool flag, const int start, const int end, const int inc);
      void setNumInc(int n);
      void getNumSteps(int& n);
      void setStateSaveFreq(const int freq);
      void getStateSaveFreq(int& freq);
      bool checkUpdateState();
      bool active();

      //bool PmSimulation::checkWhileExp();
      bool checkWhileExp();
      void setWhileExp (PmBooleanExpression *bexp);

      void getEnergy(float& val);
      void getKineticEnergy(vector<float>& energy);
      void getKineticEnergy(const int n, float& energy);

      void addRestraint(PmRestraint *res);
      void getRestraints(vector<PmRestraint*>& res);

      void addInteraction(PmInteraction *inter);
      void getInteractions(vector<PmInteraction*>& intr);

      void setSilentStep(const bool flag);
      void setWriteContact(const bool flag);
      void setWriteCurrent(const bool flag);
      void setWriteDomains(const bool flag);
      void setWriteEnergy(const bool flag);
      void setWriteJoints(const bool flag);
      void setWriteKineticEnergy(const bool flag);
      void setWriteResultsParams (string file_name, bool bformat);
      void setWriteState(const bool flag);
      void setWriteStrain(const bool flag);

      virtual void halt() = 0;
      virtual void initialize() = 0;  
      virtual void writeResults() = 0;
      virtual void readState(const string file_name) = 0;
      virtual void replay() = 0;
      virtual void step(int nsteps) = 0;

    protected:
      string name;
      PmSimulationType type;
      bool simulation_active;
      int num_steps;
      int num_inc, inc_count;
      int state_save_freq, num_state_upd;
      bool replay_active; 
      int replay_start, replay_end, replay_step, replay_inc;
      bool write_results; 
      bool silent_step_;

      PmBooleanExpression *while_exp;

      vector<PmRestraint*> restraints;
      vector<PmInteraction*> interactions;
      vector<float> energy;
      vector<float> kinetic_energy;
      vector<float> time;

      // results i/o //

      string write_file_name;
      bool write_state, write_energy, write_strain, write_domains, write_kenergy,
           write_contact_energy, write_joints;
      FILE *state_fp, *energy_fp, *strain_fp, *kenergy_fp, *joints_fp;
      bool binary_format;
      bool write_init, write_current;
      int last_written;
      PmDbInterface *pdb_db;

      void updateDisplay();
   };

}

#endif




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
//* bsim:      b i n d i n g   s i m u l a t i o n             *
//*============================================================*

#ifndef _BINDING_SIMULATION_PM_H_ 
#define _BINDING_SIMULATION_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

//*============================================================*
//*==========             PmBindingSimConf           ==========*
//*============================================================*

 class PM_EXPORT PmBindingSimConf {
    public:
      PmBindingSimConf();
      int getNextBindingSite();
      void setBindingSites(const int bind_site, const vector<int>& binding_sites);
      void setTargetSites(const int bind_site, const vector<int>&tsites);

      string name; 
      vector<string> head_names;
      string simulation_name, prev_simulation_name; 
      string energy_file_name, domain_file_name;
      int binding_site;
      int binding_site_num;
      vector<int> binding_sites;
      vector<int> target_sites;
      float energy_sd, energy_mean;
      int level;
      vector<class PmBindingSimConf*> next_head;
    };

//*============================================================*
//*==========         PmBindingSimulationResults     ==========*
//*============================================================*

typedef struct PmBindingSimulationLevel {
  float max_energy, min_energy;
  vector<PmBindingSimConf*> confs;
  } PmBindingSimulationLevel;

typedef struct PmBindingSimulationResults {
   vector<PmBindingSimulationLevel> levels;
   } PmBindingSimulationResults;

typedef struct PmBindingSimulationSummaryLevel {
  vector<string> names; 
  vector<float> energy_sd; 
  vector<float> energy_mean; 
  vector<PmMolecule*> conformations; 
  vector<PmVector3> cmap;
  } PmBindingSimulationSummaryLevel;

typedef struct PmBindingSimulationSummary {
  int num_levels, current_level;
  vector<int> level; 
  vector <PmBindingSimulationSummaryLevel> levels;
  } PmBindingSimulationSummary;

//*============================================================*
//*==========              PmBindingSimulation       ==========*
//*============================================================*

 class PM_EXPORT PmBindingSimulation {
    public:
      PmBindingSimulation();
      PmBindingSimulation(const string name);
      void getName(string& name);
      void setBindingHeads(vector<string>& heads);
      void setBindingSites(vector<string>& sites);
      void setDomainCmdFiles(vector<string>& list);
      void setFilamentName(const string val);
      void setModelCmdFiles(vector<string>& list);
      void setMotorName(const string val);
      void setPmProgramName(const string val);
      void setPotentialCmdFiles(vector<string>& list);
      void setRunCmdFile(string val);
      void setSimulationCmdFiles(vector<string>& list);
      void setNumThreads(const int num_threads);
      void writeSummary();
      void setTargetRestraintsCmdFiles(vector<string>& list);
      void setTargetName(const string val);
      void setStructureFileName(const string val);
      void execute();
      void processSimulationResults (vector<PmBindingSimConf*> new_confs, 
                                     vector<PmBindingSimConf*>& current_confs);
      void displayConfs(const bool show);
      void setSummaryFile(const string summary_file);
      void setConformationsFile(const string conf_file);
      void setIncrConformation(const bool incr_confs);
      void incConformation(const int incr);
      void setLevel(const int level);

    private:
      string name;
      string motor_name;
      string filament_name;
      string target_name;
      string structure_file_name;
      string conformations_file_name;
      string run_cmd_file_;
      string pm_binary_;
      int num_threads_;
      vector<string> binding_heads;
      vector<int> binding_sites;
      vector<string> domain_cmd_files_;
      vector<string> model_cmd_files_;
      vector<string> potential_cmd_files_;
      vector<string> simulation_cmd_files_;
      vector<string> target_restraints_cmd_files_;
      vector<string> energy_files_;
      vector<float>  energy_results_;
      int energy_window;
      PmBindingSimulationResults results;
      PmBindingSimulationSummary summary;
      bool increment_confs;
      int conformation_number, prev_conformation_number;

      void init();
      bool processStructure();
      bool checkData();
      string buildCmd(PmBindingSimConf *conf, const vector<string>& heads);
      void processEnergyFile(const string energy_file);
      bool execCmdScripts(const vector<string>& cmds, const int num_threads);
      void writeConfs();
   };
   
}

#endif




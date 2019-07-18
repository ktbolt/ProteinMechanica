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
//* bsim:    b i n d i n g       s i m u l a t i o n           *
//*============================================================*

#include "bsim.h"
#include "db.h"
#include <list>
#include "pthread.h"
#include <fstream>

namespace ProteinMechanica {

//===== typedef and function for sorting energy =====//

typedef std::pair<float, PmBindingSimConf*> EnergyConfPair;

static bool 
compareEnergy (const EnergyConfPair& l, const EnergyConfPair& r) { 
    return l.first < r.first; 
    }

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmBindingSimulation::PmBindingSimulation () 
  {
  init();
  }

PmBindingSimulation::PmBindingSimulation (const string name) {
  this->name = name;
  init();
  }

void
PmBindingSimulation::init()
  {
  //fprintf (stderr, "\n>>>>>> PmBindingSimulation::init \n");
  pm_binary_ = "pm";
  num_threads_ = 4;
  energy_window = 100;
  increment_confs = false;
  conformation_number = 0;
  prev_conformation_number = 0;

  summary.num_levels = 0;
  summary.current_level = 0;
  }

//*============================================================*
//*==========              getName                   ==========*
//*============================================================*

void 
PmBindingSimulation::getName(string& name) {
  name = this->name;
  }

//*============================================================*
//*==========              set*Name                  ==========*
//*============================================================*

void 
PmBindingSimulation::setTargetName(const string val) {
  this->target_name = val;
  }

void 
PmBindingSimulation::setStructureFileName(const string val) {
  this->structure_file_name = val;
  }

void 
PmBindingSimulation::setFilamentName(const string val) {
  filament_name = val;
  }

void 
PmBindingSimulation::setMotorName(const string val) {
  motor_name = val;
  }

//*============================================================*
//*==========              setBindingHeads           ==========*
//*============================================================*

void
PmBindingSimulation::setBindingHeads(vector<string>& heads) {
  binding_heads = heads;
  }

//*============================================================*
//*==========              setBindingSites           ==========*
//*============================================================*

void 
PmBindingSimulation::setBindingSites(vector<string>& sites)
  {
  for (unsigned int i = 0; i < sites.size(); i++) {
    this->binding_sites.push_back(atoi(sites[i].c_str()));
    }
  }

//*============================================================*
//*==========              setPmProgramName          ==========*
//*============================================================*

void
PmBindingSimulation::setPmProgramName(const string val)
  {
  this->pm_binary_ = val;
  }

//*============================================================*
//*==========              set*CmdFiles              ==========*
//*============================================================*
// set the names of the PM command files needed for a simulation.

void 
PmBindingSimulation::setDomainCmdFiles(vector<string>& list) {
  domain_cmd_files_ = list;
  }

void 
PmBindingSimulation::setModelCmdFiles(vector<string>& list) {
  model_cmd_files_ = list;
  }

void 
PmBindingSimulation::setPotentialCmdFiles(vector<string>& list) {
  potential_cmd_files_ = list;
  }

void 
PmBindingSimulation::setRunCmdFile(string val) {
  run_cmd_file_ = val;
  }

void 
PmBindingSimulation::setSimulationCmdFiles(vector<string>& list) {
  simulation_cmd_files_ = list;
  }

void 
PmBindingSimulation::setTargetRestraintsCmdFiles(vector<string>& list) {
  target_restraints_cmd_files_ = list;
  }

void 
PmBindingSimulation::setNumThreads(const int num_threads) {
  num_threads_ = num_threads;
  }

//*============================================================*
//*==========              checkData                 ==========*
//*============================================================*

bool 
PmBindingSimulation::checkData()
  {
  if (target_name.empty()) {
    pm_ErrorReport (PM, "no target name given.", "*");
    return false;
    }

  if (structure_file_name.empty()) {
    pm_ErrorReport (PM, "no structure file name given.", "*");
    return false;
    }

  if (binding_sites.empty()) {
    pm_ErrorReport (PM, "no binding sites given.", "*");
    return false;
    }

  if (target_name.empty()) {
    pm_ErrorReport (PM, "no target name given.", "*");
    return false;
    }

  if (domain_cmd_files_.empty()) {
    pm_ErrorReport (PM, "no domain cmd files given.", "*");
    return false;
    }

  if (model_cmd_files_.empty()) {
    pm_ErrorReport (PM, "no model cmd files given.", "*");
    return false;
    }

  if (potential_cmd_files_.empty()) {
    pm_ErrorReport (PM, "no potential cmd files given.", "*");
    return false;
    }

  if (simulation_cmd_files_.empty()) {
    pm_ErrorReport (PM, "no simulation cmd files given.", "*");
    return false;
    }

  if (target_restraints_cmd_files_.empty()) {
    pm_ErrorReport (PM, "no target restraints cmd files given.", "*");
    return false;
    }

  return true;
  }

//*============================================================*
//*==========              execute                   ==========*
//*============================================================*

void
PmBindingSimulation::execute()
  {
  fprintf (stderr, "\n----- execute binding simulation ----- \n");

  //===== check data =====//

  if (!checkData()) {
    return;
    }

  //===== process structure file =====//

  if (!processStructure()) {
    return;
    }

  //===== cycle through heads and binding sites =====//

  string cmd, head, sim_name, bind_name, conf_name;
  int bind_site; 
  vector<string> cmd_list;
  stringstream dss;
  bool finished; 
  int num_sims;
  int max_sims = 10;
  int num_heads = binding_heads.size();
  //num_heads = 1;
  PmBindingSimConf *conf, *new_conf;
  vector<PmBindingSimConf*> current_confs;
  vector<PmBindingSimConf*> new_confs;
  PmBindingSimulationLevel sim_level;

  // set initial conformation //

  conf = new PmBindingSimConf;
  conf->head_names.push_back(binding_heads[0]);
  conf->binding_sites = binding_sites;
  current_confs.push_back(conf);

  for (int i = 0; i < num_heads; i++) {
    cmd_list.clear();
    head = binding_heads[i];
    fprintf (stderr, "\n----- head=%s ----- \n", head.c_str());
    fprintf (stderr, ">>> num confs=%d\n", current_confs.size()); 
    num_sims = 0;
    finished = false; 
    new_confs.clear();

    while (!finished) {
      finished = true; 

      for (int j = 0; j < current_confs.size(); j++) {
        conf = current_confs[j];
        bind_site = conf->getNextBindingSite();
        bind_name = conf->name;

        if (bind_site != -1) { 
          num_sims += 1;
          finished = false; 

          if (bind_name.empty()) {
            dss << head << "Bind" << bind_site;
            }
          else {
            dss << bind_name << "_" << head << "Bind" << bind_site;
            }

          conf_name = dss.str();
          dss.str(std::string());

          dss << name << "_" << conf_name;
          sim_name = dss.str();
          dss.str(std::string());
          fprintf (stderr, ">>> sim=%s\n", sim_name.c_str());
          if (i+1 > 1) {
          fprintf (stderr, "    prev sim=%s\n", conf->simulation_name.c_str());
          }

          //fprintf (stderr, ">>> cmd=%s\n", cmd.c_str());

          new_conf = new PmBindingSimConf;
          new_conf->name = conf_name; 
          new_conf->binding_site = bind_site; 
          new_conf->level = i+1; 
          new_conf->simulation_name = sim_name; 
          new_conf->prev_simulation_name = conf->simulation_name; 
          new_conf->head_names.push_back(binding_heads[0]);
          new_conf->setBindingSites(bind_site, binding_sites);
          new_conf->setTargetSites(bind_site, conf->target_sites);
          new_confs.push_back(new_conf);
          /*
          fprintf (stderr, "    tsites=");
          for (int ts=0; ts < new_conf->target_sites.size(); ts++) {
            fprintf (stderr, "%d ", new_conf->target_sites[ts]);
            }
          fprintf (stderr, "\n");
          */

          cmd = buildCmd(new_conf, binding_heads);
          cmd_list.push_back(cmd);
          //fprintf (stderr, ">>> cmd=%s \n", cmd.c_str());
          }
        }
      }

    // execute the pm commands //

    execCmdScripts (cmd_list, num_threads_);

    // process simulation results // 

    processSimulationResults (new_confs, current_confs);

    sim_level.confs = current_confs;
    results.levels.push_back(sim_level);
    }

  fprintf (stderr, "\n\n>>> done! \n");

  //===== write conformations =====//

  writeConfs();

  //===== write summary =====//

  writeSummary();
  }

//*============================================================*
//*==========              buildCmd                  ==========*
//*============================================================*

string 
PmBindingSimulation::buildCmd(PmBindingSimConf *conf, const vector<string>& heads)
  {
  stringstream dss, dss1;
  string cmd, sim_file_name, dom_file_name;

  //===== set executable =====//

  dss << pm_binary_ << " silent=true "; 

  //===== set variables =====//

  dss << " structures_file=" << structure_file_name;
  dss << " motor_name=" << motor_name; 
  dss << " actin_name=" << filament_name; 

  for (unsigned int i = 0; i < conf->target_sites.size(); i++) { 
    dss << " target" << i+1 << "_name=" << target_name << conf->target_sites[i]; 
    dss << " binding_head" << i+1 << "=" << heads[i];
    }

  dss << " simulation_name=" << conf->simulation_name; 

  if (conf->level > 1) {
    dss << " simulation_state_file_name=" << conf->prev_simulation_name; 
    //fprintf (stderr, "    sim state file=%s \n", conf->prev_simulation_name.c_str());
    }

  dss << " simulation_file_name=" << conf->simulation_name;
  dss1 << conf->simulation_name << "_energy.pm"; 
  sim_file_name = dss1.str();
  dss1.str(std::string());
  energy_files_.push_back(sim_file_name);
  conf->energy_file_name = sim_file_name;
  dss1 << conf->simulation_name << "_domains.pdb"; 
  conf->domain_file_name = dss1.str();
  dss1.str(std::string());

  //===== set cmd files =====//

  for (unsigned int i = 0; i < domain_cmd_files_.size(); i++) { 
    dss << " " << domain_cmd_files_[i]; 
    }

  for (unsigned int i = 0; i < model_cmd_files_.size(); i++) { 
    dss << " " << model_cmd_files_[i]; 
    }

  for (unsigned int i = 0; i < simulation_cmd_files_.size(); i++) { 
    dss << " " << simulation_cmd_files_[i]; 
    }

  for (unsigned int i = 0; i < potential_cmd_files_.size(); i++) { 
    dss << " " << potential_cmd_files_[i]; 
    }

  for (unsigned int i = 0; i < target_restraints_cmd_files_.size(); i++) { 
    dss << " " << target_restraints_cmd_files_[i]; 
    }

  dss << " " << run_cmd_file_; 

  //dss << " >& log.out "; 
  cmd = dss.str();
  dss.str(std::string());
  return cmd;
  }

//*============================================================*
//*==========              writeSummary              ==========*
//*============================================================*

void
PmBindingSimulation::writeSummary()
  {
  fprintf (stderr, "\n>>>>>> PmBindingSimulation::writeSummary \n");
  string fname;
  stringstream dss;
  FILE *fp;

  vector<PmBindingSimConf*> confs;
  PmBindingSimConf *conf;
  vector<EnergyConfPair> sorted_confs;

  dss << name << "_summary.txt"; 
  fname = dss.str();
  dss.str(std::string());
  fp = fopen (fname.c_str(), "w");

  fprintf (fp, "# Protein Mechanica binding simulation summary \n");
  fprintf (fp, "simulation \"%s\" \n", name.c_str());
  fprintf (fp, "number of levels=%d \n", results.levels.size()); 

  //===== write all levels =====//

  for (unsigned int j = 0; j < results.levels.size(); j++) {
    confs = results.levels[j].confs; 
    sorted_confs.clear();

    // sort confs by energy //

    for (unsigned int i = 0; i < confs.size(); i++) {
      EnergyConfPair pair(confs[i]->energy_mean, confs[i]);
      sorted_confs.push_back(pair);
      }

    std::sort (sorted_confs.begin(), sorted_confs.end(), compareEnergy);

    for (unsigned int i = 0; i < confs.size(); i++) {
      conf = sorted_confs[i].second;
      fprintf (fp, "level=%d conformation=%s mean=%f sd=%f \n", j+1, conf->name.c_str(),
               conf->energy_mean, conf->energy_sd);
      }
    }

  fclose (fp);
  }

//*============================================================*
//*==========         processSimulationResults       ==========*
//*============================================================*

void 
PmBindingSimulation::processSimulationResults (vector<PmBindingSimConf*>new_confs,
                                               vector<PmBindingSimConf*>& current_confs)
  {
  #ifdef PmBindingSimulation_processSimulationResults 
  fprintf (stderr, "\n>>>>>> PmBindingSimulation::processSimulationResults \n");
  #endif
  fprintf (stderr, "\n>>>>>> PmBindingSimulation::processSimulationResults \n");
  PmBindingSimConf *conf;
  string fname;
  vector<string> names;
  vector<vector<float> > values;
  vector<float> total_energy;
  int num_rows, num_cols, start_num;
  float mean, sd;

  for (unsigned int c = 0; c < new_confs.size(); c++) {
    conf = new_confs[c];
    fname = conf->energy_file_name;
    PmDbPmInterface::readEnergy(fname, names, values);
    num_rows = values.size();
    num_cols = values[0].size();
    #ifdef PmBindingSimulation_processSimulationResults 
    fprintf (stderr, ">>> name=%s num_rows=%d num_cols=%d \n", conf->name.c_str(), 
             num_rows, num_cols);
    #endif

    total_energy.clear();

    if (num_rows < energy_window) {
      start_num = 0;
      }
    else {
      start_num = num_rows - energy_window;
      }

    for (int i = start_num; i < num_rows; i++) {
      total_energy.push_back(values[i][num_cols-1]);
      }

    pm_MathStandardDeviationComp (total_energy, mean, sd);
    conf->energy_mean = mean;
    conf->energy_sd = sd;
    #ifdef PmBindingSimulation_processSimulationResults 
    fprintf (stderr, "    energy mean=%f sd=%f \n", mean, sd); 
    #endif
    }

  current_confs = new_confs;
  }

//*============================================================*
//*==========         processEnergyFile              ==========*
//*============================================================*
void
PmBindingSimulation::processEnergyFile(const string fname)
  {
  }

//*============================================================*
//*==========              processStructure          ==========*
//*============================================================*

bool 
PmBindingSimulation::processStructure()
  {
  PmDbInterfaceSelect db_sel;
  PmDbType format;
  stringstream dss;
  string data_type, tname;
  vector<string> model_names;
  bool found;
  PmMolecule *mol;
  PmMoleculeType mol_type;
  PmExtent extent;
  PmDbInterface *db;

  db = db_sel.create(name, PM_DB_PDB);
  db->open(structure_file_name, PM_DB_MODE_READ, data_type);
  db->getModelNames(model_names);

  fprintf (stderr, ">>> model names=");
  for (unsigned int i = 0; i < model_names.size(); i++) {
    fprintf (stderr, "%s ", model_names[i].c_str());
    }
  fprintf (stderr, "\n");

  //===== check binding sites =====//

  for (unsigned int i = 0; i < binding_sites.size(); i++) {
    dss << target_name << binding_sites[i];
    tname = dss.str();
    dss.str(std::string());
    found = false;

    for (unsigned int j = 0; j < model_names.size(); j++) {
      if (tname == model_names[j]) {
        found = true;
        break;
        }
      }  

    if (!found) { 
      pm_ErrorReport (PM, "no binding site for %s found.", "*", tname.c_str());
      return false;
      }
    }  

  return true;
  }

//*============================================================*
//*==========              writeConfs                ==========*
//*============================================================*

void 
PmBindingSimulation::writeConfs()
  {
  fprintf (stderr, "\n>>>>>> PmBindingSimulation::writeConfs \n");
  vector<PmBindingSimConf*> confs;
  PmBindingSimConf *conf;
  string name, fname, data_type, cfname;
  PmDbInterfaceSelect rdb_sel, wdb_sel;
  PmDbInterface *rdb, *wdb;
  PmMolecule *mol;
  stringstream dss;
  vector<PmMolecule*> mols(1);
  string chain_id, sidechains, res_names, model_name;
  bool mainchain, no_hydrogen;
  vector<string> chains;
  PmAtomFilter filter;
  string one_chain, start_res;
  int model;

  //===== open database for writing structure =====//

  dss << this->name << "_conformations.pdb"; 
  cfname = dss.str();
  dss.str(std::string());
  wdb = wdb_sel.create(fname, PM_DB_PDB);

  for (unsigned int j = 0; j < results.levels.size(); j++) {
    confs = results.levels[j].confs; 

    for (unsigned int i = 0; i < confs.size(); i++) {
      conf = confs[i];
      fprintf (stderr, ">>> conf name=%s \n", conf->name.c_str()); 
      fprintf (stderr, "    dom fname=%s \n", conf->domain_file_name.c_str()); 
      name = conf->name; 
      fname = conf->domain_file_name; 
      rdb = rdb_sel.create(fname, PM_DB_PDB);
      rdb->open(fname, PM_DB_MODE_READ, data_type);
      mol = rdb->getMolecule(1, PM_MOLECULE_PROTEIN);
      rdb->close();

      model = i+1;
      model_name = conf->name;
      mols[0] = mol;

      if ((j ==0) && (i == 0)) {
        fprintf (stderr, ">>> cfname=%s \n", cfname.c_str()); 
        wdb->open(cfname, PM_DB_MODE_WRITE, data_type);
        wdb->writeDomainAtoms(mols, chains, one_chain, no_hydrogen, start_res, res_names,
                              model, model_name);
        wdb->close();
        wdb->open(cfname, PM_DB_MODE_APPEND, data_type);
        }
      else {
        wdb->writeDomainAtoms(mols, chains, one_chain, no_hydrogen, start_res, res_names,
                              model, model_name);
        }
      }
    }

  wdb->close();
  }

//*============================================================*
//*==========         setConformationsFile           ==========*
//*============================================================*

void
PmBindingSimulation::setConformationsFile(const string fname)
  {
  conformations_file_name = fname;
  }

//*============================================================*
//*==========                setLevel                ==========*
//*============================================================*

void
PmBindingSimulation::setLevel(const int level)
  {
  if ((level > 0) && (level <= summary.levels.size())) {
    summary.current_level = level-1;
    }
  }

//*============================================================*
//*==========         setSummaryFile                 ==========*
//*============================================================*

void
PmBindingSimulation::setSummaryFile(const string fname)
  {
  ifstream sum_file;
  string line, name,value;
  int bpos, pos, last_pos, level;
  char c;
  bool proc_name;
  float mean;
  PmBindingSimulationSummaryLevel sim_level;

  sum_file.open (fname.c_str());
  last_pos = 0;
  cout << "---------- summary file ----------" << "\n";

  while (!sum_file.eof()) {
    getline (sum_file, line);

    if ((line == "") || (line[0] == '#')) {
      continue;
      }

    name.clear();
    proc_name = true;
    //cout << "\n >>> line: " << line << "\n";

    for (int n = 0; n < line.size();) { 
      c = line[n++]; 

      if (c == '\n') {
        break;
        }

      if (c == '=') {
        proc_name = false;
        value.clear();
        }
      else if ((c == ' ') && !name.empty()) {
        //cout << "name=" << name << " value=" << value << "\n";

        if (name == "conformation") {
          summary.levels[level].names.push_back(value);
          cout << "level=" << level+1 << " conformation=" << value;
          }
        else if (name == "level") {
          level = atoi(value.c_str());
          summary.levels.resize(level);
          level -= 1;
          }
        else if (name == "mean") {
          summary.levels[level].energy_mean.push_back(atof(value.c_str()));
          cout << " mean=" << value;
          }
        else if (name == "sd") {
          summary.levels[level].energy_sd.push_back(atof(value.c_str()));
          cout << " sd=" << value << "\n";
          }

        proc_name = true;
        name.clear();
        }

      else if (proc_name) {
        name.push_back(c);
        }
      else if (!name.empty()) {
        value.push_back(c);
        }
      }
    }

  /*
  for (unsigned int i = 0; i < summary.names.size(); i++) {
    name = summary.names[i];
    mean = summary.energy_mean[i];
    cout << "name=" << name << " mean =" << mean << "\n";
    }
  */
  }

//*============================================================*
//*==========              setIncrConformation       ==========*
//*============================================================*

void
PmBindingSimulation::setIncrConformation(const bool incr_confs)
  {
  increment_confs = incr_confs;
  }

//*============================================================*
//*==========              incConformation           ==========*
//*============================================================*

void
PmBindingSimulation::incConformation(const int incr)
  {
  int current_level = summary.current_level;
  int n = summary.levels[current_level].names.size();

  if (n == 0) { 
    return;
    }

  PmMolecule *mol;
  string name;

  prev_conformation_number = conformation_number;
  conformation_number += incr;

  if (conformation_number == n) {
    conformation_number = 0; 
    }
  else if (conformation_number == -1) {
    conformation_number = n-1; 
    }

  name = summary.levels[current_level].names[prev_conformation_number];
  mol = summary.levels[current_level].conformations[prev_conformation_number];
  mol->displayBackbone (name, false);

  name = summary.levels[current_level].names[conformation_number];
  mol = summary.levels[current_level].conformations[conformation_number];
  mol->setColor(summary.levels[current_level].cmap[conformation_number]);
  mol->setLineWidth(2.0);
  mol->displayBackbone (name, true);

  fprintf (stderr, "\n>>> number=%d level=%d conformation=%s  energy=%f \n", 
           conformation_number+1, current_level+1, name.c_str(), 
           summary.levels[current_level].energy_mean[conformation_number]);

  // make sure the graphics image buffer is updated. //

  pmSystem.updateGraphics(true);
  pmSystem.updateGraphics();
  }

//*============================================================*
//*==========              displayConfs              ==========*
//*============================================================*

void 
PmBindingSimulation::displayConfs(const bool show)
  {
  //fprintf (stderr, "\n>>>>>> PmBindingSimulation::displayConfs \n");
  vector<PmBindingSimConf*> confs;
  PmBindingSimConf *conf;
  string name, fname, data_type;
  PmDbInterfaceSelect db_sel;
  PmDbInterface *db;
  PmMolecule *mol;
  PmBindingSimulationLevel level;
  int current_level = summary.current_level;

  if (results.levels.size()) {
    level = results.levels.back();
    confs = level.confs;

    for (unsigned int i = 0; i < confs.size(); i++) {
      conf = confs[i];
      fprintf (stderr, ">>> conf name=%s \n", conf->name.c_str()); 
      fprintf (stderr, "    dom fname=%s \n", conf->domain_file_name.c_str()); 
      name = conf->name; 
      fname = conf->domain_file_name; 
      db = db_sel.create(fname, PM_DB_PDB);
      db->open(fname, PM_DB_MODE_READ, data_type);
      mol = db->getMolecule(1, PM_MOLECULE_PROTEIN);
      mol->displayBackbone (name, show);
      }
    }

  //===== display conformations from file =====//

  else if (summary.levels[current_level].names.size() && 
           !conformations_file_name.empty()) {
    fprintf (stderr, ">>> read conformations for level=%d \n", current_level+1); 
    vector<PmVector3> cmap;
    PmGraphicsInterface::getSpectrumColorMap(summary.levels[current_level].names.size(), 
                                             cmap);
    db = db_sel.create(conformations_file_name, PM_DB_PDB);
    db->open(conformations_file_name, PM_DB_MODE_READ, data_type);

    for (unsigned int i = 0; i < summary.levels[current_level].names.size(); i++) {
      name = summary.levels[current_level].names[i];
      fprintf (stderr, ">>> conf name=%s \n", name.c_str()); 
      mol = db->getMolecule(name, PM_MOLECULE_PROTEIN);

      if (increment_confs) {
        summary.levels[current_level].conformations.push_back(mol);
        summary.levels[current_level].cmap = cmap; 

        if (i == 0) {
          mol->setColor(cmap[i]);
          mol->setLineWidth(2.0);
          mol->displayBackbone (name, show);
          }
        }
      else {
        mol->setColor(cmap[i]);
        mol->setLineWidth(2.0);
        mol->displayBackbone (name, show);
        }
      }
    }
  }

////////////////////////////////////////////////////////////////
//              t h r e a d   f u n c t i o n s              //
//////////////////////////////////////////////////////////////

#define use_threads
#ifdef use_threads
pthread_mutex_t task_queue_lock;
pthread_mutex_t finished_queue_lock;

list<string> task_queue;
list<string> finished_queue;

typedef struct ThreadData {
  int id;
  } ThreadData;

extern "C" {
  void *workerThread(void *arg);
  }
#endif

//*============================================================*
//*==========              execCmdScripts            ==========*
//*============================================================*

bool
PmBindingSimulation::execCmdScripts(const vector<string>& cmds, const int num_threads)
  {
  fprintf (stderr, "\n>>>>>> PmBindingSimulation::execCmdScripts \n"); 
  unsigned long comp_DONE = 0;
  unsigned long comp_START = 0;

  //===== set-up the mutexes =====//

  pthread_mutex_init (&task_queue_lock,  NULL);
  pthread_mutex_init (&finished_queue_lock, NULL);

  int num_tasks = cmds.size();
  fprintf (stderr, ">>> num threads=%d \n", num_threads); 
  fprintf (stderr, ">>> num tasks=%d \n", num_tasks); 
  /*
  fprintf (stderr, "------ commands ------ \n"); 
  for (unsigned long i = 0; i < num_tasks; i++) {
    fprintf (stderr, ">>> %d command=\"%s\" \n", i, cmds[i].c_str()); 
    }
  */

  //===== fill task_queue with pm cmds =====//

  task_queue.clear();
  finished_queue.clear();

  for (unsigned i = 0; i < num_tasks; i++) {
    task_queue.push_back(cmds[i]);
    }

  //===== start worker threads =====//

  list<pthread_t*>    threadIdList;
  list<ThreadData> thread_table;

  for (unsigned long i = 0; i < num_threads; i++) {
    pthread_t *tId = new pthread_t;
    threadIdList.push_back(tId);
    ThreadData Y;
    Y.id = i;
    thread_table.push_back(Y);

    int rc = pthread_create (tId, NULL, workerThread, (void*)(&(thread_table.back())));

    if (rc) {
      cout << "ERROR; return code from pthread_create() " << comp_START << "\n";
      cout.flush();
      return false; 
      }
    }

  //===== wait for threads to finish =====//

  //fprintf (stderr, ">>> PmBindingSimulation::execCmdScripts: wait for threads. \n"); 
  int lc=0;

  while (comp_DONE != num_tasks) {
    pthread_mutex_lock (&task_queue_lock);
    comp_START = num_tasks - task_queue.size();
    pthread_mutex_unlock (&task_queue_lock);

    pthread_mutex_lock (&finished_queue_lock);
    comp_DONE = finished_queue.size();
    pthread_mutex_unlock (&finished_queue_lock);
    }

  //===== call join to kill all worker threads =====//

  //fprintf (stderr, ">>> PmBindingSimulation::execCmdScripts: call join. \n"); 

  list<pthread_t*>::iterator i = threadIdList.begin();

  while (i != threadIdList.end()) {
    if (pthread_join(*(*i), NULL)!=0) {
      cout << "Thread join error!\n";
      return false; 
      }

    delete (*i);
    threadIdList.erase(i++);
    }

  pthread_mutex_destroy (&task_queue_lock);
  pthread_mutex_destroy (&finished_queue_lock);
  //fprintf (stderr, ">>> PmBindingSimulation::execCmdScripts: finished. \n"); 
  return true;
  }

//*============================================================*
//*==========              workerThread              ==========*
//*============================================================*

void*
workerThread (void *threadarg) 
  {
#ifdef use_threads
  ThreadData *my_data;
  my_data = (ThreadData*)threadarg;
  int taskid = my_data->id;

  stringstream ss; 
  ss << taskid; 
  string taskString = ss.str();
  bool more_tasks = true;

  // keep on working until task_queue is empty  //

  while (more_tasks) {
    pthread_mutex_lock (&task_queue_lock);
    string workOnMe;

    if (task_queue.size() == 0) { 
      more_tasks = false; 
      }
    else {
      workOnMe = task_queue.front();
      task_queue.pop_front();
      }

    pthread_mutex_unlock (&task_queue_lock);

    if (!more_tasks) {
      break;
      }

    // execute pm cmd // 

    system (workOnMe.c_str());

    workOnMe = "thread " + taskString + " worked on " + workOnMe;

    // let's pretend this takes some time, add a delay to the computation

    /*
    struct timeval timeout;
    timeout.tv_sec = 0;
    timeout.tv_usec = 100000; // 0.1 second delay
    select (0, NULL, NULL, NULL, & timeout);
    */

    pthread_mutex_lock (&finished_queue_lock);
    finished_queue.push_back (workOnMe);
    pthread_mutex_unlock (&finished_queue_lock);
    }

  pthread_exit (NULL);
#endif
  }

////////////////////////////////////////////////////////////////
//              c o n f o r m a t i o n s                    //
//////////////////////////////////////////////////////////////

PmBindingSimConf::PmBindingSimConf() 
  {
  level = 0;
  energy_mean = 0.0;
  energy_sd = 0.0;
  binding_site_num = 0;
  }

//*============================================================*
//*==========              getNextBindingSite        ==========*
//*============================================================*

int
PmBindingSimConf::getNextBindingSite()
  {
  int bsite;

  if (binding_site_num == binding_sites.size()) {
    return -1;
    }

  bsite = binding_sites[binding_site_num];
  binding_site_num += 1;
  return bsite;
  }

//*============================================================*
//*==========              setBindingSites           ==========*
//*============================================================*

void 
PmBindingSimConf::setBindingSites(const int bind_site, 
                                                 const vector<int>& bsites)
  {
  int n;

  for (unsigned int i = 0; i <  bsites.size(); i++) {
    n = bsites[i];

    if (n > bind_site) {
      binding_sites.push_back(n);
      }
    }
  }

//*============================================================*
//*==========              setTargetSites            ==========*
//*============================================================*

void 
PmBindingSimConf::setTargetSites(const int bind_site, 
                                                const vector<int>&tsites)
  {
  for (unsigned int i = 0; i < tsites.size(); i++) {
    target_sites.push_back(tsites[i]);
    }

  target_sites.push_back(bind_site);
  }


};



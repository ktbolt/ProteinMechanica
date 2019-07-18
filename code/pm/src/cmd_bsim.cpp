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

////////////////////////////////////////////////////////////////
//  b i n d i n g   s i m u l a t i o n   c o m m a n d s    //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pm_CmdBindSimulation      ==========*
//*============================================================*
// process binding simulation command.

void
pm_CmdBindSimulation (PmCmdDataList& dlist)
  {
  string dv, name, dname, res, type_str, sim_name, str;
  PmBindingSimulation *sim;
  PmCmdData data;

  dlist.getNext (data, dv);

  //===== create a new simulation and add it to the system =====//

  if (data.name == "create") {
    dlist.getString("name", name);
    sim = new PmBindingSimulation(name);
    pmSystem.addBindSimulation(sim);
    sim_name = name;
    }
  else {
    pmSystem.getBindSimulation (data.name, &sim);
    sim_name = data.name;
    }

  if (!sim) {
    pm_ErrorReport (PM, "no simulation named \"%s\".", "*", data.name.c_str());
    return;
    }

  string target_name, structure_fname, run_cmd_file, motor_name, fil_name;
  string program_name;
  vector<string> bind_sites, domains, models, pots, targets, sims, bind_heads;
  bool execute, print_summary, show;

  if (dlist.getString("target_name", target_name)) {
    sim->setTargetName(target_name);
    }

  if (dlist.getString("structure_file", structure_fname)) {
    sim->setStructureFileName(structure_fname);
    }

  if (dlist.getString("motor_name", motor_name)) {
    sim->setMotorName(motor_name);
    }

  if (dlist.getString("filament_name", fil_name)) {
    sim->setFilamentName(fil_name);
    }

  if (dlist.getString("pm_program_name", program_name)) {
    sim->setPmProgramName(program_name);
    }

  if (dlist.getStringList("binding_heads", bind_heads)) {
    sim->setBindingHeads(bind_heads);
    }

  if (dlist.getStringList("binding_sites", bind_sites)) {
    sim->setBindingSites(bind_sites);
    }

  if (dlist.getStringList("domain_cmd_files", domains)) {
    sim->setDomainCmdFiles(domains);
    }

  if (dlist.getStringList("model_cmd_files", models)) {
    sim->setModelCmdFiles(models);
    }

  if (dlist.getStringList("potential_cmd_files", pots)) {
    sim->setPotentialCmdFiles(pots);
    }

  if (dlist.getStringList("simulation_cmd_files", sims)) {
    sim->setSimulationCmdFiles(sims);
    }

  if (dlist.getStringList("target_restraints_cmd_files", targets)) {
    sim->setTargetRestraintsCmdFiles(targets);
    }

  if (dlist.getString("run_cmd_file", run_cmd_file)) {
    sim->setRunCmdFile(run_cmd_file);
    }

  if (dlist.getBoolean("execute", execute)) {
    if (execute) {
      sim->execute();
      }

    return;
    }

  //===== visualize conformations =====//

  if (dlist.getString("conformations", str)) { 
    string summary_file, conformations_file;
    bool incr_confs;
    int inc, level;

    if (dlist.getString("summary_file", summary_file)) {
      sim->setSummaryFile(summary_file);
      }

    if (dlist.getString("conformations_file", conformations_file)) {
      sim->setConformationsFile(conformations_file);
      }

    // increment conformations using arrow keys //

    if (dlist.getBoolean("increment_conformations", incr_confs)) { 
      string cmd;
      stringstream dss;
      sim->setIncrConformation(incr_confs);

      dss << "binding_simulation " << sim_name << " conformations increment=1 "; 
      cmd = dss.str();
      dss.str(std::string());
      pm_CmdAddKeyboardEvent (cmd, "arrow_up");

      dss << "binding_simulation " << sim_name << " conformations increment=-1 "; 
      cmd = dss.str();
      dss.str(std::string());
      pm_CmdAddKeyboardEvent (cmd, "arrow_down");
      }

    if (dlist.getInt("level", level)) { 
      sim->setLevel(level);
      }

    if (dlist.getInt("increment", inc)) { 
      sim->incConformation(inc);
      }
    else {
      dlist.getBoolean("show", show); 
      sim->displayConfs(show);
      }
    }
  }

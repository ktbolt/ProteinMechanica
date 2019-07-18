
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
//* sys:     p r o t e i n   m o d e l e r    s y s t e m      *
//*============================================================*

#ifndef _SYS_PM_H_
#define _SYS_PM_H_

#include "pm/pm.h"

using namespace std;

namespace ProteinMechanica {

class PmBindingSimulation;
class PmBody;
class PmDbInterface;
class PmForce;
class PmGeometry;
class PmGraphicsInterface;
class PmGrid;
class PmJoint;
class PmMeasurement;
class PmModel;
class PmMolecule;
class PmMotor;
class PmParticle;
class PmPotential;
class PmPotentialGeom;
class PmSimulation;
class PmSolid;
class PmSurface;
class PmTrace;

class PM_EXPORT PmSystem {
  public:
    PmSystem();
    void enableGraphics(bool flag) { use_graphics = flag; } 
    int useGraphics();  

    void addBody (PmBody *body);
    void getBody (const string name, PmBody **body);
    void getBodies (vector<PmBody*>& bodies);

    bool getCmdEcho() { return cmd_echo; }
    void setCmdEcho(const bool val) { cmd_echo = val; }
    bool getCmdError() { return cmd_error; }
    void setCmdError(const bool flag) { cmd_error = flag; }
    void addCmdVariable(const string name, const string value);
    bool getCmdVerbose() { return cmd_verbose; }
    void setCmdVerbose(const bool val) { cmd_verbose = val; }
    void setSilent(const bool val) { silent = val; }
    bool getSilent(bool val) { return silent; }

    FILE* getCmdFilePtr();
    void setCmdFilePtr(FILE *fp);
    bool getCmdWait() { return (cmd_wait); }
    void setCmdWait(const bool val);

    float getVariableValue(const string name);

    void init(const bool silent);

    void getGraphicsContext (void **context);
    void forceUpdateGraphics();
    void updateGraphics();
    void updateGraphics(const bool flag);
    PmGraphicsInterface* getGraphicsInterace();

    void addDatabase(PmDbInterface *db);
    void getDatabase(const string name, PmDbInterface **p_db);

    void getDomains(vector<PmMolecule*>& domains);
    void addDomain (PmMolecule *domain);
    void getDomain (const string name, PmMolecule **domain);
    void getDomain (PmMolecule **mol) { *mol = this->current_domain; }
    void setCurrentDomain (PmMolecule *mol) { this->current_domain = mol; }

    void procEvents (const char *script);
    void setExtent (PmExtent& extent); 
    void updateExtent(PmExtent& extent); 

    void addForce(PmForce *force);
    void getForce(const string name, PmForce **force);

    void addGeometry(PmGeometry *geom);
    void getGeometry(const string name, PmGeometry **geom);

    void addGrid(PmGrid *grid);
    void getGrid(const string name, PmGrid **grid);

    void procInterrupt();

    void addJoint (PmJoint *joint);
    void getJoint (const string name, PmJoint **joint);
    void getJoints (vector<PmJoint*>& joints);

    void addMeasurement(PmMeasurement *msr);
    void getMeasurement(const string name, PmMeasurement **msr); 

    void addModel (PmModel *model);
    void getModel (const string name, PmModel **model);
    void getModel (PmModel **model) { *model = this->current_model; }
    void setCurrentModel (PmModel *model) { this->current_model = model; }

    void addMolecule (PmMolecule *mol);
    void getMolecule (const string name, PmMolecule **mol);
    void getMolecule (PmMolecule **mol) { *mol = this->current_molecule; }
    void setCurrentMolecule (PmMolecule *mol) { this->current_molecule = mol; }
    void getCurrentMolecule (PmMolecule **mol) { *mol = this->current_molecule; }

    void addMotor(PmMotor *motor);
    void getMotor(const string name, PmMotor **motor);

    void addParticle(PmParticle *part);
    void getParticle(const string name, PmParticle **part);

    void addPotential(PmPotential *pot);
    void getPotential(const string name, PmPotential **pot);

    void addPotentialGeometry(PmPotentialGeom *pot);
    void getPotentialGeometry(const string name, PmPotentialGeom **pot);

    void procQuit();
    void procQuery(PmQuery& query);

    void addBindSimulation (PmBindingSimulation *sim);
    void getBindSimulation (const string name, PmBindingSimulation **sim);

    void addSimulation (PmSimulation *sim);
    void getSimulation (const string name, PmSimulation **sim);
    void getSimulations (vector<PmSimulation*> simulations);
    void replaySimulations();
    void setReplaySimulations(const bool flag, const int start, const int end,
                              const int inc);

    void stepSimulations();
    void setSimulationsNumInc(const int num_inc);
    void setSimulationsStateSaveFreq(const int save_freq);

    void addSolid(PmSolid *solid);
    void getSolid(const string name, PmSolid **solid);

    void addSurface (PmSurface *surf);
    void getSurface (const string name, PmSurface **surf);
    void setCurrentSurface (PmSurface *surf) { this->current_surface = surf; }

    void addTrace (PmTrace *trace);
    void getTrace (const string name, PmTrace **trace);

    void setWorkingDir(const string val);
    void getWorkingDir(string& val);

    void deriveUnits();
    void printUnits();
    float getUnitsMassScale() { return units.mass.scale; }
    void convUnitType(const string str, PmUnitType& type); 
    void setUnitsScale (PmUnitType type, const string name, float scale);

  private:

    // cmd i/o, state,
    int num_cmd_fp;
    vector<FILE*> cmd_fp;
    bool cmd_wait;
    bool cmd_echo;
    bool cmd_verbose;
    bool cmd_error;
    bool silent;

    // working directory //
    string working_directory;

    // database //
    vector<PmDbInterface*> databases;

    // graphics
    PmGraphicsInterface *graphics;
    bool use_graphics;
    bool update_graphics;

    vector<PmBody*> bodies;
    vector<PmJoint*> joints;
    vector<PmForce*> forces;
    vector<PmGeometry*> geometries;
    vector<PmMotor*> motors;
    vector<PmPotential*> potentials;
    vector<PmPotentialGeom*> potential_geometries;
    vector<PmTrace*> traces;

    vector<PmGrid*> grids;
    PmGrid *current_grid;

    vector<PmMolecule*> domains;
    PmMolecule *current_domain;

    vector<PmMeasurement*> measurements;

    vector<PmModel*> models;
    PmModel *current_model;

    vector<PmMolecule*> molecules;
    PmMolecule *current_molecule;

    vector<PmParticle *> particles;
    PmParticle *current_particle;

    vector<PmBindingSimulation*> binding_simulations;

    vector<PmSimulation*> simulations;
    PmSimulation *current_simulation;
    bool replay_simulations;

    vector<PmSolid*> solids;
    PmSolid *current_solid;

    vector<PmSurface*> surfaces;
    PmSurface *current_surface;

    bool extent_set;
    PmExtent extent;

    PmUnits units;
  };

extern PM_EXPORT PmSystem pmSystem;

}

#endif



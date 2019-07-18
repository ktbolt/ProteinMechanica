
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

/*============================================================*
 * sobj:          s i m u l a t i o n    o b j e c t          *
 *============================================================*/

#ifndef _SIMULATION_OBJECT_PM_H_
#define _SIMULATION_OBJECT_PM_H_

#include "pm/pm.h"
#include "pobj.h"
#include "mol.h"
#include "pot.h"
#include "force.h"

namespace ProteinMechanica {


// PmSimulationObj
// ---------------

 class PM_EXPORT PmSimulationObj {

    public:
       PmSimulationObj();
       ~PmSimulationObj();
       void addForce(PmVector3& pos, PmVector3& dir);
       void addForce(int index, PmVector3& dir);
       void resetForces() { forces.clear(); torques.clear(); }
       void getForces(vector<PmForceVector>& forces) { forces = this->forces; }

       void addGeometry(PmGeometry *geom);
       void getGeometries(vector<PmGeometry*>& geoms) { geoms = this->geometries; }

       void getMolecule(PmMolecule **mol);
       void getName (string& name) { name = this->name;  }
       void addXformPhysicalObj (PmPhysicalObj *pobj);
       void getXformPhysicalObjs (vector<PmPhysicalObj*>& pobjs);
       void setPhysicalObj (PmPhysicalObj *pobj);
       void getPhysicalObj (PmPhysicalObj **pobj) { *pobj = this->pobj; }
       bool hasPhysicalObj() { return this->pobj != NULL; }
       void addPotentialGeometry(PmPotentialGeom *pot);
       void displayPotentialGeometry(const string pname, const bool show);
       void getPotential(const string name, PmPotential **pot);
       bool hasPotential(const string name);
       void getPotentialGeometry(const string pname, PmPotentialGeom **pgeom);
       bool hasPotentialGeometry(const string pname);
       bool hasRegion(const string rgn_name);
       void getRegionCoords(const string rgn_name, vector<PmVector3>& coords);
       void getRegionData(const string rgn_name, vector<float>& data);
       void getRegionIndices(const string rgn_name, vector<int>& ids);
       void getRegionRadii(const string rgn_name, vector<float>& rads);
       void addTorque(PmVector3& pos, PmVector3& dir);
       void getTorques(vector<PmForceVector>& torques) { torques = this->torques; }

    protected:
       string name;
       PmPhysicalObj *pobj;
       vector<PmPhysicalObj*> xform_pobjs;
       vector<PmPotentialGeom*> potential_geometries;
       vector<PmForceVector> forces;
       vector<PmForceVector> torques;
       vector<PmGeometry*> geometries;
   };

}

#endif


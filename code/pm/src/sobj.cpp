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
//* sobj:          s i m u l a t i o n    o b j e c t          *
//*============================================================*

#include "sobj.h"

namespace ProteinMechanica {


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSimulationObj::PmSimulationObj() {
  }

PmSimulationObj::~PmSimulationObj() {
  }

//*============================================================*
//*==========            setPhysicalObj              ==========*
//*============================================================*

void 
PmSimulationObj::setPhysicalObj (PmPhysicalObj *pobj) 
  { 
  this->pobj = pobj; 

  if (pobj) {
    pobj->setSimulationObj(this);
    }
  }

//*============================================================*
//*==========            addXformPhysicalObj         ==========*
//*============================================================*

void
PmSimulationObj::addXformPhysicalObj (PmPhysicalObj *pobj)
  {
  xform_pobjs.push_back(pobj);
  }

//*============================================================*
//*==========            getXformPhysicalObjs        ==========*
//*============================================================*

void
PmSimulationObj::getXformPhysicalObjs (vector<PmPhysicalObj*>& pobjs)
  {
  pobjs = xform_pobjs;
  }

//*============================================================*
//*==========            addForce                    ==========*
//*============================================================*

void 
PmSimulationObj::addForce(PmVector3& pos, PmVector3& dir)
  {
  /*
  fprintf (stderr, ">>>>>> PmSimulationObj::addForce \n"); 
  fprintf (stderr, "   >>> pos (%g %g %g) \n", pos[0], pos[1], pos[2]); 
  fprintf (stderr, "   >>> dir (%g %g %g) \n", dir[0], dir[1], dir[2]); 
  */
  PmForceVector force;
  force.point = pos;
  force.direction = dir;
  forces.push_back(force);
  }

//*============================================================*
//*==========            addForce                    ==========*
//*============================================================*

void 
PmSimulationObj::addForce(int index, PmVector3& dir)
  {
  PmForceVector force;
  force.index = index;
  force.direction = dir;
  forces.push_back(force);
  }

//*============================================================*
//*==========            addGeometry                 ==========*
//*============================================================*

void 
PmSimulationObj::addGeometry(PmGeometry *geom)
  {
  geometries.push_back(geom);
  }

//*============================================================*
//*==========            getMolecule                 ==========*
//*============================================================*

void 
PmSimulationObj::getMolecule(PmMolecule **mol)
  {
  *mol = dynamic_cast<PmMolecule*>(pobj);
  }

//*============================================================*
//*==========            displayPotentialGeometry    ==========*
//*============================================================*

void
PmSimulationObj::displayPotentialGeometry(const string pname, const bool show)
  {
  PmPotentialGeom *pgeom;
  getPotentialGeometry(pname, &pgeom);

  if (!pgeom) {
    return;
    }

  pgeom->display(show);
  }

//*============================================================*
//*==========            getPotentialGeometry        ==========*
//*============================================================*

void
PmSimulationObj::getPotentialGeometry(const string pname, PmPotentialGeom **pot)
  {
  *pot = NULL;

  for (unsigned int i = 0; i < potential_geometries.size(); i++) {
    *pot = potential_geometries[i];
    if ((*pot)->hasName(pname)) {
      return;
      }
    }
  }

//*============================================================*
//*==========            hasPotentialGeometry        ==========*
//*============================================================*

bool 
PmSimulationObj::hasPotentialGeometry(const string pname)
  {
  for (unsigned int i = 0; i < potential_geometries.size(); i++) {
    PmPotentialGeom *pot = potential_geometries[i];
    if (pot->hasName(pname)) {
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========            getRegionCoords             ==========*
//*============================================================*

void 
PmSimulationObj::getRegionCoords(const string rgn_name, vector<PmVector3>& coords)
  {
  //fprintf (stderr, ">>>>>> PmSimulationObj::getRegionCoords  \n"); 
  coords.clear();

  if (!pobj) {
    //fprintf (stderr, "   >>> no pobj. "); 
    return;
    }

  pobj->getRegionCoords(rgn_name, coords);
  }

//*============================================================*
//*==========            getRegionData               ==========*
//*============================================================*

void
PmSimulationObj::getRegionData(const string rgn_name, vector<float>& data)
  {
  data.clear();

  if (!pobj) {
    return;
    }

  pobj->getRegionData(rgn_name, data);
  }

//*============================================================*
//*==========            getRegionRadii              ==========*
//*============================================================*

void
PmSimulationObj::getRegionRadii(const string rgn_name, vector<float>& rads)
  {
  rads.clear();

  if (!pobj) {
    return;
    }

  pobj->getRegionRadii(rgn_name, rads);
  }

//*============================================================*
//*==========            getRegionIndices            ==========*
//*============================================================*

void
PmSimulationObj::getRegionIndices(const string rgn_name, vector<int>& ids)
  {
  ids.clear();

  if (!pobj) {
    return;
    }

  pobj->getRegionIndices(rgn_name, ids);
  }

//*============================================================*
//*==========            hasRegion                   ==========*
//*============================================================*

bool 
PmSimulationObj::hasRegion(const string rgn_name)
  {
  if (!pobj) {
    return false;
    }
  return pobj->hasRegion(rgn_name);
  }

//*============================================================*
//*==========            addPotentialGeometry        ==========*
//*============================================================*

void 
PmSimulationObj::addPotentialGeometry(PmPotentialGeom *pot) 
  {
  potential_geometries.push_back(pot);
  }

//*============================================================*
//*==========            addTorque                   ==========*
//*============================================================*

void 
PmSimulationObj::addTorque(PmVector3& pos, PmVector3& dir)
  {
  /*
  fprintf (stderr, ">>>>>> PmSimulationObj::addTorque\n"); 
  fprintf (stderr, "   >>> pos (%g %g %g) \n", pos[0], pos[1], pos[2]); 
  fprintf (stderr, "   >>> dir (%g %g %g) \n", dir[0], dir[1], dir[2]); 
  */
  PmForceVector force;
  force.point = pos;
  force.direction = dir;
  torques.push_back(force);
  }

//*============================================================*
//*==========            addForce                    ==========*
//*============================================================*


}



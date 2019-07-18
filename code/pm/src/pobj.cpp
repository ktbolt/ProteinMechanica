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
//* pobj:             p h y s i c a l   o b j e c t            *
//*============================================================*

#include "pobj.h"
#include "mol.h"
#include "part.h"
#include "solid.h"
#include "surf.h"

#define ndbg_PmPhysicalObj_createRegions

namespace ProteinMechanica {

vector<PmPhysicalObj*> PmPhysicalObj::physical_objects;


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPhysicalObj::PmPhysicalObj() 
  {
  PmPhysicalObj::physical_objects.push_back(this); 
  }

PmPhysicalObj::~PmPhysicalObj() {
  }

//*============================================================*
//*==========             getName                    ==========*
//*============================================================*

bool 
PmPhysicalObj::get(string& name, PmPhysicalObj **pobj) 
  {
  string pname;
  *pobj = NULL;

  for (unsigned int i = 0; i < physical_objects.size(); i++) {
    physical_objects[i]->getName(pname);

    if (pname == name) {
      *pobj = physical_objects[i];
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========             getName                    ==========*
//*============================================================*

void 
PmPhysicalObj::getName(string& name)
  {
  PmMolecule *mol;
  PmParticle *part;
  PmSurface *surf;
  PmSolid *solid;

  switch (obj_type) {
    case PM_OBJECT_UNKNOWN:
      name = "unknown";
    break;

    case PM_OBJECT_DOMAIN:
      mol = dynamic_cast<PmMolecule*>(this);
      mol->getName(name);
    break;

    case PM_OBJECT_MODEL:
      name = "unknown";
    break;

    case PM_OBJECT_MOLECULE:
      mol = dynamic_cast<PmMolecule*>(this);
      mol->getName(name);
    break;

    case PM_OBJECT_PARTICLE:
      part = dynamic_cast<PmParticle*>(this);
      part->getName(name);
    break;

    case PM_OBJECT_SOLID:
      solid = dynamic_cast<PmSolid*>(this);
      solid->getName(name);
    break;

    case PM_OBJECT_SURFACE:
      surf = dynamic_cast<PmSurface*>(this);
      surf->getName(name);
    break;
    }
  }

//*============================================================*
//*==========             createHbondRegions         ==========*
//*============================================================*
// create a pair of regions for the hydrogen bonds between two 
// molecular physical objects.

void
PmPhysicalObj::createHbondRegions(const string rgn_prefix, PmPhysicalObj *pobj,
                                  const float cutoff, const string chain,
                                  vector<PmRegion*>& rgns)
  {
  #define ndbg_PmPhysicalObj_createHbondRegions
  #ifdef dbg_PmPhysicalObj_createHbondRegions
  fprintf (stderr, "\n>>>>>> PmPhysicalObj::createHbondRegions \n");
  fprintf (stderr, "   >>> cutoff = %f \n", cutoff);
  #endif
                            
  PmMolecule *mol1, *mol2;
  char chain_id;
  int num_res1, num_res2;
  vector<PmResidue*> rlist1, rlist2;
  PmResidue *res1, *res2; 
  int n1;
  PmAtom *n_atom, *o_atom;
  vector<int> ids1, ids2;
  //vector<PmVector3> coords1, coords2;
  string rgn_name, name1, name2;
  stringstream dss;
  PmRegion *rgn1, *rgn2;

  rgns.clear();

  mol1 = dynamic_cast<PmMolecule*>(this);
  mol2 = dynamic_cast<PmMolecule*>(pobj);

  if (!mol1 || !mol2) {
    return;
    }

  // get the residues for the given chain //

  chain_id = chain[0];
  mol1->getChainResidues(chain_id, rlist1);
  mol2->getChainResidues(chain_id, rlist2);

  num_res1 = rlist1.size();
  num_res2 = rlist2.size();

  if (!num_res1 || !num_res2) {
    return;
    }

  #ifdef dbg_PmPhysicalObj_createHbondRegions
  fprintf (stderr, "   >>> num res1 = %d \n", num_res1);
  fprintf (stderr, "   >>> num res2 = %d \n", num_res2);
  #endif

  res1 = rlist1[num_res1-1];
  res2 = rlist2[0];

  #ifdef dbg_PmPhysicalObj_createHbondRegions
  fprintf (stderr, "   >>> res1 id = %d \n", res1->id);
  fprintf (stderr, "   >>> res2 id = %d \n", res2->id);
  #endif

  // find residue in mol1 that is 4 residues   //
  // away from first residue in mo2.           //

  n1 = -1;

  for (int i = 0; i < num_res1; i++) {
    if ((res2->id - rlist1[i]->id) == 4) {
      n1 = i;
      break;
      }
    }

  if (n1 == -1) {
    return;
    }

  for (int i = n1, n2 = 0; i < num_res1 && n2 < num_res2; i++, n2++) {
    res1 = rlist1[i];
    res2 = rlist2[n2];
    #ifdef dbg_PmPhysicalObj_createHbondRegions
    fprintf (stderr, "   >>> res1 = %d   res2 = %d    ", res1->id, res2->id);
    #endif

    res1->getAtom("O", &o_atom);
    res2->getAtom("N", &n_atom);

    if (o_atom && n_atom) {
      ids1.push_back(o_atom->index);
      ids2.push_back(n_atom->index);
      #ifdef dbg_PmPhysicalObj_createHbondRegions
      fprintf (stderr, "   has bond.  O:%d N:%d \n", o_atom->index, n_atom->index);
      #endif
      }
    else {
      #ifdef dbg_PmPhysicalObj_createHbondRegions
      fprintf (stderr, "   NO BOND! \n");
      #endif
      }
    }

  if (ids1.size() == 0) { 
    return;
    }

  // create regions for each physical object //

  this->getName(name1);
  pobj->getName(name2);
  name1[0] = toupper(name1[0]);
  name2[0] = toupper(name2[0]);

  dss << rgn_prefix << name1 << '_' << name2 << "Rgn";
  rgn_name = dss.str();
  dss.str(std::string());
  this->defineRegion(rgn_name, ids1);
  this->getRegion(rgn_name, &rgn1);

  dss << rgn_prefix << name2 << '_' << name1 << "Rgn";
  rgn_name = dss.str();
  dss.str(std::string());
  pobj->defineRegion(rgn_name, ids2);
  pobj->getRegion(rgn_name, &rgn2);

  rgns.push_back(rgn1);
  rgns.push_back(rgn2);
  }

//*============================================================*
//*==========             createRegions              ==========*
//*============================================================*
// create a pair of regions between two physical objects.

void 
PmPhysicalObj::createRegions(const string rgn_prefix, PmPhysicalObj *pobj, 
                             const float cutoff, bool use_sidechains, 
                             PmAtomFilter& filter, vector<PmRegion*>& rgns)
  {
  #define dbg_PmPhysicalObj_createRegions
  #ifdef dbg_PmPhysicalObj_createRegions
  fprintf (stderr, "\n>>>>>> PmPhysicalObj::createRegions \n");
  fprintf (stderr, ">>> cutoff=%f \n", cutoff);
  fprintf (stderr, ">>> filter names size=%d \n", filter.names.size());
  fprintf (stderr, ">>> filter exclude=%d \n", filter.exclude);
  fprintf (stderr, ">>> use_sidechains=%d \n", use_sidechains); 
  #endif

  int n1, n2;
  vector<PmVector3> coords1, coords2;
  float r1, r2, d;
  vector<float> rads1, rads2;
  PmVector3 pt1, pt2, v; 
  int num_i1, num_i2;
  bool *ilist1, *ilist2;
  vector<int> ids1, ids2;
  vector<string> atom_names1, atom_names2;
  string rgn_name, name1, name2;
  stringstream dss;
  PmRegion *rgn1, *rgn2;
  PmMolecule *mol1, *mol2;
  PmMolRegionParameters mol_params;

  mol1 = mol2 = NULL;
  rgns.clear();

  //===== use centers of side chain atoms =====//

  if (use_sidechains) {
    mol1 = dynamic_cast<PmMolecule*>(this);
    mol2 = dynamic_cast<PmMolecule*>(pobj);

    if (!filter.sidechain) {
      if (mol1) {
        mol1->getSidechainGeometry("", coords1, rads1);
        #ifdef dbg_PmPhysicalObj_createRegions
        fprintf (stderr, ">>> use sidechains for pobj1 \n");
        #endif
        mol_params.use_sidechains = true;
        mol1->getBackboneAtomNames(filter.names);
        filter.exclude = true;
        }

      if (mol2) {
        mol2->getSidechainGeometry("", coords2, rads2);
        #ifdef dbg_PmPhysicalObj_createRegions
        fprintf (stderr, ">>> use sidechains for pobj2 \n");
        #endif
        }
      }

    else {
      #ifdef dbg_PmPhysicalObj_createRegions
      fprintf (stderr, ">>> use filter backbone names \n");
      #endif
      mol1->getBackboneAtomNames(filter.names);
      mol1->getAtomCoords("", filter, coords1, rads1);
      mol1->getAtomNames("", filter, atom_names1);
      mol2->getAtomCoords("", filter, coords2, rads2);
      mol2->getAtomNames("", filter, atom_names2);
      }
    }

  //===== use either side chain or main chain atoms =====//

  else if ((filter.sidechain) || (filter.mainchain)) {
    #ifdef dbg_PmPhysicalObj_createRegions
    fprintf (stderr, ">>> filter.sidechain=%d \n", filter.sidechain);
    fprintf (stderr, ">>> filter.mainchain=%d \n", filter.mainchain);
    #endif
    mol1 = dynamic_cast<PmMolecule*>(this);
    mol2 = dynamic_cast<PmMolecule*>(pobj);
    filter.exclude = filter.sidechain;
    mol1->getBackboneAtomNames(filter.names);
    mol1->getAtomCoords("", filter, coords1, rads1);
    mol1->getAtomNames("", filter, atom_names1);
    mol2->getAtomCoords("", filter, coords2, rads2);
    mol2->getAtomNames("", filter, atom_names2);
    }

  //===== use atom names in filter =====//

  else if (filter.names.size()) {
    mol1 = dynamic_cast<PmMolecule*>(this);
    mol2 = dynamic_cast<PmMolecule*>(pobj);

    // NOTE: need to check for atom radius? //

    if (mol1) {
      mol1->getAtomCoords("", filter, coords1, rads1);
      mol1->getAtomNames("", filter, atom_names1);
      }

    if (mol2) {
      mol2->getAtomCoords("", filter, coords2, rads2);
      mol2->getAtomNames("", filter, atom_names2);
      }
    }

  //===== use all the atoms =====//

  else {
    this->getCoordinates(coords1);
    this->getRadii(rads1);
    pobj->getCoordinates(coords2);
    pobj->getRadii(rads2);
    }

  n1 = coords1.size();
  n2 = coords2.size();

  if (!n1 || !n2) {
    return;
    }

  #ifdef dbg_PmPhysicalObj_createRegions
  fprintf (stderr, "   >>> num coords1=%d \n", n1); 
  fprintf (stderr, "   >>> num rads1=%d \n", rads1.size()); 
  fprintf (stderr, "   >>> num atom_names1=%d \n", atom_names1.size()); 
  fprintf (stderr, "   >>> num coords2=%d \n", n2); 
  fprintf (stderr, "   >>> num rads2=%d \n", rads2.size()); 
  fprintf (stderr, "   >>> num atom_names2=%d \n", atom_names2.size()); 
  fprintf (stderr, "   >>> comp int ... ");
  #endif

  ilist1 = new bool[n1];

  for (int i = 0; i < n1; i++) {
    ilist1[i] = false; 
    }

  ilist2 = new bool[n2];

  for (int i = 0; i < n2; i++) {
    ilist2[i] = false; 
    }

  // determine coordinates that are within cutoff //

  int m = 0;

  for (int i = 0; i < n1; i++) {
    pt1 = coords1[i];
    r1 = rads1[i];

    for (int j = 0; j < n2; j++) {
      pt2 = coords2[j];
      r2 = rads2[j];
      v = pt1 - pt2;
      d = v.length() - r1 - r2;
      //d = v.length();

      if (d <= cutoff) {
        ilist1[i] = true; 
        ilist2[j] = true;
        m += 1;
        //fprintf (stderr, " %d: %d %d  r(%f,%f) d = %g \n", m, i, j, r1, r2, d);
        }
      }
    }

  num_i1 = 0;
  num_i2 = 0;

  for (int i = 0; i < n1; i++) {
    if (ilist1[i]) {
      ids1.push_back(i+1);
      num_i1 += 1; 
      }
    }

  for (int i = 0; i < n2; i++) {
    if (ilist2[i]) {
      ids2.push_back(i+1);
      num_i2 += 1; 
      }
    }

  delete [] ilist1;
  delete [] ilist2;

  #ifdef dbg_PmPhysicalObj_createRegions
  fprintf (stderr, "done. \n\n");
  fprintf (stderr, "   >>> num int1=%d \n", num_i1); 
  fprintf (stderr, "   >>> num ids1=%d \n", ids1.size()); 
  fprintf (stderr, "\n");
  fprintf (stderr, "   >>> num int2=%d \n", num_i2);
  fprintf (stderr, "   >>> num ids2=%d \n", ids2.size()); 
  #endif

  if (!num_i1 || !num_i2) {
    return;
    }

  //===== create regions for each physical object ... =====//

  #ifdef dbg_PmPhysicalObj_createRegions
  fprintf (stderr, ">>> create regions \n");
  #endif

  this->getName(name1);
  pobj->getName(name2);
  name1[0] = toupper(name1[0]); 
  name2[0] = toupper(name2[0]); 

  dss << rgn_prefix << name1 << '_' << name2 << "Rgn";
  rgn_name = dss.str();
  dss.str(std::string());

  // 1st region //

  if (mol1) {
    mol1->defineRegion(rgn_name, ids1, rads1, coords1, atom_names1); 
    }
  else {
    this->defineRegion(rgn_name, ids1);
    }

  this->getRegion(rgn_name, &rgn1);

  // 2nd region //

  dss << rgn_prefix << name2 << '_' << name1 << "Rgn";
  rgn_name = dss.str();
  dss.str(std::string());

  if (mol2) {
    mol2->defineRegion(rgn_name, ids2, rads2, coords2, atom_names2); 
    }
  else {
    pobj->defineRegion(rgn_name, ids2);
    }

  pobj->getRegion(rgn_name, &rgn2);

  // return regions //

  rgns.push_back(rgn1);
  rgns.push_back(rgn2);
  }

//*============================================================*
//*==========             getRegion                  ==========*
//*============================================================*

void 
PmPhysicalObj::getRegion(const string rgn_name, PmRegion **rgn)
  {
  *rgn = NULL;

  for (unsigned int i = 0; i < regions.size(); i++) {
    if (regions[i]->name == rgn_name) {
      *rgn = regions[i];
      return;
      }
    }
  }

//*============================================================*
//*==========             hasRegion                  ==========*
//*============================================================*

bool 
PmPhysicalObj::hasRegion(const string rgn_name)
  {
  /*
  fprintf (stderr, "\n>>>>>> PmPhysicalObj::hasRegion \n");
  fprintf (stderr, ">>> rgn_name=%s \n", rgn_name.c_str());
  fprintf (stderr, ">>> num rgns=%d \n", regions.size()); 
  */
 
  for (unsigned int i = 0; i < regions.size(); i++) {
    if (regions[i]->name == rgn_name) {
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========             getRegionCenter            ==========*
//*============================================================*

void
PmPhysicalObj::getRegionCenter(const string rgn_name, PmVector3& center)
  {
  PmRegion *rgn;
  getRegion(rgn_name, &rgn);

  if (!rgn) {
    return;
    }

  center = rgn->center;
  }

//*============================================================*
//*==========             getRegionCoords            ==========*
//*============================================================*

void
PmPhysicalObj::getRegionCoords(const string rgn_name, vector<PmVector3>& coords)
  {
  PmRegion *rgn;
  coords.clear();
  getRegion(rgn_name, &rgn);

  if (!rgn) {
    return;
    }

  coords = rgn->coords;
  }

//*============================================================*
//*==========             getRegionData              ==========*
//*============================================================*

void
PmPhysicalObj::getRegionData(const string rgn_name, vector<float>& data)
  {
  PmRegion *rgn;
  data.clear();
  getRegion(rgn_name, &rgn);

  if (!rgn) {
    return;
    }

  data = rgn->data;
  }

//*============================================================*
//*==========             getRegionIndices           ==========*
//*============================================================*

void
PmPhysicalObj::getRegionIndices(const string rgn_name, vector<int>& ids)
  {
  PmRegion *rgn;
  ids.clear();
  getRegion(rgn_name, &rgn);

  if (!rgn) {
    return;
    }

  ids = rgn->coord_indexes;
  }

//*============================================================*
//*==========             getRegionRadii             ==========*
//*============================================================*

void
PmPhysicalObj::getRegionRadii(const string rgn_name, vector<float>& rads)
  {
  PmRegion *rgn;
  rads.clear();
  getRegion(rgn_name, &rgn);

  if (!rgn) {
    return;
    }

  rads = rgn->radii;
  }

//*============================================================*
//*==========             get/setSimulationObj       ==========*
//*============================================================*

void 
PmPhysicalObj::setSimulationObj(PmSimulationObj *sobj) {
  this->sobj = sobj;
  }

void 
PmPhysicalObj::getSimulationObj(PmSimulationObj **sobj) {
  *sobj = this->sobj;
  }

//*============================================================*
//*==========             gettype                    ==========*
//*============================================================*

void
PmPhysicalObj::getObjectType(string& str)
  {

  switch (obj_type) {
    case PM_OBJECT_UNKNOWN:
      str = "unknown";
    break;

    case PM_OBJECT_DOMAIN:
      str = "domain";
    break;

    case PM_OBJECT_MODEL:
      str = "model";
    break;

    case PM_OBJECT_MOLECULE:
      str = "molecule";
    break;

    case PM_OBJECT_PARTICLE:
      str = "particle";
    break;

    case PM_OBJECT_SOLID:
      str = "solid";
    break;

    case PM_OBJECT_SURFACE:
      str = "surface";
    break;
    }
  }


}



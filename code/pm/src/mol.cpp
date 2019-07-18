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
//* mol:                 m o l e c u l e                       *
//*============================================================*
// molecule object.

#include "mol.h"

namespace ProteinMechanica {

#define ndbg_PmMolecule
#define ndbg_getChainResidues 
#define ndbg_createDomain 
#define ndbg_displayBackbone
#define ndbg_getResidueCoords 

// for a single residue add an extra atom 
#define ncreateDomain_add_atom

// bond lengths and angles //
/*
static float c_o_l   = 0.123,   c_o_ang = 121;
static float c_n_l   = 0.133,   c_n_ang = 123;
static float n_ca_l  = 0.145,  n_ca_ang = 122;
static float ca_cb_l = 0.155,  ca_cb_ang = 110;
static float ca_c_l  = 0.152;
*/

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*
// create a molecule from a name.

PmMolecule::PmMolecule(const string name) 
  {
  this->name = name;
  this->obj_type = PM_OBJECT_MOLECULE;
  num_residues = 0;
  num_chains = 0;
  parent = NULL;
  backbone_atom = "CA";
  surface = NULL;
  datts.color.set(1,1,1);
  datts.marker_size = 0.04; 

  mode_number = 6;
  mode_color.set(1,1,1);
  mode_scale = 1.0; 

  topology = NULL;
  residue_contacts = NULL;
  }

PmMolecule::~PmMolecule() {
  }

//*============================================================*
//*==========              copy.                     ==========*
//*============================================================*
// copy a molecule.                   

void 
PmMolecule::copy(const string name, PmMolecule **p_mol)
  {
  PmMolecule *mol = new PmMolecule(name);

  // copy atoms //

  for (int i = 0; i < number_of_atoms; i++) {
    PmAtom *atom = &atoms[i];
    mol->addAtom(atom);
    }

  mol->type = this->type;
  mol->backbone_atom = this->backbone_atom;
  mol->parent = this->parent;
  mol->connect = this->connect;
  *p_mol = mol;
  }

//*============================================================*
//*==========              convMoleculeType          ==========*
//*============================================================*
// convert molecule string type to enum. 

void 
PmMolecule::convMoleculeType(const string str, PmMoleculeType& type)
  {
  type = PM_MOLECULE_UNKNOWN;

  if (str == "protein") {
    type = PM_MOLECULE_PROTEIN;
    }
  else if (str == "miscellaneous") {
    type = PM_MOLECULE_MISCELLANEOUS;
    }
  else if (str == "nucleotide") {
    type = PM_MOLECULE_NUCLEOTIDE;
    }
  else if (str == "peptidoglycan") {
    type = PM_MOLECULE_PEPTIDOGLYCAN;
    }
  else if (str == "water") {
    type = PM_MOLECULE_WATER;
    }
  }

//*============================================================*
//*==========              getAtoms                  ==========*
//*============================================================*
// get the atoms for the molecule.     

void 
PmMolecule::getAtoms(vector<PmAtom>& mol_atoms)
  {
  mol_atoms.clear();

  for (int i = 0; i < number_of_atoms; i++) {
    mol_atoms.push_back (atoms[i]);
    }
  }

void
PmMolecule::getAtoms(PmAtomFilter& filter, vector<PmAtom>& mol_atoms)
  {
  mol_atoms.clear();

  for (int i = 0; i < number_of_atoms; i++) {
    if (filter.hasName(atoms[i].name)) {
      mol_atoms.push_back (atoms[i]);
      }
    }
  }

//*============================================================*
//*==========          getBackboneAtomNames          ==========*
//*============================================================*
// get the backbone atom names for the molecule.

void
PmMolecule::getBackboneAtomNames(vector<string>& names)
  {
  if (type == PM_MOLECULE_NUCLEIC_ACID) {
    names.push_back("OP1");
    names.push_back("OP2");
    names.push_back("P");
    names.push_back("O5'");
    names.push_back("O5'");
    names.push_back("C5'");
    names.push_back("C4'");
    names.push_back("O4'");
    names.push_back("C2'");
    names.push_back("C3'");
    names.push_back("O3'");
    }
  else {
    names.push_back("N");
    names.push_back("CA");
    names.push_back("C");
    names.push_back("O");
    }
  }

//*============================================================*
//*==========              getAtomCoords             ==========*
//*============================================================*
// get the coordinates for the atoms of the given sequence.
// ------------------------------------------------------------
// desc: molecular descriptor
// filter: atom filter

void 
PmMolecule::getAtomCoords(const string desc, PmAtomFilter& filter, 
                          vector<PmVector3>& coords)
  {
  #ifdef dbg_PmMolecule_getAtomCoords
  fprintf (stderr, ">>>>>> PmMolecule::getAtomCoords:  desc = %s \n", desc.c_str());
  #endif
  coords.clear();
  vector<PmResidue*> rlist;

  //===== get the residues from the descriptor =====//

  getResidues(desc, rlist);

  #ifdef dbg_PmMolecule_getAtomCoords
  fprintf (stderr, ">>> num res[%d] \n", rlist.size());
  #endif

  if (!rlist.size()) {
    return;
    }

  //===== get the coordinates for the atoms of the residues =====//

  PmResidue *res;

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];

      if (filter.hasName(atom->name)) {
        coords.push_back(atom->pos);
        //fprintf (stderr, ">>>> atom=%s\n", atom->name);
        }
      }
    }
  }

//*============================================================*
//*==========              getAtomCoords             ==========*
//*============================================================*
// get the coordinates for the atoms of the given sequence.
// ------------------------------------------------------------
// desc: molecular descriptor
// filter: atom filter

void
PmMolecule::getAtomCoords(const string desc, PmAtomFilter& filter,
                          vector<PmVector3>& coords, vector<float>& rads)
  {
  #ifdef dbg_PmMolecule_getAtomCoords
  fprintf (stderr, ">>>>>> PmMolecule::getAtomCoords:  desc = %s \n", desc.c_str());
  #endif
  coords.clear();
  rads.clear();
  vector<PmResidue*> rlist;

  //===== get the residues from the descriptor =====//

  getResidues(desc, rlist);

  #ifdef dbg_PmMolecule_getAtomCoords
  fprintf (stderr, ">>> num res[%d] \n", rlist.size());
  #endif

  if (!rlist.size()) {
    return;
    }

  //===== get the coordinates for the atoms of the residues =====//

  PmResidue *res;

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];

      if (filter.hasName(atom->name)) {
        coords.push_back(atom->pos);
        rads.push_back(atom->getRadius());
        //fprintf (stderr, "#### atom=%s id=%d \n", atom->name, atom->id); 
        }
      }
    }
  }

//*============================================================*
//*==========              getAtomNames              ==========*
//*============================================================*

void
PmMolecule::getAtomNames(const string desc, PmAtomFilter& filter,
                          vector<string>& names)
  {
  #ifdef dbg_PmMolecule_getAtomCoords
  fprintf (stderr, ">>>>>> PmMolecule::getAtomCoords:  desc = %s \n", desc.c_str());
  #endif
  names.clear();
  vector<PmResidue*> rlist;

  //===== get the residues from the descriptor =====//

  getResidues(desc, rlist);

  #ifdef dbg_PmMolecule_getAtomCoords
  fprintf (stderr, ">>> num res[%d] \n", rlist.size());
  #endif

  if (!rlist.size()) {
    return;
    }

  //===== get the coordinates for the atoms of the residues =====//

  PmResidue *res;

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];

      if (filter.hasName(atom->name)) {
        names.push_back(atom->name);
        }
      }
    }
  }

//*============================================================*
//*==========              xformAtoms                ==========*
//*============================================================*
// modify the atom coordinates of the molecule using the given
// transformation.

void 
PmMolecule::xformAtoms(PmXform& xform)
  {
  PmMatrix3x3 mat = xform.matrix;
  PmVector3 pos, xpos;

  for (int i = 0; i < number_of_atoms; i++) {
    pos = atoms[i].pos;
    xpos = mat*(pos-xform.center) + xform.translation + xform.center;
    atoms[i].pos = xpos;
    }
  }

//*============================================================*
//*==========              getAtoms                  ==========*
//*============================================================*
// get the two atoms described in <atoms_desc>.
// atoms_desc = { resName atomName resName atomName }

bool 
PmMolecule::getAtoms(vector<string>atoms_desc, PmAtom& atom1, PmAtom& atom2)
  {
  if (atoms_desc.size() != 4) {
    return false;
    }

  vector<PmResidue*> rlist;

  for (int i = 0; i < 2; i++) {
    string desc = atoms_desc[2*i];
    getResidues(desc, rlist);

    if (!rlist.size()) {
      return false;
      }

    PmResidue *res = rlist[0];
    PmAtom *atom;
    res->getAtom(atoms_desc[2*i+1], &atom);

    if (!atom) {
      return false;
      }
    else if (i == 0) {
      atom1 = *atom;
      }
    else if (i == 1) {
      atom2 = *atom;
      }
    }

  return true;
  }

//*============================================================*
//*==========              getCenter                 ==========*
//*============================================================*
// get the centroid for a group of residues. 

void
PmMolecule::getCenter(const string desc, PmAtomFilter& filter, PmVector3& center)
  {
  vector<PmResidue*> rlist;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    return;
    }

  PmResidue *res;
  center.set(0,0,0);
  int num = 0;

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];
    int n = res->num_atoms;
    num += n;

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];
      center = center + atom->pos;
      }
    }

  if (num != 0) {
    center = (1.0/num)*center;
    }
  }

//*============================================================*
//*==========              check                     ==========*
//*============================================================*
// check a molecule for a non-continguous backbone. 

void 
PmMolecule::check()
  {
  if (!num_chains) buildChains(); 
  
  int num_coords; 
  PmVector3 *coords, v;
  int n1, num_breaks;
  float tol = 0.5, d;
  vector<int> ends;
  vector<PmResidue*> rlist; 
  PmResidue *res1, *res2;

  // check each chain in the molecule //

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    fprintf (stderr, "---- chain %c ---- \n", chain->id);
    getChainAtomCoords (chain->id, backbone_atom, &num_coords, &coords);
    getChainResidues (chain->id, rlist);

    // test for non-contiguous backbone //

    if (num_coords) {
      num_breaks = 0; 
      n1 = 0;

      for (int j = 0; j < num_coords-1; j++) {
        v = coords[j] - coords[j+1];
        d = v.length();
 
        if (d > tol) {
          res1 = rlist[n1];
          res2 = rlist[j];
          num_breaks += 1; 
          fprintf (stderr, "fragment %d:  %d - %d \n", num_breaks, res1->id, res2->id);
          n1 = j+1;
          }
        }

      if (num_breaks) {
        res1 = rlist[n1];
        res2 = rlist[num_coords-1];
        num_breaks += 1; 
        fprintf (stderr, "fragment %d:  %d - %d \n", num_breaks, res1->id, res2->id);
        }
      else {
        fprintf (stderr, "contiguous chain \n");
        }
      }
    }
  }

//*============================================================*
//*==========              print                     ==========*
//*============================================================*
// print the information for a molecule.

void
PmMolecule::print()
  {
  if (!num_chains) buildChains();

  vector<PmAtom> atoms;
  string name, res_desc;
  vector<PmResidue*> rlist;
  char chain_id;
  int id, n, num_res;

  this->getName(name);
  this->getAtoms(atoms);
  this->getResidues(res_desc, rlist);

  fprintf (stderr, "\n---------- \"%s\" ---------- \n", name.c_str());
  fprintf (stderr, "number of atoms = %d \n", atoms.size());
  fprintf (stderr, "number of residues = %d \n", rlist.size());
  fprintf (stderr, "\n------ chain:residues ------\n");

  num_res = rlist.size();
  chain_id = rlist[0]->atoms[0]->chain;
  id = rlist[0]->id;
  n = 1;
  fprintf (stderr, "%c: ", chain_id);
  fprintf (stderr, "%d", id);

  for (int i = 0; i < num_res-1; i++) {
    if (rlist[i]->atoms[0]->chain != rlist[i+1]->atoms[0]->chain) { 
      fprintf (stderr, "-%d", rlist[i]->id);
      fprintf (stderr, "\n");

      fprintf (stderr, "%c: ", rlist[i+1]->atoms[0]->chain);
      fprintf (stderr, "%d", rlist[i+1]->id);
      id = rlist[i+1]->id;
      n = 1;
      }

    else if (rlist[i]->id+1 != rlist[i+1]->id) {
      fprintf (stderr, "-%d", rlist[i]->id);
      n = 1;
      fprintf (stderr, ",%d", rlist[i+1]->id);
      }
 
    else {
      n += 1; 
      }
    }

  fprintf (stderr, "-%d", rlist[num_res-1]->id);
  fprintf (stderr, "\n");
  fprintf (stderr, "\n------ atoms ------\n");

  for (int i = 0; i < atoms.size(); i++) {
    fprintf (stderr, "id=%d chain=%c resid=%d name=%s\n", atoms[i].id, atoms[i].chain, 
             atoms[i].seq, atoms[i].name);
    }
  }

//*============================================================*
//*==========        printSecondaryStructure         ==========*
//*============================================================*
// print secondary structure.            

void
PmMolecule::printSecondaryStructure()
  {
  int n;
  n = helices.size() + sheets.size();
  fprintf (stderr, "\n---------- \"%s\" secondary structure ---------- \n", name.c_str());

  if (n == 0) {
    fprintf (stderr, "*** no secondary structure defined. \n");
    return;
    }

  fprintf (stderr, "number of helices = %d \n", helices.size());
  fprintf (stderr, "number of sheets  = %d \n", sheets.size());
  //fprintf (stderr, "number of loops   = %d \n", loops.size());

  fprintf (stderr, "\n------ helices:num id chain[sequence] ------\n");
  for (int i = 0; i < helices.size(); i++) {
    fprintf (stderr, "%d:  %s  %c[%d-%d] \n",  
             helices[i].num,
             helices[i].id.c_str(),
             helices[i].init_chain_id, 
             helices[i].init_seq_num, 
             helices[i].term_seq_num);
    }

  fprintf (stderr, "\n------ sheets:num id chain[sequence] ------\n");
  for (int i = 0; i < sheets.size(); i++) {
    fprintf (stderr, "%d:  %s  %c[%d-%d] \n", 
             sheets[i].num, 
             sheets[i].id.c_str(), 
             sheets[i].init_chain_id, 
             sheets[i].init_seq_num, 
             sheets[i].term_seq_num);
    }

  /*
  fprintf (stderr, "\n------ loops:chain[sequence] ------\n");
  for (int i = 0; i < loops.size(); i++) {
    fprintf (stderr, "%d: %c[%d-%d] \n", i+1, loops[i].init_chain_id, 
             loops[i].init_seq_num, loops[i].term_seq_num);
    }
  */
  }

//*============================================================*
//*==========              checkDescriptor           ==========*
//*============================================================*
// check if a descriptor returns any residues.

bool 
PmMolecule::checkDescriptor(const string desc)
  {
  vector<PmResidue*> rlist;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    return false;
    }

  return true;
  }

//*============================================================*
//*==========              getChain                  ==========*
//*============================================================*
// get a chain.

PmChain *
PmMolecule::getChain (char chain_id)
  {
  buildChains();
  PmChain *chain = NULL;

  for (int i = 0; i < num_chains; i++) {
    if (chains[i].id == chain_id) {
      chain = &chains[i];
      break;
      }
    }

  return (chain);
  }

//*============================================================*
//*==========              getChainAtoms             ==========*
//*============================================================*
// get the atoms for a chain.

void 
PmMolecule::getChainAtoms(char chain_id, vector<PmAtom>& catoms)
  {
  vector<PmAtom*> chain_atoms;
  catoms.clear(); 
  this->getChainAtoms(chain_id, "", chain_atoms);

  if (!chain_atoms.size()) {
    return;
    }

  for (unsigned int i = 0; i < chain_atoms.size(); i++) {
    catoms.push_back(*chain_atoms[i]);
    }
  }

//*============================================================*
//*==========              getChains                 ==========*
//*============================================================*
// get all the chains.

void 
PmMolecule::getChains (vector<PmChain*>& clist)
  {

  buildChains();

  if (num_chains == 0) {
    clist.clear();
    return;
    }

  clist.resize(num_chains);

  for (int i = 0; i < num_chains; i++) {
    clist[i] = &chains[i];
    }
  }

//*============================================================*
//*==========              getChainAtoms             ==========*
//*============================================================*
// get the atoms for a chain. use an atom name (e.g. CA).

void
PmMolecule::getChainAtoms (char chain_id, string name, vector<PmAtom*>& chain_atoms)
  {

  chain_atoms.clear();
  PmChain *chain = getChain (chain_id);

  if (!chain) {
    return;
    }

  bool no_name = name.empty();
  int num_res = chain->num_residues;
  int m = chain->res_index;

  for (int i = 0; i < num_res; i++) {
    PmResidue *res = &residues[m++];
    int num_atoms = res->num_atoms;

    for (int j = 0; j < num_atoms; j++) {
      PmAtom *atom = res->atoms[j];

      if (no_name) {
        chain_atoms.push_back (atom); 
        }
      else if (name == atom->name) {
        chain_atoms.push_back (atom); 
        }
      }
    }
  }

//*============================================================*
//*==========              getChainAtoms             ==========*
//*============================================================*
// get the atoms for a chain. use an atom filter.

void
PmMolecule::getChainAtoms (char chain_id, PmAtomFilter& filter, 
                           vector<PmAtom*>& chain_atoms)
  {

  PmChain *chain = getChain (chain_id);

  if (!chain) {
    return;
    }

  int num_res = chain->num_residues;
  int m = chain->res_index;

  for (int i = 0; i < num_res; i++) {
    PmResidue *res = &residues[m++];
    int num_atoms = res->num_atoms;

    for (int j = 0; j < num_atoms; j++) {
      PmAtom *atom = res->atoms[j];

      if (filter.hasName(atom->name)) {
        chain_atoms.push_back (atom);
        }
      }
    }
  }

//*============================================================*
//*==========              getChainAtomBonds         ==========*
//*============================================================*
// get the atom bonds for a chain.

void
PmMolecule::getChainAtomBonds (const char chain_id, vector<PmBond>& bonds)
  {
  bonds.clear();
  if (!num_chains) buildChains();
  PmChain *chain = getChain (chain_id);
  if (!chain) return;

  #define ndbg_getChainAtomBonds 
  #ifdef dbg_getChainAtomBonds 
  fprintf (stderr, "\n>>>>>> PmMolecule::getChainAtomBonds: \n");
  fprintf (stderr, "   >>> chain id [%c] \n", chain_id);
  fprintf (stderr, "   >>> num chains [%d] \n", num_chains);
  #endif

  // if we are processing an uknown compound then    //
  // check to see if there is connectivity for it.   //

  if (atoms[0].ctype == PM_COMPOUND_UNKNOWN) {
    #ifdef dbg_getChainAtomBonds 
    fprintf (stderr, "   >>> unknown compound \n");
    fprintf (stderr, "   >>> connect size %d\n", connect.size());
    #endif

    if (!connect.size()) {
      return;
      }

    int n, id1, id2;
    PmAtom *atom1, *atom2;
    int num_bonds = 0;

    for (unsigned int i = 0; i < connect.size();) {
      n = connect[i++];
      id1 = connect[i++];
      getAtom (id1, &atom1);

      if (!atom1) {
        fprintf (stderr, "   >>> atom id %d not found. \n", id1); 
        continue;
        }

      //fprintf (stderr, "   >>> %d:  %d ", n, id1); 

      for (int j = 0; j < n-1; j++) {
        id2 = connect[i++];
        getAtom (id2, &atom2);

        if (!atom2) {
          fprintf (stderr, "   >>> atom id %d not found. \n", id2); 
          continue;
          }

        addBond(atom1, atom2, num_bonds, bonds);
        //fprintf (stderr, " %d ", id2); 
        }

      //fprintf (stderr, "\n");
      }

    return;
    }

  // this is a known compound so assume    //
  // we know how it is bonded.             //

  vector<PmResidue*> rlist;
  PmResidue *res, *nres;
  getChainResidues (chain_id, rlist);
  int n = rlist.size();
  #ifdef dbg_getChainAtomBonds 
  fprintf (stderr, "   >>> number of residues %d \n", n);
  #endif

  for (int i = 0; i < n; i++) {
    res = rlist[i];
    //fprintf (stderr, "   >>> res[%d] id[%d]  mtype[%d]\n", i, res->id, res->mtype);

    if (i < n-1) {
      nres = rlist[i+1];
      }
    else {
      nres = NULL;
      }

    addResBond (res, nres, bonds);
    }

  // for a single residue add a bond to the next //
  // residue to provide continuity.              //

#ifdef createDomain_add_atom
  if (n == 1) {
    PmAtom *c, *n;
    int num = 0;
    res = rlist[0];
    res->getAtom ("C", &c);
    n = &atoms[num_atoms-1];

    if (c) {
      addBond (c, n, num, bonds);
      }
    }
#endif
  }

//*============================================================*
//*==========          getChainAtomColors            ==========*
//*============================================================*
// get the atom colors for a chain.

void
PmMolecule::getChainAtomColors (const char chain_id, PmVector3 **atom_colors)
  {

  if (!num_chains) buildChains();
  PmChain *chain = getChain (chain_id);
  if (!chain) return;

  vector<PmAtom*> chain_atoms;
  string atom_name;
  getChainAtoms (chain->id, atom_name, chain_atoms);
  int num = chain_atoms.size();
  PmVector3 *colors = new PmVector3 [num];

  for (int j = 0; j < num; j++) {
    PmAtom *atom = chain_atoms[j];
    atom->getColor (colors[j]);
    }

  *atom_colors = colors;
  }

//*============================================================*
//*==========              getChainAtomCoords        ==========*
//*============================================================*
// get the atom coordinates for a chain.

void
PmMolecule::getChainAtomCoords (char chain_id, string name, int *num_coords, 
                                PmVector3 **coords) 
  {

  //fprintf (stderr, "\n>>>>>> PmMolecule::getChainAtomCoords: \n");
  *num_coords = 0;
  *coords = NULL;
  PmChain *chain = getChain (chain_id);
  if (!chain) return;

  vector<PmAtom*> chain_atoms;
  getChainAtoms (chain_id, name, chain_atoms);
  int num_atoms = chain_atoms.size();
  if (!num_atoms) return;
  PmVector3 *vec = new PmVector3 [num_atoms];

  for (int i = 0; i < num_atoms; i++) {
    PmAtom *atom = chain_atoms[i];
    vec[i] = atom->pos;
    }

  *num_coords = num_atoms;
  *coords = vec;
  }

//*============================================================*
//*==========              getChainAtomRadii         ==========*
//*============================================================*
// get the atom radii for a chain.

void
PmMolecule::getChainAtomRadii (char chain_id, string name, int *num, float **rads)
  {
  *num = 0;
  *rads = NULL;
  PmChain *chain = getChain (chain_id);
  if (!chain) return;

  vector<PmAtom*> chain_atoms;
  getChainAtoms (chain_id, name, chain_atoms);
  int num_atoms = chain_atoms.size();
  if (!num_atoms) return;
  float *r = new float [num_atoms];

  for (int i = 0; i < num_atoms; i++) {
    PmAtom *atom = chain_atoms[i];
    r[i] = atom->getRadius();
    }

  *num = num_atoms;
  *rads = r;
  }

//*============================================================*
//*==========              getChainResidues          ==========*
//*============================================================*
// get the residues for a chain.

void
PmMolecule::getChainResidues (char chain_id, vector<PmResidue*>& rlist)
  {
  PmChain *chain = getChain (chain_id);
  if (!chain) return;

  int num_res = chain->num_residues;
  int n = chain->res_index;
  rlist.resize(num_res);

  for (int i = 0; i < num_res; i++) {
    rlist[i] = &residues[n++];
    }
  }

//*============================================================*
//*==========              getChainResidues          ==========*
//*============================================================*
// get the residues for a chain using a string.

void
PmMolecule::getChainResidues (char chain_id, char **desc, vector<PmResidue*>& rlist)
  {
  #define ndbg_getChainResidues 
  #ifdef dbg_getChainResidues 
  fprintf (stderr, "\n>>>>>> PmMolecule::getChainResidues:  chain [%c] \n", chain_id);
  #endif

  int i, j, n, m;
  int num, begin, end, rindex, pnum, size;
  char *str, val[2][80], name[20];
  PmCompoundType ctype;
  PmMoleculeType mtype;
  PmResidue *res, *list;

  PmChain *chain = getChain (chain_id);
  if (!chain) return;
  num = chain->num_residues;
  rindex = chain->res_index;
  list = &residues[rindex];

  #ifdef dbg_getChainResidues 
  fprintf (stderr, "   >>> desc [%s] \n", *desc);
  fprintf (stderr, "   >>> num res [%d] \n", num);
  #endif

  str = *desc;
  n = 0;
  m = 0;
  pnum = 0;
  size = 100;
  rlist.resize(size);
  #ifdef dbg_getChainResidues 
  fprintf (stderr, "   >>> parse: ");
  #endif

  while (*str != '\0') {
    char c = *str;
    char nc = *(str+1);
    #ifdef dbg_getChainResidues 
    fprintf (stderr, " %c ", c);
    #endif

    if (!isdigit(c) || (nc == '\0')) {
      if (nc == '\0') {
        val[m][n++] = c;
        }

      val[m][n] = '\0';

      if (c == '-') {
        n = 0;
        m++;
        }

      else if ((c == ',') || (c == ']') || (nc == '\0')) {
        begin = atoi (val[0]); 

        if (!m) {
          end = begin; 
          }
        else {
          end = atoi (val[1]); 
          }

        #ifdef dbg_getChainResidues 
        fprintf (stderr, "\n");
        fprintf (stderr, "begin [%d]   end [%d] \n", begin, end);
        fprintf (stderr, "   >>> search residues: ");
        #endif

        for (j = 0, i = begin; i <= end && j < num; j++) {
          res = &list[j];
          #ifdef dbg_getChainResidues 
          fprintf (stderr, " id[%d]", res->id);
          #endif

          if ((res->id >= begin) && (res->id <= end)) {
            #ifdef dbg_getChainResidues 
            fprintf (stderr, ":found");
            #endif

            if (pnum == size) {
              size += 100;
              rlist.resize(size);
              }

            rlist[pnum++] = res;
            }
          }

        m = 0;
        n = 0;

        if (c == ']') {
          break;
          }
        }

 
      // parse a compund name //

      else { 
        n = 0;

        while (1) { 
          c = *str;

          if ((c == ',') || (c == ']')) {
            name[n++] = '\0';
            pm_CmpdConvType(name, ctype, mtype);

            for (i = 0; i < num; i++) {
              res = &list[i];

              if (res->ctype == ctype) { 
                if (pnum == size) {
                  size += 100;
                  rlist.resize(size);
                  }

                rlist[pnum++] = res;
                }
              }

            break;
            }
          
          name[n++] = c;
          str++;
          }

        if (c == ']') {
          break;
          }

        str++;
        }
      }
    else {
      val[m][n++] = c;
      }
   
    str++;
    }

  str++;
  *desc = str;
  rlist.resize(pnum);
  #ifdef dbg_getChainResidues 
  fprintf (stderr, "   >>> pnum [%d] \n", pnum);
  #endif
  }

//*============================================================*
//*==========          getChainPlaneGeometry         ==========*
//*============================================================*
// get the plane geometry for a chain.

void
PmMolecule::getChainPlaneGeometry(char chain_id, PmResidue *aux_res,
                                  int& num_planes, PmConn **p_conn, int& num_coords, 
                                  PmVector3 **p_coords)
  {
  PmChain *chain = getChain (chain_id);

  if (!chain) {
    return;
    }

  //fprintf (stderr, "\n>>>>>> PmMolecule_getChainPlaneGeometry \n");
  vector<PmMolPlane> planes;
  vector<PmResidue*> rlist;
  PmVector3 *coords;
  PmResidue *res, *nres;
  int num_verts;
  getChainResidues (chain_id, rlist);
  int n = rlist.size();

  for (int i = 0; i < n; i++) {
    res = rlist[i];

    if (i == n-1) {
      nres = aux_res;
      }
    else {
      nres = rlist[i+1];
      }

    addResPlane (res, nres, planes);
    }

  num_coords = 0;
  *p_coords = NULL;
  num_planes = planes.size();
  //fprintf (stderr, "   >>> num planes = %d \n", num_planes);

  for (int i = 0; i < num_planes; i++) {
    num_coords += planes[i].vertices.size();
    }

  //fprintf (stderr, "   >>> num coords = %d \n", num_coords); 

  if (num_coords == 0) {
    return;
    }

  PmConn *conn = new PmConn(num_planes + num_coords);
  n = 0;
  coords = new PmVector3[num_coords];
  num_coords = 0;

  for (unsigned int i = 0; i < planes.size(); i++) {
    num_verts = planes[i].vertices.size(); 
    (*conn)[n++] = num_verts;

    for (int j = 0; j < num_verts; j++) {
      coords[num_coords] = planes[i].vertices[j]; 
      (*conn)[n++] = num_coords;
      num_coords += 1;
      }
    }

  *p_conn = conn;
  *p_coords = coords;
  }

//*============================================================*
//*==========               checkClash               ==========*
//*============================================================*
// check if two domains have atoms that clash.

void 
PmMolecule::checkClash(PmMolecule* dom, bool show)
  {
  }

//*============================================================*
//*==========               addConnect               ==========*
//*============================================================*
// add atom connectivity.

void 
PmMolecule::addConnect(const int n, int conn[])
  {
  int id;
  connect.push_back(n);

  for (int i = 0; i < n; i++) {
    id = conn[i];
    connect.push_back(id);
    }

  /*
  for (int i = 0; i < n; i++) {
    id = conn[i];
    id_found = false;
    fprintf (stderr, "\n id %d: ", id);

    for (int j = 0; j < num_atoms; j++) {
      fprintf (stderr, "%d ", atoms[j].id);

      if (id == atoms[j].id) {
        connect.push_back(j+1);
        id_found = true;
        break;
        }
      }

    if (!id_found) {
      pm_ErrorReport (PM, "connect atom id %d not found in \"%s\".", "*", id, 
                      name.c_str());
      }
    }
  */
  }

//*============================================================*
//*==========              getDimensions             ==========*
//*============================================================*
// get molecule dimensions.

void
PmMolecule::getDimensions(vector<float>& dims)
  {
  if (dimensions.size() == 0) {
    PmPcaResults pca;
    PmAtomFilter filter;
    this->compPrincipalComponents(filter, pca);

    dimensions.push_back(pca.s1);
    dimensions.push_back(pca.s2);
    dimensions.push_back(pca.s3);
    }

  dims = dimensions;
  }

//*============================================================*
//*==========              getRadii                  ==========*
//*============================================================*
// get the radii for the atoms.

void 
PmMolecule::getRadii(vector<float>& rads)
  {
  rads.clear();
  vector<PmResidue*> rlist;
  string desc;
  PmResidue *res;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    return;
    }

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];
      rads.push_back(atom->getRadius());
      }
    }
  }

//*============================================================*
//*==========              getCoordinates            ==========*
//*============================================================*
// get the coordinates for the atoms.

void
PmMolecule::getCoordinates(vector<PmVector3>& coords)
  {
  coords.clear();
  vector<PmResidue*> rlist;
  string desc;
  PmResidue *res;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    return;
    }

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];
      coords.push_back(atom->pos);
      }
    }
  }

//*============================================================*
//*==========            checkDomainContact          ==========*
//*============================================================*
// check contact between domains. 

void 
PmMolecule::checkDomainContact(vector<PmMolecule*>& dom_list, const float tol,
                               PmVector3& color, const bool show)
  {
  string desc, dname;
  PmResidue *res, *tres;
  PmAtomFilter filter;
  PmAtom *atom, *tatom;
  vector<PmAtom*> res_atoms;
  vector<PmVector3>center, tcenter; 
  PmVector3 v;
  vector<float> radius, tradius; 
  int num_res, num_tres;
  vector<PmResidue*> rlist, trlist;
  int num_contact;
  float d, r, tr; 
  PmMolecule *dom;
  vector<PmVector3> vertices; 

  // get the residues for this molecule //

  this->getResidues(desc, rlist);
  this->getResidueBounds(center, radius);
  num_res = center.size();

  fprintf (stderr, "\n------ check for domain-domain contact ------ \n");
  fprintf (stderr, ">>> tolerance = %f \n", tol);
  fprintf (stderr, ">>> source domain %s: number of residues = %d \n", 
           this->name.c_str(), num_res);

  for (unsigned int di = 0; di < dom_list.size(); di++) {
    dom = dom_list[di];
    tcenter.clear();
    tradius.clear();
    dom->getName(dname);
    dom->getResidueBounds(tcenter, tradius);
    num_tres = tcenter.size();
    fprintf (stderr, ">>> checking %s: number of residues = %d  ", dname.c_str(),
             num_tres);
    dom->getResidues(desc, trlist);
    num_contact = 0;
    vertices.clear();

    for (int i = 0; i < num_res; i++) {
      res = rlist[i];

      for (int j = 0; j < num_tres; j++) {
        tres = trlist[j];
        v = center[i] - tcenter[j];
        d = v.length();

        if (d <= (radius[i] + tradius[j])) {
          for (int ia = 0; ia < res->num_atoms; ia++) {
            atom = res->atoms[ia];
            r = atom->getRadius();

            for (int ib = 0; ib < tres->num_atoms; ib++) {
              tatom = tres->atoms[ib];
              v = atom->pos - tatom->pos;
              tr = tatom->getRadius();
              d = v.length();

              if (d <= (r + tr - tol)) {
                vertices.push_back(atom->pos);
                vertices.push_back(tatom->pos);
                num_contact += 1;
                //fprintf (stderr, "%d: d[%f]   r[%f] tr[%f] \n", num_contact, d, r, tr); 
                //fprintf (stderr, "    res[%d] tres[%d] \n", atom->seq, tatom->seq); 
                }
              }
            }
          }
        }
      }

    fprintf (stderr, "  number of contacts = %d  \n", num_contact); 

    if (vertices.size()) {
      PmGraphicsLine *line;
      PmGraphicsAttributes atts;
      PmGraphicsGeometry *geom;
      int num_verts = vertices.size();
      PmVector3 *gverts = new PmVector3[num_verts];
      string geom_name = "contact";
      //geom_name = geom_name + '[' + this->name + ':' + dname + ']' + "lines";
      PmVector3 pt1, pt2; 

      for (int i = 0; i < num_verts; i++) {
        gverts[i] = vertices[i];
        }

      /*
      for (int i = 0; i < num_verts/2; i++) {
        pt1 = gverts[2*i];
        pt2 = gverts[2*i+1];
        v = pt2 - pt1;
        d = v.length();
        fprintf (stderr, "%d  %f  \n", i, d); 
        }
      */

      getGraphicsGeometry(geom_name, &geom);

      if (!geom) {
        line = new PmGraphicsLine(geom_name, num_verts, gverts);
        addGraphicsGeometry (line);
        }
      else {
        line = dynamic_cast<PmGraphicsLine*>(geom);
        }

      atts.setColor(color);
      atts.disjoint = true;
      line->setAttributes(atts);
      line->display();
      }
    }
  }

//*============================================================*
//*==========                displayResContact       ==========*
//*============================================================*
// show residue-residue contact.

void 
PmMolecule::displayResContact(const float cutoff, const bool use_radii, 
                              PmAtomFilter& filter, PmVector3& color, bool show)
  {
  int num_res, n;
  vector<PmResidue*> rlist;
  PmResidue *res1, *res2;
  PmVector3 v, center1, center2, *center;
  float d, *radius;
  int num_contact;
  PmAtom *atom, *atom1, *atom2;
  PmAtomFilter sc_filter;
  PmExtent extent;
  float dx, dy, dz, r1, r2;

  vector<PmVector3> vertices;
  string geom_name, atom_name;
  PmGraphicsLine *line;
  PmGraphicsGeometry *geom;
  int num_verts;
  PmVector3 *gverts;
  PmGraphicsAttributes atts;
 
  // get the residues for this molecule //

  this->getResidues("", rlist);
  num_res = rlist.size();
  num_contact = 0;

  // use atomic radii plus a cutoff //

  if (use_radii) {
    center = new PmVector3[num_res];
    radius = new float[num_res];
    this->getBackboneAtomNames(filter.names);
    sc_filter.exclude = true;

    // first compute the center's and extent for each residue //

    for (int i = 0; i < num_res; i++) {
      res1 = rlist[i];
      center[i].set(0,0,0);
      n = 0;
      extent.min[0] =  1e6; extent.min[1] =  1e6; extent.min[2] = 1e6;
      extent.max[0] = -1e6; extent.max[1] = -1e6; extent.max[2] = -1e6;

      for (int ia = 0; ia < res1->num_atoms; ia++) {
        atom = res1->atoms[ia];

        if (filter.hasName(atom->name)) {
          center[i] = center[i] + atom->pos;
          extent.update(atom->pos);
          n += 1;
          }
        }

      center[i] = center[i] / n;
      dx = extent.max[0] - extent.min[0];
      dy = extent.max[1] - extent.min[1];
      dz = extent.max[2] - extent.min[2];

      if (dx > dy) {
        radius[i] = dx;
        }
      else {
        radius[i] = dy;
        }

      if (dz > radius[i]) {
        radius[i] = dz;
        }

      radius[i] = radius[i] / 2.0 + 0.2;
      }

    num_contact = 0;

    for (int i = 0; i < num_res; i++) {
      res1 = rlist[i];

      for (int j = i+1; j < num_res; j++) {
        res2 = rlist[j];

        if (i != j) {
          v = center[i] - center[j];
          d = v.length();

          if (d <= (radius[i] + radius[j])) {
            for (int ia = 0; ia < res1->num_atoms; ia++) {
              atom1 = res1->atoms[ia];
              r1 = atom1->getRadius();

              if (sc_filter.hasName(atom1->name)) {
                for (int ib = 0; ib < res2->num_atoms; ib++) {
                  atom2 = res2->atoms[ib];

                  if (sc_filter.hasName(atom2->name)) {
                    v = atom1->pos - atom2->pos;
                    r2 = atom2->getRadius();
                    d = v.length();

                    if (d <= (r1 + r2 + cutoff)) {
                      num_contact += 1;
                      vertices.push_back(atom1->pos);
                      vertices.push_back(atom2->pos);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


  // use residues //

  else {

    for (int i = 0; i < num_res; i++) {
      res1 = rlist[i];
      center1.set(0,0,0);
      n = 0;

      for (int ia = 0; ia < res1->num_atoms; ia++) {
        atom = res1->atoms[ia];

        if (filter.hasName(atom->name)) {
          center1 = center1 + atom->pos;
          n += 1;
          }
        }

      if (!n) {
        continue;
        }

      center1 = center1 / n;

      for (int j = i+1; j < num_res; j++) {
        res2 = rlist[j];
        center2.set(0,0,0);
        n = 0;

        for (int ia = 0; ia < res2->num_atoms; ia++) {
          atom = res2->atoms[ia];

          if (filter.hasName(atom->name)) {
            center2 = center2 + atom->pos;
            n += 1;
            }
          }
    
        if (n) {
          center2 = center2 / n;
          v = center1 - center2;
          d = v.length();

          if (d <= cutoff) { 
            num_contact += 1;
            vertices.push_back(center1);
            vertices.push_back(center2);
            }
          }
        }
      }
    }

  if (!vertices.size()) {
    return;
    }

  num_verts = vertices.size();
  gverts = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gverts[i] = vertices[i];
    }

  buildGeomName ('\0', "residueContact", "", geom_name);
  getGraphicsGeometry(geom_name, &geom);

  if (!geom) {
    line = new PmGraphicsLine(geom_name, num_verts, gverts);
    addGraphicsGeometry (line);
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(geom);
    }

  fprintf (stderr, "    >>> number of contacts = %d \n", num_contact);
  atts.setColor(color);
  atts.disjoint = true;
  line->setAttributes(atts);
  line->display();
  }

//*============================================================*
//*==========                getResContact           ==========*
//*============================================================*
// compute residue-residue contact.

void
PmMolecule::getResContact()
  {

  if (this->residue_contacts) { 
    return;
    }

  //fprintf (stderr, "\n>>>>>> PmMolecule_getResContact \n");
  string desc;
  PmResidue *res;
  PmAtomFilter filter;
  PmAtom *atom;
  vector<PmAtom*> res_atoms;
  PmVector3 *center, v;
  float *radius, d;
  int num_res;
  vector<PmResidue*> rlist;
  int n, num_contact;
  PmExtent extent;
  float dx, dy, dz;
  PmResidueContact *res_con;
  PmResidue **rclist, *list[100]; 

  // get the residues for this molecule //

  this->getResidues(desc, rlist);

  // just use sidechain atoms //

  this->getBackboneAtomNames(filter.names);
  filter.exclude = true;
  num_res = rlist.size();
  center = new PmVector3[num_res];
  radius = new float[num_res];

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];
    center[i].set(0,0,0);
    n = 0;
    extent.min[0] =  1e6; extent.min[1] =  1e6; extent.min[2] = 1e6;
    extent.max[0] = -1e6; extent.max[1] = -1e6; extent.max[2] = -1e6;

    for (int ia = 0; ia < res->num_atoms; ia++) {
      atom = res->atoms[ia];

      if (filter.hasName(atom->name)) {
        res_atoms.push_back(atom);
        center[i] = center[i] + atom->pos;
        extent.update(atom->pos);
        n += 1;
        }
      }
 
    center[i] = center[i] / n;
    dx = extent.max[0] - extent.min[0];
    dy = extent.max[1] - extent.min[1];
    dz = extent.max[2] - extent.min[2];

    if (dx > dy) {
      radius[i] = dx;
      }
    else {
      radius[i] = dy;
      }

    if (dz > radius[i]) {
      radius[i] = dz;
      }

    radius[i] = radius[i] / 2.0 + 0.2;
    }

  num_contact = 0;
  res_con = new PmResidueContact[num_res]; 

  for (int i = 0; i < num_res; i++) {
    res_con[i].num = 0; 
    res_con[i].list = NULL; 
    }

  for (int i = 0; i < num_res; i++) {
    n = 0;

    for (int j = 0; j < num_res; j++) {
      if (i != j) {
        v = center[i] - center[j];
        d = v.length();

        if (d < (radius[i] + radius[j])) {
          num_contact += 1;
          list[n] = rlist[j];
          n += 1;
          }
        }
      }

    if (n) {
      rclist = new PmResidue*[n]; 

      for (int j = 0; j < n; j++) {
        rclist[j] = list[j]; 
        }

      res_con[i].num = n;
      res_con[i].list = rclist;
      }
    }

  delete[] center;
  delete[] radius;
  this->residue_contacts = res_con;

#define nprint_PmMolecule_getResContact
#ifdef print_PmMolecule_getResContact

  char ic, jc;

  fprintf (stderr, "\n\n--------- residue contacts = %d -------- \n", num_contact);

  for (int i = 0; i < num_res; i++) {
    ic = rlist[i]->getChainId();
    fprintf (stderr, "%d [%c%d]: ", i+1, ic, rlist[i]->id); 

    for (int j = 0; j < res_con[i].num; j++) {
      jc = res_con[i].list[j]->getChainId();
      fprintf (stderr, "%c%d ", jc, res_con[i].list[j]->id); 
      }

    fprintf (stderr, "\n");
    }
#endif
  }

//*============================================================*
//*==========                joinDomains             ==========*
//*============================================================*
// join the src domain onto this domain.

PmMolecule *
PmMolecule::joinDomains(const string name, const string desc, PmMolecule *src_dom, 
                        const string src_desc, const bool nterm,
                        const vector<string>& chain_ids)
  {
  vector<PmResidue*> res, src_res;
  int num_atoms;
  PmAtom **atoms, *atom;
  string atom_names, chain_sub;
  vector<PmAtom*> chain_atoms; 
  char chain_id, src_chain_id, sub_chain_id; 

  #define dbg_joinDomains
  #ifdef dbg_joinDomains
  fprintf (stderr, "\n>>>>>> PmMolecule::joinDomain:  desc[%s] src_desc[%s] \n",
           desc.c_str(), src_desc.c_str());
  fprintf (stderr, ">>> nterm = %d \n", nterm);
  #endif

  vector<PmAtom> dst_atoms;
  src_dom->getAtoms(dst_atoms);

  this->getResidues(desc, res);
  chain_id = res[0]->getChainId();
  src_dom->getResidues(src_desc, src_res);
  src_chain_id = src_res[0]->getChainId();

  // create domain //

  PmMolecule *domain = new PmMolecule(name);

  // first add atoms that are not from the joining chains //

  for (int i = 0; i < num_chains; i++) {
    if (chains[i].id != chain_id) {
      getChainAtoms(chains[i].id, atom_names, chain_atoms);

      for (unsigned int j = 0; j < chain_atoms.size(); j++) { 
        atom = chain_atoms[j];
        domain->addAtom (atom);
        }
      }
    }

  for (int i = 0; i < src_dom->num_chains; i++) {
    if (src_dom->chains[i].id != src_chain_id) {
      src_dom->getChainAtoms(src_dom->chains[i].id, atom_names, chain_atoms);
      sub_chain_id = '\0';

      // check for renaming of chain ids //

      for (unsigned int j = 0; j < chain_ids.size(); j++) {
        chain_sub = chain_ids[j];
   
        if (chain_sub[0] == src_dom->chains[i].id) {
          sub_chain_id = chain_sub[2];
          break;
          }
        }

      if (sub_chain_id != '\0') {
        for (unsigned int j = 0; j < chain_atoms.size(); j++) { 
          atom = chain_atoms[j];
          domain->addAtom(atom->id, atom->name, atom->type, atom->pos, sub_chain_id, 
                          atom->ctype, atom->mtype, atom->seq); 
          }
        }
      else {
        for (unsigned int j = 0; j < chain_atoms.size(); j++) { 
          atom = chain_atoms[j];
          domain->addAtom (atom);
          }
        }
      }
    }

  // now add atoms from both domains that are to be joined //

  // add to nterm //

  if (nterm) {
    int start_res, res_seq, atom_id, end_res;

    for (unsigned int i = 0; i < src_dom->residues.size(); i++) {
      if ((src_dom->residues[i].atoms[0]->chain == src_chain_id) &&
          (src_dom->residues[i].id == src_res[0]->id)) {
        end_res = i;
        break;
        }
      }

    for (unsigned int i = 0; i < residues.size(); i++) {
      if ((residues[i].atoms[0]->chain == chain_id) &&
          (residues[i].id == res[0]->id)) { 
        start_res = i;
        break;
        }
      }

    res_seq = residues[start_res].id - end_res;

    if (res_seq < 0) {
      pm_ErrorReport (PM, "can't renumber nterm sequence.", "*");
      return NULL;
      }

    atom_id = 1;
    #ifdef dbg_joinDomains
    fprintf (stderr, "  res size = %d \n", src_dom->residues.size());
    fprintf (stderr, "  end res #  = %d \n", end_res);
    fprintf (stderr, "  end res id = %d \n", src_dom->residues[end_res].id);
    fprintf (stderr, "  res seq   = %d \n", res_seq);
    fprintf (stderr, "  start res id = %d \n", residues[start_res].id);
    #endif 

    // add src domain's atoms //

    for (unsigned int i = 0; i < src_dom->residues.size(); i++) {
      if (src_dom->residues[i].id == src_res[0]->id) {
        break;
        }

      num_atoms = src_dom->residues[i].num_atoms;
      atoms = src_dom->residues[i].atoms;

      if (atoms[0]->chain != chain_id) {
        break;
        }

      for (int j = 0; j < num_atoms; j++) {
        atom = atoms[j];
        domain->addAtom(atom_id, atom->name, atom->type, atom->pos, chain_id,
                        atom->ctype, atom->mtype, res_seq);
        atom_id += 1;
        }

      res_seq += 1;
      }

    // add this domain's atoms //

    for (unsigned int i = start_res; i < residues.size(); i++) {
      //fprintf (stderr, "  res[%d] \n", residues[i].id);
      num_atoms = residues[i].num_atoms;
      atoms = residues[i].atoms;

      for (int j = 0; j < num_atoms; j++) {
        domain->addAtom (atoms[j]);
        atom_id = atoms[j]->id;
        //fprintf (stderr, "  atom[%d] \n", atoms[j]->id);
        }
      }
    }

  // add to cterm //

  else {
    int start_res, res_seq, atom_id;

    for (unsigned int i = 0; i < residues.size(); i++) {
      if (residues[i].atoms[0]->chain == chain_id) {
        start_res = i;
        break;
        }
      }

    // add this domain's atoms //

    for (unsigned int i = start_res; i < residues.size(); i++) {
      if (residues[i].id == res[0]->id) { 
        break;
        }

      //fprintf (stderr, "  res[%d] \n", residues[i].id);
      num_atoms = residues[i].num_atoms;
      atoms = residues[i].atoms;

      for (int j = 0; j < num_atoms; j++) {
        domain->addAtom (atoms[j]);
        atom_id = atoms[j]->id;
        //fprintf (stderr, "  atom[%d] \n", atoms[j]->id);
        }
      }

    // add source domain's atoms //

    #ifdef dbg_joinDomains
    fprintf (stderr, ">>> add src atoms. num res = %d \n", src_dom->residues.size());
    fprintf (stderr, "    src chain = %c \n", src_chain_id);
    fprintf (stderr, "    src res id = %d \n", src_res[0]->id);
    #endif
    res_seq = res[0]->id;
    atom_id += 1;

    for (unsigned int i = 0; i < src_dom->residues.size(); i++) {
      if ((src_dom->residues[i].atoms[0]->chain == src_chain_id) &&
          (src_dom->residues[i].id == src_res[0]->id)) {
        start_res = i;
        break;
        }
      }

    #ifdef dbg_joinDomains
    fprintf (stderr, ">>> start res id = %d \n", src_dom->residues[start_res].id);
    #endif

    for (unsigned int i = start_res; i < src_dom->residues.size(); i++) {
      //fprintf (stderr, "  res[%d] \n", src_dom->residues[i].id);
      num_atoms = src_dom->residues[i].num_atoms;
      atoms = src_dom->residues[i].atoms;

      if (atoms[0]->chain != src_chain_id) { 
        break;
        }

      for (int j = 0; j < num_atoms; j++) {
        atom = atoms[j];
        domain->addAtom(atom_id, atom->name, atom->type, atom->pos, chain_id, 
                        atom->ctype, atom->mtype, res_seq); 
        atom_id += 1;
        }

      res_seq += 1;
      }
    }

  domain->type = type;
  domain->backbone_atom = backbone_atom;
  return domain;
  }

//*============================================================*
//*==========                Helix                   ==========*
//*============================================================*
// add a helix.

void 
PmMolecule::addHelix(PmHelix& helix) 
  {
  helices.push_back(helix);
  }

void
PmMolecule::getHelices(vector<PmHelix>& helices)
  {
  helices = this->helices;
  }

//*============================================================*
//*==========                getHelix                ==========*
//*============================================================*
// get a helix.

bool
PmMolecule::getHelix(const string id, PmHelix& helix)
  {

  // first search for helix by id //

  for (unsigned int i = 0; i < helices.size(); i++) {
    if (id == helices[i].id) { 
      helix = helices[i];
      return true;
      }
    }

  if (id.size() == 1) {
    return false;
    }

  // try searching for helix using id as a residue //

  char chain_id;
  int res;
  istringstream inputString(id);
  inputString >> chain_id >> res; 
  //fprintf (stderr, "chain id[%c]  res[%d] \n", chain_id, res);

  for (unsigned int i = 0; i < helices.size(); i++) {
    if ((chain_id == helices[i].init_chain_id) && (res >= helices[i].init_seq_num) && 
        (res <= helices[i].term_seq_num)) {
      helix = helices[i];
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========            displayHelixProps           ==========*
//*============================================================*
// display helix properties.

void
PmMolecule::displayHelixProps(vector<PmHelixProps>& helix_props, 
                              PmGraphicsAttributes& atts)
  {
  fprintf (stderr, ">>>> PmMolecule::displayHelixProps \n");
  fprintf (stderr, ">>> size = %d \n", helix_props.size());
  char chain_id = '\0';
  string geom_name, pgeom_name;
  PmGraphicsLine *line;
  PmGraphicsPoint *points;
  PmGraphicsGeometry *geom, *pgeom;
  int num;
  PmVector3 *verts;

  for (int i = 0; i < helix_props.size(); i++) {
    num = helix_props[i].origin.size();

    if (num == 0) {
      return;
      }

    chain_id = helix_props[i].chain_id;
    buildGeomName (chain_id, "helixPropsLine", "", geom_name);
    getGraphicsGeometry(geom_name, &geom);

    buildGeomName (chain_id, "helixPropsPoint", "", pgeom_name);
    getGraphicsGeometry(pgeom_name, &pgeom);

    if (!geom) {
      verts = new PmVector3[num];

      for (int j = 0; j < num; j++) {
        verts[j] = helix_props[i].origin[j];
        }

      line = new PmGraphicsLine(geom_name, num, verts);
      addGraphicsGeometry (line);

      points = new PmGraphicsPoint(pgeom_name, num, verts);
      addGraphicsGeometry (points);
      }
    else {
      line = dynamic_cast<PmGraphicsLine*>(geom);
      points = dynamic_cast<PmGraphicsPoint*>(pgeom);
      }

    line->setAttributes(atts);
    line->display();

    points->setAttributes(atts);
    points->display();
    }
  }

//*============================================================*
//*==========            getHelixProps               ==========*
//*============================================================*
// get helix properties for a residue.

void
PmMolecule::getHelixProps(const string desc, PmHelixProps& props)
  {
  //fprintf (stderr, "\n>>>>>> PmMolecule::getHelixProps \n");
  //fprintf (stderr, ">>> desc = %s \n", desc.c_str());
  vector<PmHelixProps> helix_props;
  vector<PmResidue*> rlist;
  PmResidue *res;
  PmVector3 origin;

  // get residues from desc  //

  getResidues(desc, rlist);

  if (rlist.size() == 0) {
    return;
    }

  this->getHelixProps(helix_props);

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (unsigned int j = 0; j < helix_props[0].res_id.size(); j++) {
      //fprintf (stderr, ">>> res = %d \n", helix_props[0].res_id[j]); 

      if (res->id == helix_props[0].res_id[j]) {
        props.chain_id = helix_props[0].chain_id;
        props.res_id.push_back(helix_props[0].res_id[j]);
        props.origin.push_back(helix_props[0].origin[j]);
        props.axis.push_back(helix_props[0].axis[j]);
        props.radius.push_back(helix_props[0].radius[j]);
        return;
        }
      }
    }
  }

//*============================================================*
//*==========            getHelixProps               ==========*
//*============================================================*
// get helix properties. 

void 
PmMolecule::getHelixProps(vector<PmHelixProps>& helix_props)
  {
  //fprintf (stderr, ">>>> PmMolecule::getHelixProps \n");
  buildChains();
  PmHelixProps props;
  PmVector3 ca_pos[4], v10, v21, v32, dv1, dv2, axis, origin;
  PmAtom *atom;
  float dmag, emag, ctheta, twist, res_per_turn, radius, height;
  //fprintf (stderr, ">>>> num chains [%d] \n", num_chains);

  for (int ci = 0; ci < num_chains; ci++) {
    PmChain chain = chains[ci];
    char chain_id = chains[ci].id;
    //fprintf (stderr, "\n------ chain [%c] ------ \n", chain_id);

    int num_res = chain.num_residues;
    int m = chain.res_index;
    helix_props.push_back(props);
    helix_props[ci].chain_id = chain_id; 
    //fprintf (stderr, ">>> num res[%d] \n", num_res);

    for (int i = 0; i < num_res-3; i++) {
      PmResidue *res = &residues[m];

      for (int j = 0; j < 4; j++) {
        residues[m+j].getAtom("CA", &atom);
        ca_pos[j] = atom->pos;
        }

      // compute vectors joining ca atoms // 
     
      v10 = ca_pos[1] - ca_pos[0];
      v21 = ca_pos[2] - ca_pos[1];
      v32 = ca_pos[3] - ca_pos[2];
      dv1 = v10 - v21;
      dv2 = v21 - v32;

      // helix axis  //
      axis = dv1.cross(dv2);
      axis.normalize();

      dmag = dv1.length();
      emag = dv2.length();
      ctheta = dv1*dv2 / (dmag*emag);

      twist = (180.0/M_PI) * acos(ctheta);
      res_per_turn = 360.0 / twist;
      radius = sqrt(dmag*emag) / (2.0*(1.0 - ctheta));
      height = v21 * axis;
      origin = ca_pos[1] - (radius/dmag)*dv1;

      helix_props[ci].res_id.push_back(res->id); 
      helix_props[ci].origin.push_back(origin); 
      helix_props[ci].axis.push_back(axis); 
      helix_props[ci].ca_pos.push_back(ca_pos[0]); 
      helix_props[ci].twist_angle.push_back(twist); 
      helix_props[ci].res_per_turn.push_back(res_per_turn); 
      helix_props[ci].radius.push_back(radius); 
      m += 1;
      }
    }
  }

//*============================================================*
//*==========              getHelixResidues          ==========*
//*============================================================*
// get the resides for a molecule from a helix descriptor.

void
PmMolecule::getHelixResidues(const string dstr, vector<PmResidue*>& rlist)
  {
  //fprintf (stderr, "\n>>>>>> PmMolecule::getHelixResidues desc [%s] \n", dstr.c_str());
  char str[20], *desc;
  PmHelix helix;
  string helix_id, res_desc;

  // parse and process helix descriptor. the descriptor  //
  // will be of the form helix:<id>                      //

  desc = (char*)dstr.c_str();
  for (int i = 0; i < 5; i++) desc++;

  while (desc && (*desc != '\n') && (*desc != ' ') && (*desc != '\0')) {
    helix_id.push_back(*desc);
    /*
    n += 1;

    if (n > 3) {
      pm_ErrorReport (PM, "helix id too long (> 3 characters).", "*");
      return;
      }
    */
    desc++;
    }

  if (!getHelix(helix_id, helix)) {
    pm_ErrorReport (PM, "helix \"%s\" not found.", "*", helix_id.c_str());
    return;
    }

  sprintf (str, "%c[%d-%d]", helix.init_chain_id, helix.init_seq_num, helix.term_seq_num);
  res_desc = str;
  //fprintf (stderr, "    >>>> helix residues %s \n", str);
  getResidues(res_desc, rlist);
  }

//*============================================================*
//*==========                Loop                    ==========*
//*============================================================*
// add a loop. 

void
PmMolecule::addLoop(PmLoop& loop)
  {
  loops.push_back(loop);
  }

void
PmMolecule::getLoops(vector<PmLoop>& loops)
  {
  loops = this->loops;
  }

//*============================================================*
//*==========              getLoopResidues           ==========*
//*============================================================*
// get the resides for a molecule from a loop descriptor.

void
PmMolecule::getLoopResidues(const string dstr, vector<PmResidue*>& rlist)
  {
  //fprintf (stderr, "\n>>>>>> PmMolecule::getLoopResidues desc [%s] \n", dstr.c_str());
  char str[20], *desc;
  PmLoop loop;
  string loop_id, res_desc, edesc;
  vector<PmResidue*> crlist;
  int init_seq_num, term_seq_num;

  // parse and process loop descriptor. the descriptor  //
  // will be of the form loop:<id>                      //

  desc = (char*)dstr.c_str();
  for (int i = 0; i < 4; i++) desc++;
  int n = 0;

  while (desc && (*desc != '\n') && (*desc != ' ') && (*desc != '\0')) {
    loop_id.push_back(*desc);
    desc++;
    }

  if (!getLoop(loop_id, loop)) {
    pm_ErrorReport (PM, "loop \"%s\" not found.", "*", loop_id.c_str());
    return;
    }

  getChainResidues (loop.init_chain_id, crlist);
  n = crlist.size();
  init_seq_num = loop.init_seq_num;
  term_seq_num = loop.term_seq_num;

  /*
  if (loop.init_seq_num > crlist[0]->id) {
    init_seq_num -= 1;
    }

  if (loop.term_seq_num < crlist[n-1]->id) {
    term_seq_num += 1;
    }
  */

  sprintf (str, "%c[%d-%d]", loop.init_chain_id, init_seq_num, term_seq_num);
  //sprintf (str, "%c[%d-%d]", loop.init_chain_id, loop.init_seq_num, loop.term_seq_num);
  res_desc = str;
  //fprintf (stderr, "    >>>> loop residues %s \n", str);
  getResidues(res_desc, rlist);
  }

//*============================================================*
//*==========                getLoop                 ==========*
//*============================================================*
// get a loop. 

bool
PmMolecule::getLoop(const string id, PmLoop& loop)
  {

  // first search for loop by id //

  for (unsigned int i = 0; i < loops.size(); i++) {
    if (id == loops[i].id) {
      loop = loops[i];
      return true;
      }
    }

  if (id.size() == 1) {
    return false;
    }

  // try searching for helix using id as a residue //

  char chain_id;
  int res;
  istringstream inputString(id);
  inputString >> chain_id >> res;
  //fprintf (stderr, "chain id[%c]  res[%d] \n", chain_id, res);

  for (unsigned int i = 0; i < helices.size(); i++) {
    if ((chain_id == loops[i].init_chain_id) && (res >= loops[i].init_seq_num) &&
        (res <= loops[i].term_seq_num)) {
      loop = loops[i];
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========              setMode                   ==========*
//*============================================================*
// set mode display attributes/

void 
PmMolecule::setModeColor(const PmVector3& color){
  this->mode_color = color;
  }

void 
PmMolecule::setModeNumber(const int number) 
  {
  if ((number < 1) || (number > (int)modes.size())) {
    return;
    }

  this->mode_number = number;
  }

void 
PmMolecule::setModeScale(const float scale) {
  this->mode_scale = scale;
  }

//*============================================================*
//*==========          displayModeVector             ==========*
//*============================================================*
// display normal mode vectors for the molecule.

void 
PmMolecule::displayModeVector(const string gname, PmAtomFilter& filter, 
                              const bool show)
  {
  if (!modes.size()) {
    fprintf (stderr, "    >>> no mode data. \n");
    return;
    }

  if (!pmSystem.useGraphics()) return;
  if (!num_chains) buildChains(); 

  int num_coords;
  PmVector3 *coords; 
  string geom_name, atom_name;
  PmGraphicsLine *line;
  PmGraphicsGeometry *geom;
  vector<PmAtom*> chain_atoms;
  PmAtom *atom;
  PmMode mode = modes[mode_number-1];

  #define ndbg_displayModeVectors
  #ifdef dbg_displayModeVectors
  fprintf (stderr, "\n>>>>>> PmMolecule::displayModeVectors: \n");
  fprintf (stderr, "   >>> num chains [%d] \n", num_chains);
  #endif

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    string name;
    #ifdef dbg_displayModeVectors
    fprintf (stderr, "---- chain %c ---- \n", chain->id);
    #endif
    buildGeomName (chain->id, "modeVectors", gname, geom_name);
    getGraphicsGeometry(geom_name, &geom);

    if (!geom) {
      //getChainAtomCoords(chain->id, atom_name, &num_coords, &coords);
      getChainAtoms(chain->id, filter, chain_atoms);

      if (chain_atoms.size()) {
        num_coords = chain_atoms.size();
        coords = new PmVector3[2*num_coords]; 

        for (int j = 0; j < num_coords; j++) {
          atom = chain_atoms[j]; 
          coords[2*j] = atom->pos; 
          coords[2*j+1] = atom->pos + mode_scale*mode.eigen_vectors[atom->index]; 
          }

        line = new PmGraphicsLine(geom_name, 2*num_coords, coords);
        addGraphicsGeometry (line);
        }
      }
    else {
      line = dynamic_cast<PmGraphicsLine*>(geom);
      }

    PmGraphicsAttributes atts;
    //atts.setLineWidth (this->datts.line_width);
    atts.setDisjoint(true);
    atts.setColor(this->mode_color);
    line->setAttributes(atts);
    line->display();
    }
  }

//*============================================================*
//*==========              getMode                   ==========*
//*============================================================*
// get a normal mode for the molecule.

void 
PmMolecule::getMode(const int id, PmMode& mode)
  {
  mode.size = 0;
  mode.eigen_vectors = NULL;
  
  if ((id < 1) || (id > (int)modes.size())) {
    return;
    }

  mode = modes[id-1];
  }

//*============================================================*
//*==========              setModes                  ==========*
//*============================================================*
// set normal mode data for the molecule.                       

void
PmMolecule::setModes(vector<PmMode>& modes) {
  this->modes = modes;
  }

//*============================================================*
//*==========              procQuery                 ==========*
//*============================================================*
// process a query on a molecule.

void 
PmMolecule::procQuery(PmQuery& query)
  {
  //fprintf (stderr, ">>>>>> PmMolecule::procQuery  \"%s\" \n", query.name.c_str());
  PmMolecule *mol;
  string str, mol_name, gtype;
  char chain_id;
  unsigned int n, n1, n2, n3;
  n = query.name.length();
  n1 = query.name.find("[");
  n2 = query.name.find("]");

  if ((n1 == string::npos) || (n2 == string::npos)) {
    return;
    }

  str.assign(query.name, n1+1, n2-n1-1);
  n3 = query.name.find(":");
  gtype.assign(query.name,n2+1,n-n2+1);

  if (n3 == string::npos) {
    mol_name = str;
    chain_id = '\0';
    }
  else {
    mol_name.assign(query.name, n1+1, n3-n1-1);
    chain_id = query.name[n3+1];
    }

  /*
  fprintf (stderr, "   >>> mol name \"%s\" \n", mol_name.c_str());
  fprintf (stderr, "   >>> chain id[%c] \n", chain_id); 
  fprintf (stderr, "   >>> gtype    \"%s\" \n", gtype.c_str());
  */
  pmSystem.getMolecule (mol_name, &mol);

  if (!mol) {
    pmSystem.getDomain(mol_name, &mol);

    if (!mol) {
      return;
      }
    }

  mol->procQuery(query, chain_id, gtype);
  }

//*============================================================*
//*==========              procQuery                 ==========*
//*============================================================*
// process a query on a molecule.

void
PmMolecule::procQuery(PmQuery& query, char chain_id, string gtype)
  {
  #define ndbg_PmMolecule_procQuery
  #ifdef dbg_PmMolecule_procQuery
  fprintf (stderr, ">>>>>> PmMolecule::procQuery  \"%s\" \n", query.name.c_str());
  fprintf (stderr, "   >>> chain = %c \n", chain_id); 
  fprintf (stderr, "   >>> entity = %d \n", query.entity); 
  #endif
  PmChain *chain;
  vector<PmResidue*> rlist;
  PmResidue *res, *res1, *res2, *res3;
  unsigned int entity; 
  char res_name[3];
  string mobj, atom_names;
  PmAtom *atom1, *atom2, *atom3;
  PmVector3 vec;
  float d1, d2, d3;

  // get the chain //

  chain = getChain (chain_id);

  if (!chain) {
    return;
    }

  if (this->parent) {
    mobj = "domain";
    }
  else {
    mobj = "molecule";
    }

  entity = query.entity;

  // process entity depending on geometry picked //

  if (gtype == "backbone") {
    getChainResidues (chain_id, rlist);

    if (entity >= rlist.size()) {
      return;
      }

    res1 = rlist[entity-1];
    res1->getAtom(this->backbone_atom, &atom1);
    vec = atom1->pos - query.point;
    d1 = vec.length();

    res2 = rlist[entity];
    res2->getAtom(this->backbone_atom, &atom2);
    vec = atom2->pos - query.point;
    d2 = vec.length();

    if (entity+1 < rlist.size()) {
      res3 = rlist[entity+1];
      res3->getAtom(this->backbone_atom, &atom3);
      vec = atom3->pos - query.point;
      d3 = vec.length();
      }
    else {
      d3 = 1e6; 
      }

    if (d1 < d2) {
      if (d1 < d3) {
        res = res1;
        }
      else {
        res = res3;
        }
      }
    else {
      if (d2 < d3) {
        res = res2;
        }
      else {
        res = res3;
        }
      }

    pm_CmpdGetName (res->ctype, res_name);
    fprintf (stderr, "    %s %s chain %c  residue %d %s \n", mobj.c_str(), 
             name.c_str(), chain_id, res->id, res_name); 
    }

  else if ((gtype == "atoms") || (gtype == "bonds_atoms")) {
    vector<PmAtom*> chain_atoms;
    PmAtom *atom;
    getChainAtoms(chain_id, atom_names, chain_atoms);

    if (entity >= chain_atoms.size()) {
      return;
      }

    #ifdef dbg_PmMolecule_procQuery
    fprintf (stderr, "   >>> entity = %d \n", entity); 
    #endif

    atom = chain_atoms[entity-1];
    pm_CmpdGetName (atom->ctype, res_name);
    fprintf (stderr, "    %s %s chain %c  residue %d %s  atom %s \n", mobj.c_str(), 
             name.c_str(), chain_id, atom->seq, res_name, atom->name); 
    }

  else if (gtype == "bonds") { 
    vector<PmBond> bonds;
    PmAtom *atom1, *atom2;
    getChainAtomBonds (chain_id, bonds);

    if (entity >= bonds.size()) {
      return;
      }

    atom1 = bonds[entity-1].atom1;
    atom2 = bonds[entity-1].atom2;
    pm_CmpdGetName (atom1->ctype, res_name);

    fprintf (stderr, "    %s %s chain %c  residue %d %s  bond atoms %s %s \n", 
             mobj.c_str(), name.c_str(), chain_id, atom1->seq, res_name, atom1->name,
             atom2->name); 
    }
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a molecule using atom ids.

void
PmMolecule::defineRegion(const string name, const vector<int>& atom_ids) 
  {
  //fprintf (stderr, "\n ######### PmMolecule::defineRegion  name=\"%s\" \n", name.c_str());
  vector<PmAtom> mol_atoms;
  PmAtom atom;
  PmVector3 center; 
  int n, natoms, id;

  PmMolRegion *region = new PmMolRegion;
  region->name = name;
  this->getAtoms(mol_atoms);
  natoms = mol_atoms.size();

  n = 0;
  center.set(0,0,0);

  for (unsigned int i = 0; i < atom_ids.size(); i++) {
    id = atom_ids[i]-1;

    if ((id < 0) || (id >= natoms)) {
      continue;
      }

    atom = mol_atoms[id];
    region->coords.push_back(atom.pos);
    region->coord_indexes.push_back(atom.id);
    region->radii.push_back(atom.getRadius());
    region->atom_names.push_back(atom.name);
    center = center + atom.pos;
    n += 1;
    //fprintf (stderr, "id %d  atom id %d   seq %d  \n", id, atom.id, atom.seq); 
    }

  if (!n) {
    return;
    }

  n = region->coords.size();
  region->center = center / n;
  regions.push_back(region);
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a molecule using ids and coordinates.

void
PmMolecule::defineRegion(const string name, const vector<int>& ids, 
                         const vector<float>& rads, const vector<PmVector3>& coords, 
                         const vector<string>& atom_names) 
  {
  /*
  fprintf (stderr, "\n>>>>>> PmMolecule::defineRegion  name=\"%s\" \n", name.c_str());
  fprintf (stderr, ">>> ids size=%d \n", ids.size()); 
  fprintf (stderr, ">>> coords size=%d \n", coords.size()); 
  fprintf (stderr, ">>> rads size=%d \n", rads.size()); 
  fprintf (stderr, ">>> atom names size=%d \n", atom_names.size()); 
  */
  PmVector3 center; 
  int n, id, num_coords;

  PmMolRegion *region = new PmMolRegion;
  region->name = name;

  n = 0;
  center.set(0,0,0);
  num_coords = coords.size();

  for (unsigned int i = 0; i < ids.size(); i++) {
    id = ids[i]-1;

    if ((id < 0) || (id >= num_coords)) {
      continue;
      }
      
    region->coords.push_back(coords[id]);
    region->coord_indexes.push_back(id);
    region->radii.push_back(rads[id]);

    if (!atom_names.empty()) {
      region->atom_names.push_back(atom_names[id]);
      }

    center = center + coords[id]; 
    n += 1;
    //fprintf (stderr, ">>> add atom=\"%s\" \n", atom_names[id].c_str());
    }

  if (!n) {
    return;
    }

  n = region->coords.size();
  region->center = center / n;
  regions.push_back(region);
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a molecule.

void 
PmMolecule::defineRegion(const string name, PmAtomFilter& filter, const string desc,
                         PmMolRegionParameters& params)
  {
  /*
  fprintf (stderr, "\n>>>>>> PmMolecule::defineRegion  name=\"%s\" \n", name.c_str());
  fprintf (stderr, "   >>> desc=%s \n", desc.c_str()); 
  fprintf (stderr, "   >>> use surface[%d] \n", params.use_surface); 
  fprintf (stderr, "   >>> use sidechains[%d] \n", params.use_sidechains); 
  fprintf (stderr, "   >>> no  mainchain[%d] \n", params.no_mainchain); 
  */

  // define the region from a surface  //

  if (params.use_surface) {
    defineSurfRegion(name, filter, desc, params);
    return;
    }
  else if (params.use_pca) {
    definePcaRegion(name, filter, desc, params);
    return;
    }

  // get the residues specified in <desc>  //

  vector<PmResidue*> rlist;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    pm_ErrorReport (PM, "no residues found for region \"%s\".", "*", name.c_str());
    return;
    }

  //fprintf (stderr, ">>> rlist size=%d \n", rlist.size()); 
  PmMolRegion *region = new PmMolRegion;
  region->name = name;
  region->descriptor = desc;

  // compute the center of the residues and create //
  // a list of atom indexes.                       //

  PmResidue *res;
  PmVector3 center, res_center, pos;
  unsigned int n = 0;
  center.set(0,0,0);
  float r, max_r, min_r;
  vector<PmVector3> res_coords; 

  //====== use centers of sidechains  ======//

  if (params.use_sidechains) {
    for (unsigned int i = 0; i < rlist.size(); i++) {
      res = rlist[i];
      res_center.set(0,0,0);
      res_coords.clear();

      for (int j = 0; j < res->num_atoms; j++) {
        PmAtom *atom = res->atoms[j];

        if (filter.hasName(atom->name)) {
          res_coords.push_back(atom->pos);
          res_center = res_center + atom->pos;
          }
        }

      n = res_coords.size();

      if (n) {
        res_center = res_center / n;
        max_r = 0.0;
        min_r = 1e6;

        for (unsigned int j = 0; j < res_coords.size(); j++) {
          pos = res_center - res_coords[j];
          r = pos.length();
          if (r > max_r) max_r = r;
          if (r < min_r) min_r = r;
          }

        region->coords.push_back(res_center);
        center = center + res_center;
        //r = (max_r + min_r) / 2.0;
        r = max_r;
        region->radii.push_back(r);
        }
      }

    n = region->coords.size();

    if (!n) {
      pm_ErrorReport (PM, "no atoms found for region \"%s\".", "*", name.c_str());
      return;
      }
    }

  // use specified atoms of each residue  //

  else {
    for (unsigned int i = 0; i < rlist.size(); i++) {
      res = rlist[i];

      for (int j = 0; j < res->num_atoms; j++) {
        PmAtom *atom = res->atoms[j];
        //fprintf (stderr, ">>> atom name=%s \n", atom->name); 

        if (filter.hasName(atom->name)) {
          //fprintf (stderr, "    atom name in filter \n"); 
          region->coords.push_back(atom->pos);
          region->coord_indexes.push_back(atom->id);
          region->radii.push_back(atom->getRadius());
          region->atom_names.push_back(atom->name);
          center = center + atom->pos;
          n += 1;
          //fprintf (stderr, ">>> atom=%s ", atom->name);
          //fprintf (stderr, "atom id %d   seq %d  \n", atom->id, atom->seq);
          }
        }
      }

    if (!n) {
      pm_ErrorReport (PM, "no atoms found for region \"%s\".", "*", name.c_str());
      return;
      }
    }

  n = region->coords.size();
  region->center = center / n;

  // add to regions defined for the molecule //

  regions.push_back(region);
  }

//*============================================================*
//*==========           defineSurfRegion             ==========*
//*============================================================*
// define a surface region for a molecule.

void
PmMolecule::defineSurfRegion(const string name, PmAtomFilter& filter, const string desc,
                             PmMolRegionParameters& params)
  {

  if (!surface) {
    return;
    }

  // get the residues specified in <desc>  //

  vector<PmResidue*> rlist;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    return;
    }

  float tol, dist, r;
  PmResidue *res;
  PmVector3 center, coord, v;
  vector<PmVector3> coords;
  PmAtomFilter afilter;

  PmMolRegion *region = new PmMolRegion;
  region->name = name;
  region->descriptor = desc;

  // use only sidechain atoms for test //

  getBackboneAtomNames(afilter.names);
  afilter.exclude = true;

  // get the surface coordinats;

  //surface->getCoordinates(coords);
  surface->getCellCenters(coords);

  // find the surface coordinates closest //
  // to each atonm.                       //

  center.set(0,0,0);
  tol = params.tolerance*params.tolerance;

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    for (int j = 0; j < res->num_atoms; j++) {
      PmAtom *atom = res->atoms[j];

      if (afilter.hasName(atom->name)) {
        r = atom->getRadius();

        for (unsigned int k = 0; k < coords.size(); k++) {
          coord = coords[k];
          //fprintf (stderr, "   >>> coord %g %g %g \n", coord[0], coord[1], coord[2]); 
          v = coord - atom->pos;
          dist = v * v; 

          if (fabs(dist - r) < tol) {
            region->coords.push_back(coord);
            region->coord_indexes.push_back(k);
            region->radii.push_back(0.05);
            center = center + coord;
            }
          }
        }
      }
    }

  int n = region->coords.size();
  //fprintf (stderr, "   >>> num points[%d] \n", n); 

  if (!n) {
    return;
    }

  region->center = center / n;

  // add to regions defined for the molecule //

  regions.push_back(region);
  }

//*============================================================*
//*==========           definePcaRegion              ==========*
//*============================================================*
// define a region for a molecule using principal components.

void
PmMolecule::definePcaRegion(const string name, PmAtomFilter& filter, const string desc,
                            PmMolRegionParameters& params)
  {
  PmPcaResults pca;
  PmVector3 pt;
  PmAtomFilter pca_filter;

  pca_filter.names.push_back("CA");
  this->compPrincipalComponents(pca_filter, pca);

  PmMolRegion *region = new PmMolRegion;
  region->name = name;
  region->descriptor = desc;
  region->center = pca.com;

  region->coords.push_back(pca.com);
  region->coord_indexes.push_back(0);
  region->radii.push_back(0.05);

  pt = pca.com + pca.s1*pca.axis1;
  region->coords.push_back(pt);
  region->coord_indexes.push_back(1);
  region->radii.push_back(0.05);


  // add to regions defined for the molecule //

  regions.push_back(region);
  }

//*============================================================*
//*==========              displayRegion             ==========*
//*============================================================*
// display a region for a molecule.

void 
PmMolecule::displayRegion(const string name, PmVector3 color, PmMoleculeRenderType rtype,
                          bool use_spheres)
  {
  PmRegion *rgn;

  getRegion(name, &rgn);

  if (!rgn) {
    return;
    }

  string geom_name;
  PmGraphicsPoint *points;
  PmGraphicsSphere *spheres;
  PmGraphicsGeometry *geom;
  int num;
  PmVector3 *verts;
  float *rads;
  PmGraphicsAttributes atts;

  num = rgn->coords.size();

  if (!num) {
    return;
    }

  buildGeomName ('\0', "region", rgn->name, geom_name);
  getGraphicsGeometry(geom_name, &geom);

  if (!geom) {
    verts = new PmVector3[num];

    for (int i = 0; i < num; i++) {
      verts[i] = rgn->coords[i];
      }

    if (use_spheres) {
      rads = new float[num];

      for (int i = 0; i < num; i++) {
        rads[i] = rgn->radii[i];
        }

      spheres = new PmGraphicsSphere(geom_name, num, verts, rads);
      addGraphicsGeometry (spheres);
      geom = spheres;
      }
    else {
      points = new PmGraphicsPoint(geom_name, num, verts);
      addGraphicsGeometry (points);
      geom = points;
      }
    }
  else {
    points = dynamic_cast<PmGraphicsPoint*>(geom);
    }

  atts.setDisplayType(PmGraphicsInterface::mapRenderType(rtype));
  atts.setMarker(false);
  atts.setScale(2.0);
  atts.setColor(color);

  geom->setAttributes(atts);
  geom->display();
  }

//*============================================================*
//*==========              getRegion                 ==========*
//*============================================================*
// get a region for a molecule.

void 
PmMolecule::getRegion(const string name, PmRegion **rgn)
  {
  fprintf (stderr, "\n>>>>>> PmMolecule::getRegion [%s] \n", name.c_str());
  *rgn = NULL;

  if (!regions.size()) {
    return;
    }

  for (unsigned int i = 0; i < regions.size(); i++) {
    if (regions[i]->name == name) {
      *rgn = regions[i];
      return;
      }
    }
  }

//*============================================================*
//*==========              hasRegion                 ==========*
//*============================================================*
// check if a region is defined for a molecule.

bool
PmMolecule::hasRegion(const string name)
  {
  if (!regions.size()) {
    return false;
    }

  for (unsigned int i = 0; i < regions.size(); i++) {
    if (regions[i]->name == name) {
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========              hasResidue                ==========*
//*============================================================*
// check if a molecule has the given residue. 

bool 
PmMolecule::hasResidue(const char chain_id, const int id)
  {
  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *first, *last; 

  buildChains();

  // get the residues for the given chain id //

  getChainResidues (chain_id, rlist);
  num_res = rlist.size();

  if (num_res == 0) {
    return false;
    }

  first = rlist[0];
  last = rlist[num_res-1];

  if ((id < first->id) || (id > last->id)) {
    return false;
    }

  for (unsigned int i = 0; i < rlist.size(); i++) {
    if (rlist[i]->id == id) {
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========            getResidueBoundSpheres      ==========*
//*============================================================*
// get the bounding spheres for each residue.

void 
PmMolecule::getResidueBounds(vector<PmVector3>& center, vector<float>& radius) 
  {
  string desc;
  PmResidue *res;
  PmAtomFilter filter;
  PmAtom *atom;
  vector<PmAtom*> res_atoms;
  PmVector3 res_center, v;
  float rad;
  int num_res;
  vector<PmResidue*> rlist;
  int n;
  PmExtent extent;
  float dx, dy, dz;

  // get the residues for this molecule //

  this->getResidues(desc, rlist);

  // just use sidechain atoms //

  this->getBackboneAtomNames(filter.names);
  filter.exclude = true;
  num_res = rlist.size();
  center.clear();
  radius.clear(); 

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];
    res_center.set(0,0,0);
    n = 0;
    extent.min[0] =  1e6; extent.min[1] =  1e6; extent.min[2] = 1e6;
    extent.max[0] = -1e6; extent.max[1] = -1e6; extent.max[2] = -1e6;

    for (int ia = 0; ia < res->num_atoms; ia++) {
      atom = res->atoms[ia];

      if (filter.hasName(atom->name)) {
        res_atoms.push_back(atom);
        res_center = res_center + atom->pos;
        extent.update(atom->pos);
        n += 1;
        }
      }
 
    res_center = res_center / n;
    dx = extent.max[0] - extent.min[0];
    dy = extent.max[1] - extent.min[1];
    dz = extent.max[2] - extent.min[2];

    if (dx > dy) {
      rad = dx;
      }
    else {
      rad = dy;
      }

    if (dz > rad) {
      rad = dz;
      }

    radius.push_back(rad / 2.0 + 0.2);
    center.push_back(res_center);
    }
  }

//*============================================================*
//*==========              getResidues               ==========*
//*============================================================*
// get the residues for a molecule from a string descriptor.
// the desciptor will be of the form: <chain>[<reslist>]

void 
PmMolecule::getResidues(const string dstr, vector<PmResidue*>& rlist)
  {
  #define ndbg_getResidues  
  #ifdef dbg_getResidues  
  fprintf (stderr, "\n>>>>>> PmMolecule::getResidues desc [%s] \n", dstr.c_str());
  #endif
  rlist.clear();

  // build residue list //

  buildChains();

  // if no descriptor return all residues //

  if (dstr.empty()) {
    for (unsigned int i = 0; i < residues.size(); i++) {
      rlist.push_back(&residues[i]);
      }
    return;
    }

  // check for helix //

  if (!dstr.compare(0,5,"helix")) {
    getHelixResidues(dstr, rlist);
    return;
    }

  // check for sheet //

  if (!dstr.compare(0,5,"sheet")) {
    getSheetResidues(dstr, rlist);
    return;
    }

  // check for loop //

  if (!dstr.compare(0,4,"loop")) {
    getLoopResidues(dstr, rlist);
    return;
    }

  char chain_id, *s, *desc;
  PmResidue *res;
  vector<PmResidue*> chain_rlist;
  PmChain *chain;

  // first check that all chain ids are valid //

  s = (char*)dstr.c_str();

  while (s && (*s != '\n') && (*s != ' ') && (*s != '\0')) {
    if (isalpha(*s)) {
      chain = getChain (*s);

      if (!chain) {
        #ifdef dbg_getResidues  
        fprintf (stderr, "   >>> chain id %c not valid. \n", *s); 
        #endif
        return;
        }
      }

    s++;
    }

  // parse and process domain descriptor. the descriptor will be of //
  // the form A[2-10]B[5-10].                                       //

  desc = (char*)dstr.c_str();

  while (desc && (*desc != '\n') && (*desc != ' ') && (*desc != '\0')) {
    chain_id = *desc;
    desc++;
    chain = getChain (chain_id);
    #ifdef dbg_getResidues  
    fprintf (stderr, "  >>>> chain id [%c] \n", chain_id);
    fprintf (stderr, "  >>>> desc = %c \n", *desc);
    #endif

    if (!chain) {
      return;
      }

    if (*desc == '[') {
      *desc++;
      getChainResidues (chain_id, &desc, chain_rlist);
      }
    /*  note: remove syntax <chain><res> since <chain> can be integer.
    else if (isdigit(*desc)) {
      getChainResidues (chain_id, &desc, chain_rlist);
      }
    */
    else {
      getChainResidues (chain_id, chain_rlist);
      }

    #ifdef dbg_getResidues  
    fprintf (stderr, "  >>>> num res [%d] \n", chain_rlist.size());
    #endif

    for (unsigned int i = 0; i < chain_rlist.size(); i++) {
      res = chain_rlist[i];
      //fprintf (stderr, "   #### res id[%d]  mtype[%d]\n", res->id, res->mtype);
      rlist.push_back(res);
      }
    }
  }

//*============================================================*
//*==========              getResidueCoords          ==========*
//*============================================================*
// get the center of mass for the residues specified in the 
// string res.

bool 
PmMolecule::getResidueCoords (const string res, const string sel, PmVector3& jpos)
  {
  #define ndbg_getResidueCoords 
  #ifdef dbg_getResidueCoords 
  fprintf (stderr, ">>>>>> PmMolecule::getResidueCoords  name=%s \n", name.c_str());
  fprintf (stderr, ">>> res=\"%s\" \n", res.c_str());
  fprintf (stderr, ">>> sel=\"%s\" \n", sel.c_str());
  #endif

  jpos.set(0,0,0);
  int n = res.size();

  if (n == 0) {
    return false;
    }

  char *str = new char[n+1];
  memset(str, '\0', n+1);
  res.copy(str, n);
  char *desc = str;

  vector<PmResidue*> rlist;
  char chain_id = *desc++;
  PmChain *chain = getChain(chain_id);

  if (!chain) {
    #ifdef dbg_getResidueCoords 
    fprintf (stderr, ">>> no chain %c \n", chain_id); 
    #endif
    return false;
    }

  desc++;
  getChainResidues (chain_id, &desc, rlist);

  #ifdef dbg_getResidueCoords 
  fprintf (stderr, ">>> rlist size [%d] \n", rlist.size());
  #endif

  if (rlist.size() == 0) {
    return false;
    }

  n = 0;
  int num_res = rlist.size();
  PmVector3 pos(0,0,0);

  // compute center of residues  //

  for (int i = 0; i < num_res; i++) {
    PmResidue *res = rlist[i];
    int num_atoms = res->num_atoms;
    PmAtom *atom;

    #ifdef dbg_getResidueCoords 
    for (int j = 0; j < num_atoms; j++) {
      fprintf (stderr, ">>> atom=%s \n", res->atoms[j]->name); 
      }
    #endif

    if (sel != "") {
      res->getAtom (sel, &atom);

      if (atom) {
        pos = pos + atom->pos;
        n += 1;
        }
      }
    else {
      for (int j = 0; j < num_atoms; j++) {
        atom = res->atoms[j];
        pos = pos + atom->pos;
        n += 1;
        }
      }
    }

  if (!n) {
    delete[] str;
    return false;
    }

  delete[] str;
  jpos = (1.0/n)*pos;
  return true;
  }

//*============================================================*
//*==========              getNumResidues            ==========*
//*============================================================*
// get the number of residues.

int
PmMolecule::getNumResidues() {
  buildChains();
  return (num_residues);
  }

//*============================================================*
//*==========                Sheet                   ==========*
//*============================================================*
// add a sheet.

void
PmMolecule::addSheet(PmSheet& sheet) {
  sheets.push_back(sheet);
  }

void
PmMolecule::getSheets(vector<PmSheet>& sheets) {
  sheets = this->sheets;
  }

//*============================================================*
//*==========                getSheets               ==========*
//*============================================================*
// get sheets using id.

bool
PmMolecule::getSheets(const string id, const int num, vector<PmSheet>& sheet_list)
  {

  // first search for sheets by id //

  for (unsigned int i = 0; i < sheets.size(); i++) {
    if (id == sheets[i].id) {
      if (num) {
        if (num == sheets[i].num) {
          sheet_list.push_back(sheets[i]);
          }
        }
      else {
        sheet_list.push_back(sheets[i]);
        }
      }
    }

  if (sheet_list.size()) {
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getSheetResidues          ==========*
//*============================================================*
// get the resides for a molecule from a sheet descriptor.

void
PmMolecule::getSheetResidues(const string dstr, vector<PmResidue*>& rlist)
  {
  //fprintf (stderr, "\n>>>>>> PmMolecule::getSheetResidues desc [%s] \n", dstr.c_str());
  char *desc;
  string sheet_id, num_str, res_desc;
  vector<PmSheet> sheet_list;
  stringstream dss;
  bool has_num = false;
  int n, sheet_num;

  // parse and process sheet descriptor. the descriptor  //
  // will be of the form sheet:<id>                      //

  desc = (char*)dstr.c_str();
  for (int i = 0; i < 5; i++) desc++;

  while (desc && (*desc != '\n') && (*desc != ' ') && (*desc != '\0')) {
    if (*desc == ':') { 
      has_num = true;
      }
    else if (has_num) { 
      num_str.push_back(*desc);
      }
    else {
      sheet_id.push_back(*desc);
      }

    desc++;
    }

  if (has_num) {
    sscanf (num_str.c_str(), "%d", &sheet_num); 
    //fprintf (stderr, ">>> num  = %d \n", sheet_num);
    }
  else {
    sheet_num = 0;
    }

  if (!getSheets(sheet_id, sheet_num, sheet_list)) {
    pm_ErrorReport (PM, "sheet \"%s\" not found.", "*", sheet_id.c_str());
    return;
    }

  n = sheet_list.size();
  dss << sheet_list[0].init_chain_id << '['; 

  for (int i = 0; i < n; i++) {
    dss << sheet_list[i].init_seq_num << '-' << sheet_list[i].term_seq_num; 

    if (i != n-1) {
      dss << ','; 
      }
    }

  dss << ']'; 
  res_desc = dss.str();
  //fprintf (stderr, "    >>> sheet residues %s \n", res_desc.c_str()); 
  getResidues(res_desc, rlist);
  }

//*============================================================*
//*==========              getSidechainGeometry      ==========*
//*============================================================*
// get the center of mass and radius for sidechains.           

void 
PmMolecule::getSidechainGeometry(const string desc, vector<PmVector3>& coords,
                                 vector<float>& radii)
  {
  vector<PmResidue*> rlist; 
  PmResidue *res; 
  PmVector3 center;
  PmExtent extent;
  PmAtom *atom;
  PmAtomFilter filter;
  int n;
  float dx, dy, dz;

  // get residues for descriptor //

  getResidues(desc, rlist);

  if (!rlist.size()) {
    return;
    }

  // set exlude backbone atoms //

  filter.exclude = true;
  this->getBackboneAtomNames(filter.names);

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];
    center.set(0,0,0);
    n = 0;
    extent.min[0] =  1e6; extent.min[1] =  1e6; extent.min[2] = 1e6;
    extent.max[0] = -1e6; extent.max[1] = -1e6; extent.max[2] = -1e6;

    for (int ia = 0; ia < res->num_atoms; ia++) {
      atom = res->atoms[ia];

      if (filter.hasName(atom->name)) {
        center = center + atom->pos;
        extent.update(atom->pos);
        n += 1;
        }
      }

    center = center / n;
    dx = extent.max[0] - extent.min[0];
    dy = extent.max[1] - extent.min[1];
    dz = extent.max[2] - extent.min[2];

    if (dx > dy) {
      radii.push_back(dx);
      }
    else {
      radii.push_back(dy);
      }

    if (dz > radii[i]) {
      radii[i] = dz;
      }

    radii[i] = radii[i] / 2.0 + 0.2;
    coords.push_back(center);
    }
  }

//*============================================================*
//*==========              displayTopology           ==========*
//*============================================================*
// display topology of secondary structures.

void
PmMolecule::displayTopology(bool show_chain, bool show_contact, PmVector3& color, 
                            bool show)
  {
  //fprintf (stderr, "\n>>>>>> PmMolecule::displayTopology \n");
  PmTopology *topo;;

  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *first_res, *last_res;
  PmAtomFilter filter;
  string desc;
  stringstream dss;

  int num_chains;
  vector<PmChain*> chains;
  PmChain *chain;
  char chain_id;

  PmStructureList *slist; 
  int init_snum, term_snum;
  //char stype[] = { 'H', 'L', 'S' };

  PmVector3 pos;
  vector<PmVector3> vertices; 
  string geom_name, atom_name;
  PmGraphicsLine *line;
  PmGraphicsGeometry *geom;
  int num_verts;
  PmVector3 *gverts;
  PmGraphicsAttributes atts;

  if (this->parent) {
    this->parent->getTopology(&topo);
    }
  else {
    this->getTopology(&topo);
    }

  this->getChains(chains);
  num_chains = chains.size();
  fprintf (stderr, "   >>> num chains = %d \n", num_chains);

  // show chain topology for the domain //

  if (show_chain) {
    for (int i = 0; i < num_chains; i++) {
      chain = chains[i];
      chain_id = chain->id;

      this->getChainResidues (chain_id, rlist);
      num_res = rlist.size();
      first_res = rlist[0];
      last_res = rlist[num_res-1];
      vertices.clear();

      for (slist = topo[i].list; slist != NULL; slist = slist->next) {
        init_snum = slist->init_seq_num;
        term_snum = slist->term_seq_num;

        if ((init_snum < first_res->id) || (term_snum > last_res->id)) {
          continue;
          }

        /*
        fprintf (stderr, "   >>> %c%s %d-%d \n", stype[slist->type], slist->id.c_str(),
                 init_snum, term_snum); 
        */

        // for a loop add one residue to the end to  //
        // ensure the domains overlap by one residue //

        if (slist->type == PM_SECONDARY_STRUCTURE_LOOP) {
          if (init_snum != first_res->id) init_snum -= 1;
          if (term_snum != last_res->id)  term_snum += 1;
          }

        dss << chain_id << '[' << init_snum << ']';
        desc = dss.str();
        dss.str(std::string());
        this->getResidueCoords(desc, "CA", pos);
        vertices.push_back(pos);

        dss << chain_id << '[' << term_snum << ']';
        desc = dss.str();
        dss.str(std::string());
        this->getResidueCoords(desc, "CA", pos);
        vertices.push_back(pos);
        }

      if (vertices.size()) {
        num_verts = vertices.size();
        gverts = new PmVector3[num_verts];
    
        for (int i = 0; i < num_verts; i++) {
          gverts[i] = vertices[i];
          }

        buildGeomName (chain_id, "chainTopology", "", geom_name);
        getGraphicsGeometry(geom_name, &geom);

        if (!geom) {
          line = new PmGraphicsLine(geom_name, num_verts, gverts);
          addGraphicsGeometry (line);
          }
        else {
          line = dynamic_cast<PmGraphicsLine*>(geom);
          }

        atts.setColor(color);
        atts.disjoint = true;
        line->setAttributes(atts);
        line->display();
        }
      }
    }


  // show contact topology for the domain //

  if (show_contact) {
    PmStructureListPtr **contact, *cnlist;
    PmStructureList *s1, *s2;
    PmVector3 pos1;
    int rid;

    for (int i = 0; i < num_chains; i++) {
      chain = chains[i];
      chain_id = chain->id;

      this->getChainResidues (chain_id, rlist);
      num_res = rlist.size();
      first_res = rlist[0];
      last_res = rlist[num_res-1];
      vertices.clear();

      contact = topo[i].contact;
      s1 = topo[i].list;

      for (int j = 0; j < topo[i].num_structures; j++, s1 = s1->next) {
        //fprintf (stderr, "<----%c%s----> \n", stype[s1->type], s1->id.c_str());
        //fprintf (stderr, "<----%c %c[%d-%d]----> \n", stype[s1->type], chain_id, 
        //         s1->init_seq_num, s1->term_seq_num); 
        init_snum = s1->init_seq_num;
        term_snum = s1->term_seq_num;
        //dss << chain_id << '[' << init_snum << '-' << term_snum << ']';
        rid = (init_snum + term_snum) / 2;
        dss << chain_id << '[' << rid << ']';
        desc = dss.str();
        dss.str(std::string());

        if (this->getResidueCoords(desc, "CA", pos1)) {
          cnlist = contact[j];
          //fprintf (stderr, "   >>> %s pos1 = %g %g %g \n", desc.c_str(), pos1[0], 
          //         pos1[1], pos1[2]); 

          while (cnlist) {
            s2 = cnlist->structure;
            //fprintf (stderr, "%c%s\n", stype[s2->type], s2->id.c_str());
            init_snum = s2->init_seq_num;
            term_snum = s2->term_seq_num;
            //dss << chain_id << '[' << init_snum << '-' << term_snum << ']';
            rid = (init_snum + term_snum) / 2;
            dss << chain_id << '[' << rid << ']';
            desc = dss.str();
            dss.str(std::string());

            if (this->getResidueCoords(desc, "CA", pos)) {
              vertices.push_back(pos1);
              vertices.push_back(pos);
              //fprintf (stderr, "   >>> pos = %g %g %g \n", pos[0], pos[1], pos[2]); 
              //fprintf (stderr, "   >>> %s pos = %g %g %g \n", desc.c_str(), pos[0], 
              //         pos[1], pos[2]); 
              }

            cnlist = cnlist->next;
            }
          }
        }

      if (vertices.size()) {
        num_verts = vertices.size();
        gverts = new PmVector3[num_verts];

        for (int i = 0; i < num_verts; i++) {
          gverts[i] = vertices[i];
          }

        buildGeomName (chain_id, "contactTopology", "", geom_name);
        getGraphicsGeometry(geom_name, &geom);

        if (!geom) {
          line = new PmGraphicsLine(geom_name, num_verts, gverts);
          addGraphicsGeometry (line);
          }
        else {
          line = dynamic_cast<PmGraphicsLine*>(geom);
          }

        atts.setColor(color);
        atts.disjoint = true;
        line->setAttributes(atts);
        line->display();
        }
      }
    }
  }

//*============================================================*
//*==========              getTopology               ==========*
//*============================================================*
// build topology of secondary structures.

void
PmMolecule::getTopology(PmTopology **topo)
  {
  if (!this->topology) {
    this->buildTopology();
    }

  *topo = this->topology;
  }

//*============================================================*
//*==========              buildTopology             ==========*
//*============================================================*
// build topology of secondary structures.

void 
PmMolecule::buildTopology()
  {
  char chain_id, str[20];
  vector<PmChain*> clist;
  PmHelix helix;
  PmSheet sheet;
  PmLoop loop;
  PmStructureList *slist, *sptr, *last;
  vector<PmResidue*> rlist;
  int n, init_seq_num, term_seq_num, next_seq_num, res_id;
  bool add_loop;

  // there may be no secondary structures defined //

  if ((helices.size() + sheets.size()) == 0) {
    return;
    }


  // initialize topology for each chain //

  getChains (clist);
  if (!clist.size()) return; 
  topology = new PmTopology[num_chains];

  for (int i = 0; i < num_chains; i++) {
    topology[i].chain_id = clist[i]->id;
    topology[i].list = NULL; 
    topology[i].num_structures = 0; 
    }


  // add helices //

  for (unsigned int i = 0; i < helices.size(); i++) {
    helix = helices[i];
    //fprintf (stderr, ">>> helix [%s] \n", helix.id.c_str());

    for (int j = 0; j < num_chains; j++) {
      if (topology[j].chain_id == helix.init_chain_id) {
        n = j;
        break;
        }
      }

    slist = topology[n].list; 
    last = NULL;

    while (slist) {
      if (helix.init_seq_num < slist->init_seq_num) {
        break;
        }

      last = slist;
      slist = slist->next;
      }

    sptr = new PmStructureList;
    sptr->id = helix.id; 
    sptr->type = PM_SECONDARY_STRUCTURE_HELIX;
    sptr->init_seq_num = helix.init_seq_num; 
    sptr->term_seq_num = helix.term_seq_num; 
    topology[n].num_structures += 1; 

    if (last) {
      sptr->next = slist; 
      last->next = sptr;
      }
    else {
      sptr->next = topology[n].list;
      topology[n].list = sptr; 
      }
    }

  // add sheets //

  for (unsigned int i = 0; i < sheets.size(); i++) {
    sheet = sheets[i];
    //fprintf (stderr, ">>> sheet [%s] \n", sheet.id.c_str());

    for (int j = 0; j < num_chains; j++) {
      if (topology[j].chain_id == sheet.init_chain_id) {
        n = j;
        break;
        }
      }

    slist = topology[n].list;
    last = NULL;

    while (slist) {
      if (sheet.init_seq_num < slist->init_seq_num) {
        break;
        }

      last = slist;
      slist = slist->next;
      }

    sptr = new PmStructureList;
    sptr->id = sheet.id; 
    sptr->type = PM_SECONDARY_STRUCTURE_SHEET;
    sptr->init_seq_num = sheet.init_seq_num;
    sptr->term_seq_num = sheet.term_seq_num;
    topology[n].num_structures += 1; 

    if (last) {
      sptr->next = slist;
      last->next = sptr;
      }
    else {
      sptr->next = topology[n].list;
      topology[n].list = sptr;
      }
    }

  // generate loops //

  int num_loops = 0;

  for (int j = 0; j < num_chains; j++) {
    chain_id = topology[j].chain_id;
    getChainResidues (chain_id, rlist);
    slist = topology[j].list;

    if (slist->init_seq_num != rlist[0]->id) {
      num_loops += 1;
      loop.id = "1"; 
      loop.num = 1; 
      loop.init_chain_id = chain_id; 
      loop.init_seq_num = rlist[0]->id; 
      loop.term_chain_id = chain_id; 
      loop.term_seq_num = slist->init_seq_num-1; 
      addLoop(loop);

      sptr = new PmStructureList;
      sptr->id = loop.id; 
      sptr->type = PM_SECONDARY_STRUCTURE_LOOP;
      sptr->init_seq_num = loop.init_seq_num;
      sptr->term_seq_num = loop.term_seq_num;
      topology[j].num_structures += 1; 

      sptr->next = topology[j].list;
      topology[j].list = sptr;
      }

    while (slist->next) {
      if (slist->term_seq_num != slist->next->init_seq_num-1) {
        init_seq_num = -1; 
        add_loop = false;
        next_seq_num = slist->next->init_seq_num;

        // find first residue //

        for (unsigned int k = 0; k < rlist.size(); k++) {
          if (rlist[k]->id == slist->term_seq_num) { 
            n = k+1;
            break;
            }
          }

        for (unsigned int k = n; k < rlist.size(); k++) { 
          res_id = rlist[k]->id;

          if (res_id != next_seq_num) { 
            if (init_seq_num == -1) {
              init_seq_num = res_id; 
              term_seq_num = res_id; 
              }
            else if (res_id > term_seq_num+1) {
              add_loop = true;
              }
            else {
              term_seq_num += 1;
              }
            }
          else { 
            add_loop = true;
            }

          if (add_loop) { 
            sprintf (str, "%d", num_loops); 
            loop.id.clear(); 
            loop.id = str; 
            loop.num = num_loops+1; 
            loop.init_seq_num = init_seq_num;
            loop.term_seq_num = term_seq_num;
            loop.init_chain_id = chain_id; 
            loop.term_chain_id = chain_id; 
            addLoop(loop);

            sptr = new PmStructureList;
            sptr->id = loop.id; 
            sptr->type = PM_SECONDARY_STRUCTURE_LOOP;
            sptr->init_seq_num = loop.init_seq_num;
            sptr->term_seq_num = loop.term_seq_num;
            topology[j].num_structures += 1; 

            sptr->next = slist->next;
            slist->next = sptr;
            slist = slist->next;
            num_loops += 1;

            if (res_id == next_seq_num) { 
              break;
              }

            init_seq_num = res_id; 
            term_seq_num = res_id; 
            }
          }
        }

      slist = slist->next;
      }
    }


  // print topology //

#define dbg_PmMolecule_buildTopology
#ifdef dbg_PmMolecule_buildTopology

  fprintf (stderr, "\n------- topology -------\n");
  for (int j = 0; j < num_chains; j++) {
    chain_id = topology[j].chain_id;
    slist = topology[j].list;
    fprintf (stderr, "\n>>>> chain [%c]: ", chain_id); 

    while (slist) {
      if (slist->type == PM_SECONDARY_STRUCTURE_HELIX) {
        fprintf (stderr, "H[%d:%d]->", slist->init_seq_num, slist->term_seq_num);
        //fprintf (stderr, "H[%s]->", slist->helix->id.c_str()); 
        }
      else if (slist->type == PM_SECONDARY_STRUCTURE_SHEET) {
        //fprintf (stderr, "s[%s]->", slist->sheet->id.c_str()); 
        fprintf (stderr, "S[%d:%d]->", slist->init_seq_num, slist->term_seq_num);
        }
      else if (slist->type == PM_SECONDARY_STRUCTURE_LOOP) { 
        fprintf (stderr, "L[%d:%d]->", slist->init_seq_num, slist->term_seq_num);
        }

      slist = slist->next;
      }

    fprintf (stderr, "\n");
    }
#endif


#define compContactTopo_PmMolecule_buildTopology
#ifdef compContactTopo_PmMolecule_buildTopology

  // now compute contact topology //

  PmResidueContact *res_con;
  PmStructureListPtr **contact, *cnlist, *cnptr;
  PmStructureList *s1, *s2;
  vector<PmResidue*> rlist1, rlist2;
  stringstream dss;
  string desc;
  //fprintf (stderr, "\n---------- compute contact ---------- \n");

  PmResidue *res;
  PmAtomFilter filter;
  int init_snum1, init_snum2, term_snum1,  term_snum2;
  //char stype[] = { 'H', 'L', 'S'};
  char res_chain;
  int scount;

  // compute residue-residue contact for this molecule //

  this->getResContact();
  res_con = this->residue_contacts;
  this->getResidues("", rlist);

  // use contact information to compute contact topology //

  for (int i = 0; i < num_chains; i++) {
    chain_id = topology[i].chain_id;
    s1 = topology[i].list;
    contact = new PmStructureListPtr*[topology[i].num_structures];
    //fprintf (stderr, "\n>>>> chain %c (%d): ", chain_id, topology[i].num_structures); 

    for (int j = 0; j < topology[i].num_structures; j++) {
      contact[j] = NULL;
      }

    /*
    for (int j = 0; j < rlist.size(); j++) {
      if (rlist[j]->getChainId() == chain_id) {
        res_con = &this->residue_contacts[j];
        break;
        }
      }
    */

    scount = 0;

    while (s1) {
      init_snum1 = s1->init_seq_num;
      term_snum1 = s1->term_seq_num;
      //fprintf (stderr, "<----%c%s----> \n", stype[s1->type], s1->id.c_str()); 
      //fprintf (stderr, "<----%c[%d-%d]----> \n", chain_id, s1->init_seq_num, s1->term_seq_num); 
      n = -1;

      for (unsigned int j = 0; j < rlist.size(); j++) {
        if (rlist[j]->id == init_snum1) {
          res_chain = rlist[j]->getChainId();

          if (res_chain == chain_id) { 
            n = j;
            break;
            }
          }
        }

      /*
      for (int j = 0; j < res_con[n].num; j++) {
        res = res_con[n].list[j];
        fprintf (stderr, "%d,", res->id); 
        }
      */

      for (s2 = topology[i].list; s2 != NULL && n != -1; s2 = s2->next) {
        if (s1 == s2) {
          continue;
          }

        init_snum2 = s2->init_seq_num;
        term_snum2 = s2->term_seq_num;

        for (int j = 0; j < res_con[n].num; j++) {
          res = res_con[n].list[j];
          res_chain = res->getChainId();

          if ((res->id >= init_snum2) && (res->id <= term_snum2) && 
              (res_chain == chain_id)) {
            //fprintf (stderr, "%c%s\n", stype[s2->type], s2->id.c_str()); 
            //fprintf (stderr, "add %c[%d-%d] \n", res_chain, s2->init_seq_num, s2->term_seq_num); 
            cnlist = contact[scount];

            while (cnlist) { 
              if (cnlist->structure == s2) {
                break;
                }

              cnlist = cnlist->next; 
              }

            if (!cnlist) { 
              cnptr = new PmStructureListPtr; 
              cnptr->structure = s2; 
              cnptr->next = contact[scount];
              contact[scount] = cnptr;
              }
            }
          }
        }

      scount += 1;
      s1 = s1->next; 
      }

    topology[i].contact = contact;

/*
    s1 = topology[i].list;

    for (int j = 0; j < topology[i].num_structures; j++) {
      fprintf (stderr, "<----%c%s----> \n", stype[s1->type], s1->id.c_str()); 
      cnlist = contact[j];

      while (cnlist) { 
        s2 = cnlist->structure;
        fprintf (stderr, "%c%s\n", stype[s2->type], s2->id.c_str()); 
        cnlist = cnlist->next; 
        }

      s1 = s1->next; 
      }
*/
    }

#endif
  }

//*============================================================*
//*==========              createDomain              ==========*
//*============================================================*
// create a domain from two or more domains. 

PmMolecule *
PmMolecule::createDomain (const string name, const vector<string>& domain_names) 
  {
  pm_ErrorWarnReport (PM, "PmMolecule::createDomain: not implemented for lists", "*");
  return NULL;
  }

//*============================================================*
//*==========              createDomain              ==========*
//*============================================================*
// create a domain from a string descriptor.

PmMolecule * 
PmMolecule::createDomain (const string name, const string desc, PmAtomFilter& filter)
  {
  #ifdef dbg_createDomain 
  fprintf (stderr, "\n>>>>>> PmMolecule::createDomain:  name [%s] desc [%s] \n", 
           name.c_str(), desc.c_str());
  #endif

  vector<PmResidue*> rlist, nrlist;
  getResidues(desc, rlist);

  if (!rlist.size()) {
    return NULL;
    }

  // for a domain of a single residue get the      //
  // next residue. this is needed to add an extra  //
  // atom from the next residue.                   //

#ifdef createDomain_add_atom
  if (rlist.size() == 1) {
    char chain_id, str[10];
    int res_id;
    string ndesc;

    if (desc.find('[') != string::npos) {
      sscanf (desc.c_str(), "%c[%d]", &chain_id, &res_id);
      }
    else {
      sscanf (desc.c_str(), "%c%d", &chain_id, &res_id);
      }

    sprintf (str, "%c%d", chain_id, res_id+1);
    ndesc = str;
    getResidues(ndesc, nrlist);
    }
#endif

  PmMolecule *domain = new PmMolecule(name);
  PmResidue *res, *nres;
  PmAtom *atom;
  char *aname;

  //===== add atoms =====//

  for (unsigned int i = 0; i < rlist.size(); i++) {
    res = rlist[i];

    if (i == rlist.size()-1) {
      nres = NULL; 
      }
    else { 
      nres = rlist[i+1];
      }

    if (filter.peptide && nres) { 
      for (int j = 0; j < res->num_atoms; j++) {
        aname = res->atoms[j]->name;

        if ((strlen(aname) == 1) && ((*aname == 'C') || 
            (*aname == 'O'))) {
          domain->addAtom (res->atoms[j]);
          //fprintf (stderr, "add atom [%s] \n", res->atoms[j]->name);
          }
        }

      for (int j = 0; j < nres->num_atoms; j++) {
        aname = nres->atoms[j]->name;

        if ((strlen(aname) == 1) && (*aname == 'N')) {  
          domain->addAtom (nres->atoms[j]);
          //fprintf (stderr, "add atom [%s] \n", nres->atoms[j]->name);
          }
        }
      }

    else if (filter.last_nca_atoms && i == rlist.size()-1) { 
      res->getAtom("N", &atom);
      domain->addAtom (atom);
      res->getAtom("CA", &atom);
      domain->addAtom (atom);
      }

    // res sidechain and c-alpha, nres c-alpha //

    else if (filter.sidechain_ca) { 
      for (int j = 0; j < res->num_atoms; j++) {
        aname = res->atoms[j]->name;

        if (strlen(aname) > 1) { 
          domain->addAtom (res->atoms[j]);
          }
        }

      if (nres) {
        nres->getAtom("CA", &atom);
        domain->addAtom (atom);
        }
      break;
      }

    else if (filter.sidechain_group) { 
      for (int j = 0; j < res->num_atoms; j++) {
        aname = res->atoms[j]->name;

        if (strlen(aname) == 1) {
          if ((*aname != 'O') && (*aname != 'N') && (*aname != 'C')) {
            domain->addAtom (res->atoms[j]);
            //fprintf (stderr, "add atom [%s] \n", res->atoms[j]->name);
            }
          }
        else {
          domain->addAtom (res->atoms[j]);
          //fprintf (stderr, "add atom [%s] \n", res->atoms[j]->name);
          }
        }
      }

    else if (!filter.peptide) {
      for (int j = 0; j < res->num_atoms; j++) {
        if (filter.hasName(res->atoms[j]->name)) {
          domain->addAtom (res->atoms[j]);
          }
        }
      }
    }

  // for a domain with a single residue add  //
  // the next resiudes N atom.               //

#ifdef createDomain_add_atom
  if (nrlist.size()) {
    res = nrlist[0];
    res->getAtom ("N", &atom);

    if (atom) {
      domain->addAtom (atom);
      }
    }
#endif

  // add connect //

  if (connect.size()) {
    int n, dn, id, conn[100];

    for (unsigned int i = 0; i < connect.size();) {
      n = connect[i++];
      id = connect[i];
      domain->getAtom(id, &atom);

      if (!atom) {
        i += n;
        }
      else {
        dn = 0;

        for (int j = 0; j < n; j++, i++) {
          id = connect[i];
          domain->getAtom(id, &atom);

          if (atom) {
            conn[dn++] = id;
            }
          }

        if (dn) {
          domain->addConnect(dn, conn);
          }
        }
      }
    }

  // add modes //

  //fprintf (stderr, "   >>> modes.size[%d] \n", modes.size());

  for (unsigned int m = 0; m < modes.size(); m++) {
    PmMode mode;
    int k, num;
    mode.eigen_vectors = new PmVector3[domain->number_of_atoms];
    mode.eigen_value = modes[m].eigen_value; 
    num = 0;

    for (unsigned int i = 0; i < rlist.size(); i++) {
      res = rlist[i];
  
      for (int j = 0; j < res->num_atoms; j++) {
        k = res->atoms[j]->index;
        mode.eigen_vectors[num] = modes[m].eigen_vectors[k]; 
        num += 1;
        }
      }

    domain->modes.push_back(mode);
    }

  domain->type = type;
  domain->parent = this;
  domain->backbone_atom = backbone_atom;
  pmSystem.addDomain (domain);

  //===== check for atoms =====//

  if (!domain->atoms.size()) { 
    pm_ErrorReport (PM, "domain \"%s\" has no atoms", "*", this->name.c_str());
    return NULL;
    }

  return (domain);
  }

//*============================================================*
//*==========              createDomain              ==========*
//*============================================================*
// create a domain with just a name. 

PmMolecule * 
PmMolecule::createDomain (const string name)
  {
  PmMolecule *domain = new PmMolecule(name);
  domain->parent = this;
  pmSystem.addDomain (domain);
  return (domain);
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get the mass properties for the molecule.

void 
PmMolecule::getMassProps (PmMassProperties& props)
  {
  float mass;

  PmMatrix3x3 inertia;

  float cmx, cmy, cmz, cx, cy, cz;
  float x, y, z, xx, xy, xz, yy, yz, zz;
  float scale = 1.0 / pmSystem.getUnitsMassScale();

  #define ndbg_PmMolecule_getMassProps
  #ifdef dbg_PmMolecule_getMassProps
  fprintf (stderr, ">>>>>> PmMolecule::getMassProps %s  \n", this->name.c_str());
  fprintf (stderr, ">>> mass scale=%g \n", pmSystem.getUnitsMassScale());
  fprintf (stderr, ">>> scale=%g \n", scale); 
  fprintf (stderr, ">>> num atoms=%d \n", number_of_atoms); 
  #endif

  float total_mass = 0.0;
  cmx = cmy = cmz = 0.0;
  cx = cy = cz = 0.0;

  //=====  first compute center of mass =====// 

  for (int i = 0; i < number_of_atoms; i++) {
    x = atoms[i].pos[0];
    y = atoms[i].pos[1];
    z = atoms[i].pos[2];
    mass = scale * atoms[i].getMass();
    total_mass += mass;
    cmx += x*mass;
    cmy += y*mass;
    cmz += z*mass;
    cx += x;
    cy += y;
    cz += z;
    }

  cmx /= total_mass; 
  cmy /= total_mass; 
  cmz /= total_mass;
  cx  /= number_of_atoms;  
  cy  /= number_of_atoms;  
  cz  /= number_of_atoms;

  #ifdef dbg_PmMolecule_getMassProps
  fprintf (stderr, ">>> total mass=%g \n", total_mass); 
  #endif

  //===== compute inertia tensor =====// 

  if (number_of_atoms == 1) {
    float r = atoms[0].getRadius();
    mass = scale * atoms[0].getMass();
    float d = (2.0 / 5.0)*mass*r*r;
    inertia(0,0) = d; 
    inertia(1,1) = d; 
    inertia(2,2) = d; 
    }
  else {
    xx = xy = xz = yy = yz = zz = 0.0;

    for (int i = 0; i < number_of_atoms; i++) {
      x = atoms[i].pos[0] - cmx;
      y = atoms[i].pos[1] - cmy;
      z = atoms[i].pos[2] - cmz;
      mass = scale * atoms[i].getMass();
      xx += mass*x*x;
      xy += mass*x*y;
      xz += mass*x*z;
      yy += mass*y*y;
      yz += mass*y*z;
      zz += mass*z*z;
      }

    inertia(0,0) = yy + zz; 
    inertia(0,1) = -xy; 
    inertia(0,2) = -xz;
    inertia(1,0) = -xy; 
    inertia(1,1) = xx + zz; 
    inertia(1,2) = -yz; 
    inertia(2,0) = -xz; 
    inertia(2,1) = -yz; 
    inertia(2,2) = xx + yy;
    }

  #ifdef dbg_PmMolecule_getMassProps
  fprintf (stderr, ">>> Ixx=%g \n", inertia(0,0)); 
  fprintf (stderr, ">>> Iyy=%g \n", inertia(1,1)); 
  fprintf (stderr, ">>> Izz=%g \n", inertia(2,2)); 
  #endif

  props.mass = total_mass;
  props.com.set(cmx, cmy, cmz);
  props.inertia = inertia;
  }

//*============================================================*
//*==========              addResBond                ==========*
//*============================================================*
// get the atom coordinates for a chain.

void
PmMolecule::addResBond (PmResidue *res, PmResidue *nres, vector<PmBond>& bonds)
  {
  PmAtom *c, *ca, *cb, *cg1, *cg2, *cg, *cd, *cd1, *cz, *cz2, *cz3, *cd2, *ce,
         *ce1, *ce2, *ce3, *ch2, *se;

  PmAtom *n, *ne, *ne1, *ne2, *nh1, *nh2, *nd2, *nd1, *nz;

  PmAtom *o, *og, *og1, *od1, *od2, *oe1, *oe2, *oh;

  PmAtom *sg, *sd;

  PmAtom *n1, *n2, *n3, *n4, *n5, *n6, *n7, *n8, *n9,
         *c1, *c2, *c3, *c4, *c5, *c6, *c7, *c8,
         *c1s, *c2s, *c3s, *c4s, *c5s,
         *o2, *o3, *o4, *o5, *o6, *o1p, *o2p,
         *o2s, *o3s, *o4s, *o5s,
         *p, *np, *o1;

  PmCompoundType ctype = res->ctype;
  int num = bonds.size();

  /*
  fprintf (stderr, ">>>>>> PmMolecule::addResBond:   res id [%d] \n", res->id);
  fprintf (stderr, "   >>> ctype [%d] \n", res->ctype);
  fprintf (stderr, "   >>> mtype [%d] \n", res->mtype);
  */

  if ((res->mtype == PM_MOLECULE_NUCLEOTIDE) || 
      (res->mtype == PM_MOLECULE_NUCLEIC_ACID)) {
    n1 = n2 = n3 = n4 = n5 = n6 = n7 = n8 = n9 = NULL;
    c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = NULL;
    c1s = c2s = c3s = c4s = c5s = NULL;
    o2 = o3 = o4 = o5 = o6 = o1p = o2p = NULL;
    o2s = o3s = o4s = o5s = NULL;
    p = np = NULL;

    res->getAtom ("N1", &n1);
    res->getAtom ("N3", &n3);
    res->getAtom ("C2", &c2);
    res->getAtom ("C4", &c4);
    res->getAtom ("C5", &c5);
    res->getAtom ("C6", &c6);
    addBond (n1, c2, num, bonds);
    addBond (c2, n3, num, bonds); 
    addBond (n3, c4, num, bonds);
    addBond (c4, c5, num, bonds);
    addBond (c5, c6, num, bonds);
    addBond (c6, n1, num, bonds);

    res->getAtom ("P", &p);
    if (nres) nres->getAtom ("P", &np);
    res->getAtom ("O1P", &o1p);
    res->getAtom ("O2P", &o2p);

    if (!o1p) {
      res->getAtom ("OP1", &o1p);
      res->getAtom ("OP2", &o2p);
      }

    res->getAtom ("O5*", &o5s);

    if (o5s) {
      res->getAtom ("O3*", &o3s);
      res->getAtom ("O4*", &o4s);
      res->getAtom ("C3*", &c3s);
      res->getAtom ("C4*", &c4s);
      res->getAtom ("C1*", &c1s);
      res->getAtom ("C2*", &c2s);
      res->getAtom ("C5*", &c5s);
      }
    else {
      res->getAtom ("O5'", &o5s);

      if (o5s) {
        res->getAtom ("O3'", &o3s);
        res->getAtom ("O4'", &o4s);
        res->getAtom ("C3'", &c3s);
        res->getAtom ("C4'", &c4s);
        res->getAtom ("C1'", &c1s);
        res->getAtom ("C2'", &c2s);
        res->getAtom ("C5'", &c5s);
        }
      else {
        res->getAtom ("O3", &o3s);
        res->getAtom ("O4", &o4s);
        res->getAtom ("C3", &c3s);
        res->getAtom ("C4", &c4s);
        res->getAtom ("C1", &c1s);
        res->getAtom ("C2", &c2s);
        res->getAtom ("C5", &c5s);
        }
      }

    res->getAtom ("O2", &o2s);

    if (!o2s) {
      res->getAtom ("O2*", &o2s);

      if (!o2s) {
        res->getAtom ("O2'", &o2s);
        }
      }

    addBond (p, o1p, num, bonds);
    addBond (p, o2p, num, bonds);
    addBond (p, o5s, num, bonds);
    addBond (o5s, c5s, num, bonds);
    addBond (c5s, c4s, num, bonds);
    addBond (c4s, c3s, num, bonds);
    addBond (c4s, o4s, num, bonds);
    addBond (o4s, c1s, num, bonds);
    addBond (c1s, c2s, num, bonds);
    addBond (c2s, c3s, num, bonds);
    addBond (c3s, o3s, num, bonds);
    //addBond (c2s, o2s, num, bonds);
    addBond (np, o3s, num, bonds);

    switch (ctype) {
      case PM_COMPOUND_ADENOSINE:
        res->getAtom ("N6", &n6);
        res->getAtom ("N7", &n7);
        res->getAtom ("N9", &n9);
        res->getAtom ("C8", &c8);
        addBond (c6, n6, num, bonds);
        addBond (c5, n7, num, bonds);
        addBond (n7, c8, num, bonds);
        addBond (c8, n9, num, bonds);
        addBond (n9, c4, num, bonds);
        addBond (c1s, n9, num, bonds);
        //addBond (c1, n9, num, bonds);
        //fprintf (stderr, "   >>> c1 [%x] \n", c1);
      break;

      case PM_COMPOUND_GUANOSINE:
        res->getAtom ("N2", &n2);
        res->getAtom ("N7", &n7);
        res->getAtom ("N9", &n9);
        res->getAtom ("C8", &c8);
        res->getAtom ("O6", &o6);
        addBond (c5, n7, num, bonds);
        addBond (n7, c8, num, bonds);
        addBond (c8, n9, num, bonds);
        addBond (n9, c4, num, bonds);
        addBond (c2, n2, num, bonds);
        addBond (c6, o6, num, bonds);
        addBond (c1s, n9, num, bonds);
      break;

      case PM_COMPOUND_CYTOSINE:
        res->getAtom ("N4", &n4);
        res->getAtom ("O2", &o2);
        addBond (c4, n4, num, bonds);
        addBond (c2, o2, num, bonds);
        addBond (c1s, n1, num, bonds);
      break;

      case PM_COMPOUND_URIDINE:
        res->getAtom ("O4", &o4);
        res->getAtom ("O2", &o2);
        addBond (c2, o2, num, bonds);
        addBond (c4, o4, num, bonds);
        addBond (c1s, n1, num, bonds);
      break;

      case PM_COMPOUND_THYMINE:
        res->getAtom ("C7", &c7);
        res->getAtom ("O2", &o2);
        res->getAtom ("O4", &o4);
        addBond (c2, o2, num, bonds);
        addBond (c4, o4, num, bonds);
        addBond (c7, c5, num, bonds);
        addBond (c1s, n1, num, bonds);
      break;

      default:
      break;
      }
    }

  else if ((res->mtype == PM_MOLECULE_AMINO_ACID) ||
           (res->mtype == PM_MOLECULE_PROTEIN)) {
    res->getAtom("CA", &ca);
    res->getAtom("N",  &n);
    res->getAtom("O",  &o);
    res->getAtom("C",  &c);
    addBond (n, ca, num, bonds);
    addBond (ca, c, num, bonds);
    addBond (c,  o, num, bonds);

    switch (ctype) {
      case PM_COMPOUND_ALANINE:
        //fprintf (stderr, "   >>> alanine   ca [%d]  num [%d] \n", ca->id, num);
        res->getAtom("CB", &cb);
        addBond (ca, cb, num, bonds);
      break;

      case PM_COMPOUND_ARGININE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD", &cd);
        res->getAtom ("NE", &ne);
        res->getAtom ("CZ", &cz);
        res->getAtom ("NH1", &nh1);
        res->getAtom ("NH2", &nh2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd, num, bonds);
        addBond (cd, ne, num, bonds);
        addBond (ne, cz, num, bonds);
        addBond (cz, nh1, num, bonds);
        addBond (cz, nh2, num, bonds);
      break;

      case PM_COMPOUND_ASPARAGINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("OD1", &od1);
        res->getAtom ("ND2", &nd2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, od1, num, bonds);
        addBond (cg, nd2, num, bonds);

      break;

      case PM_COMPOUND_ASPARTIC_ACID:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("OD1", &od1);
        res->getAtom ("OD2", &od2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, od1, num, bonds);
        addBond (cg, od2, num, bonds);
        //fprintf (stderr, "   >>> asp  c[%x]  n[%x] \n", c, n);
      break;

      case PM_COMPOUND_CYSTEINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("SG", &sg);
        addBond (ca, cb, num, bonds);
        addBond (cb, sg, num, bonds);
      break;

      case PM_COMPOUND_GLUTAMIC_ACID:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD", &cd);
        res->getAtom ("OE1", &oe1);
        res->getAtom ("OE2", &oe2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd, num, bonds);
        addBond (cd, oe1, num, bonds);
        addBond (cd, oe2, num, bonds);
      break;

      case PM_COMPOUND_GLUTAMINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD", &cd);
        res->getAtom ("OE1", &oe1);
        res->getAtom ("NE2", &ne2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd, num, bonds);
        addBond (cd, oe1, num, bonds);
        addBond (cd, ne2, num, bonds);
      break;

      case PM_COMPOUND_HISTIDINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("ND1", &nd1);
        res->getAtom ("CD2", &cd2);
        res->getAtom ("CE1", &ce1);
        res->getAtom ("NE2", &ne2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds); 
        addBond (cg, nd1, num, bonds);
        addBond (nd1, ce1, num, bonds);
        addBond (ce1, ne2, num, bonds);
        addBond (ne2, cd2, num, bonds);
        addBond (cd2, cg, num, bonds);
      break;

      case PM_COMPOUND_ISOLEUCINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG1", &cg1);
        res->getAtom ("CG2", &cg2);
        res->getAtom ("CD1", &cd1);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg1, num, bonds);
        addBond (cb, cg2, num, bonds);
        addBond (cg1, cd1, num, bonds);
      break;

      case PM_COMPOUND_LEUCINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD1", &cd1);
        res->getAtom ("CD2", &cd2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd1, num, bonds);
        addBond (cg, cd2, num, bonds);
      break;

      case PM_COMPOUND_LYSINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD", &cd);
        res->getAtom ("CE", &ce);
        res->getAtom ("NZ", &nz);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd, num, bonds);
        addBond (cd, ce, num, bonds);
        addBond (ce, nz, num, bonds);
      break;

      case PM_COMPOUND_METHIONINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("SD", &sd);
        res->getAtom ("CE", &ce);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, sd, num, bonds);
        addBond (sd, ce, num, bonds);
      break;

      case PM_COMPOUND_SELMETHIONINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("SD", &sd);
        res->getAtom ("SE", &se);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, sd, num, bonds);
        addBond (sd, se, num, bonds);
      break;

      case PM_COMPOUND_PHENYLALANINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD1", &cd1);
        res->getAtom ("CD2", &cd2);
        res->getAtom ("CE1", &ce1);
        res->getAtom ("CE2", &ce2);
        res->getAtom ("CZ", &cz);
        res->getAtom ("O1", &o1);
        res->getAtom ("O2", &o2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd1, num, bonds);
        addBond (cd1, ce1, num, bonds);
        addBond (ce1, cz, num, bonds);
        addBond (cz, ce2, num, bonds);
        addBond (ce2, cd2, num, bonds);
        addBond (cg, cd2, num, bonds);

        if (o1) {
          addBond (c, o1, num, bonds);
          }

        if (o2) {
          addBond (c, o2, num, bonds);
          }
      break;

      case PM_COMPOUND_PROLINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD", &cd);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cd, cg, num, bonds);
        addBond (cd, n, num, bonds);
      break;

      case PM_COMPOUND_SERINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("OG", &og);
        addBond (ca, cb, num, bonds);
        addBond (cb, og, num, bonds);
      break;

      case PM_COMPOUND_THREONINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("OG1", &og1);
        res->getAtom ("CG2", &cg2);
        addBond (ca, cb, num, bonds);
        addBond (cb, og1, num, bonds);
        addBond (cb, cg2, num, bonds);
      break;

      case PM_COMPOUND_TRYPTOPHAN:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD1", &cd1);
        res->getAtom ("CD2", &cd2);
        res->getAtom ("NE1", &ne1);
        res->getAtom ("CE2", &ce2);
        res->getAtom ("CE3", &ce3);
        res->getAtom ("CZ2", &cz2);
        res->getAtom ("CZ3", &cz3);
        res->getAtom ("CH2", &ch2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd1,num, bonds);
        addBond (cd1, ne1, num, bonds);
        addBond (ne1, ce2, num, bonds);
        addBond (ce2, cz2, num, bonds);
        addBond (cz2, ch2, num, bonds);
        addBond (ch2, cz3, num, bonds);
        addBond (cz3, ce3, num, bonds);
        addBond (ce3, cd2, num, bonds);
        addBond (cd2, ce2, num, bonds);
        addBond (cd2, cg, num, bonds);
      break;

      case PM_COMPOUND_TYROSINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG", &cg);
        res->getAtom ("CD1", &cd1);
        res->getAtom ("CD2", &cd2);
        res->getAtom ("CE1", &ce1);
        res->getAtom ("CE2", &ce2);
        res->getAtom ("CZ", &cz);
        res->getAtom ("OH", &oh);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg, num, bonds);
        addBond (cg, cd1, num, bonds);
        addBond (cd1, ce1, num, bonds);
        addBond (ce1, cz, num, bonds);
        addBond (cz, oh, num, bonds);
        addBond (cz, ce2, num, bonds);
        addBond (ce2, cd2, num, bonds);
        addBond (cd2, cg, num, bonds);
      break;

      case PM_COMPOUND_VALINE:
        res->getAtom ("CB", &cb);
        res->getAtom ("CG1", &cg1);
        res->getAtom ("CG2", &cg2);
        res->getAtom ("O1", &o1);
        res->getAtom ("O2", &o2);
        addBond (ca, cb, num, bonds);
        addBond (cb, cg1, num, bonds);
        addBond (cb, cg2, num, bonds);

        if (o1) {
          addBond (c, o1, num, bonds);
          }

        if (o2) {
          addBond (c, o2, num, bonds);
          }
      break;

      default:
      break;
      }

    if (nres) {
      nres->getAtom ("N", &n);
      addBond (c, n, num, bonds);
      }
    }
  }

//*============================================================*
//*==========              addResPlane               ==========*
//*============================================================*
// get the atom plane for a chain.

void
PmMolecule::addResPlane(PmResidue *res, PmResidue *nres, vector<PmMolPlane>& planes)
  {
  //fprintf(stderr, "\n>>>>>> PmMolecule::addResPlane  \n");
  PmCompoundType ctype = res->ctype;

  if ((res->mtype == PM_MOLECULE_NUCLEOTIDE) ||
      (res->mtype == PM_MOLECULE_NUCLEIC_ACID)) {

    switch (ctype) {
      case PM_COMPOUND_ADENOSINE:
      break;

      default:
      break;
      }
    }

  else if ((res->mtype == PM_MOLECULE_AMINO_ACID) ||
           (res->mtype == PM_MOLECULE_PROTEIN)) {
    PmAtom *c, *ca, *n, *o, *nca;
    PmMolPlane plane;
    PmVector3 u, v, w, h, ru, rw, norm, nca_pos, ca_pos, n_pos;
    PmMatrix3x3 rmat;

    res->getAtom("CA", &ca);
    res->getAtom("N",  &n);
    res->getAtom("C",  &c);
    res->getAtom("O",  &o);

    /*
    fprintf(stderr, "   >>> nres %x \n", nres);
    fprintf(stderr, "   >>> ca %x \n", ca);
    fprintf(stderr, "   >>> c  %x \n", c);
    fprintf(stderr, "   >>> o  %x \n", o);
    fprintf(stderr, "   >>> n  %x \n", n);
    */

    // if no ca then this is a special peptide residue //

    if (ca == NULL) {
      if (!nres) {
        return;
        }

      nres->getAtom("N", &n);
      u = o->pos - c->pos;
      u.normalize();
      v = n->pos - c->pos;
      v.normalize();
      norm = u.cross(v);
      norm.normalize();

      // compute ca atom position //

      pm_MathRotationAroundAxis (121.0, norm, rmat);
      ru = rmat*u;
      ca_pos = c->pos + 0.153*ru;

      // compute next ca atom position //
 
      w = c->pos - n->pos;
      w.normalize();
      pm_MathRotationAroundAxis (-123.0, norm, rmat);
      rw = rmat*w;
      nca_pos = n->pos + 0.147*rw;

      plane.vertices.push_back(ca_pos);
      plane.vertices.push_back(o->pos);
      plane.vertices.push_back(nca_pos);
      v = ca_pos - o->pos;
      h = nca_pos + v;
      plane.vertices.push_back(h);
      planes.push_back(plane);
      return;
      }


    // if no oxygen then this is a special sidechain residue //

    if (o == NULL) {
      if (!nres) {
        return;
        }

      nres->getAtom("O", &o);
      u = ca->pos - c->pos;
      u.normalize();
      v = o->pos - c->pos;
      v.normalize();
      norm = u.cross(v);
      norm.normalize();

      pm_MathRotationAroundAxis (114.0, norm, rmat);
      ru = rmat*u;
      n_pos = c->pos + 0.132*ru;

      w = c->pos - n_pos;
      w.normalize();
      pm_MathRotationAroundAxis (-123.0, norm, rmat);
      rw = rmat*w;
      nca_pos = n_pos + 0.147*rw;

      plane.vertices.push_back(ca->pos);
      plane.vertices.push_back(o->pos);
      plane.vertices.push_back(nca_pos);
      v = ca->pos - o->pos;
      h = nca_pos + v;
      plane.vertices.push_back(h);
      planes.push_back(plane);
      return;
      }

    plane.atoms.push_back(ca);
    plane.atoms.push_back(c);
    plane.atoms.push_back(o);
    plane.atoms.push_back(n);

    if (nres) {
      nres->getAtom("CA", &nca);
      plane.vertices.push_back(ca->pos);
      plane.vertices.push_back(o->pos);
      plane.vertices.push_back(nca->pos);
      v = ca->pos - o->pos;
      h = nca->pos + v;
      plane.vertices.push_back(h);
      planes.push_back(plane);
      }
    }
  }

//*============================================================*
//*==========              getType                   ==========*
//*============================================================*
// get the type of the molecule.

PmMoleculeType 
PmMolecule::getType()
  {
  if (type != PM_MOLECULE_UNKNOWN) {
    return type;
    }

  if ((atoms[0].mtype == PM_MOLECULE_AMINO_ACID) ||
      (atoms[0].mtype == PM_MOLECULE_PROTEIN)) { 
    type = PM_MOLECULE_PROTEIN;
    }

  else if (atoms[0].mtype == PM_MOLECULE_NUCLEOTIDE) { 
    type = PM_MOLECULE_NUCLEOTIDE;
    }

  else if (atoms[0].mtype == PM_MOLECULE_NUCLEIC_ACID) { 
    type = PM_MOLECULE_NUCLEIC_ACID;
    }

  else {
    type = atoms[0].mtype;
    }

  return type;
  }

//*============================================================*
//*==========                  setType               ==========*
//*============================================================*
// set the molecule type.                 

void 
PmMolecule::setType (PmMoleculeType type) 
  { 
  //fprintf (stderr, "\n>>> PmMolecule::set %d \n", type);
  this->type = type; 

  if ((type == PM_MOLECULE_NUCLEIC_ACID) || (type == PM_MOLECULE_NUCLEOTIDE)) {
    backbone_atom = "P";
    }
  else { 
    backbone_atom = "CA";
    }

  for (int i = 0; i < number_of_atoms; i++) {
    atoms[i].mtype = type; 
    }
  }

////////////////////////////////////////////////////////////////
//                    g r a p h i c s                        //
//////////////////////////////////////////////////////////////

void 
PmMolecule::convMoleculeRenderType(const string str, PmMoleculeRenderType& type)
  {
  type = PM_MOLECULE_RENDER_UNKNOWN;

  if (str == "point") {
    type = PM_MOLECULE_RENDER_POINT;
    }
  else if (str == "line") {
    type = PM_MOLECULE_RENDER_LINE;
    }
  else if (str == "solid") {
    type = PM_MOLECULE_RENDER_SOLID;
    }
  }

void
PmMolecule::convMoleculeDisplayType(const string str, PmMoleculeDisplayType& type)
  {
  type = PM_MOLECULE_DISPLAY_UNKNOWN;

  if (str == "point") {
    type = PM_MOLECULE_DISPLAY_POINT;
    }
  else if (str == "cross") {
    type = PM_MOLECULE_DISPLAY_CROSS;
    }
  else if (str == "sphere") {
    type = PM_MOLECULE_DISPLAY_SPHERE;
    }
  }

//*============================================================*
//*==========          addGraphicsGeometry           ==========*
//*============================================================*
// add graphics geometry object.

void
PmMolecule::addGraphicsGeometry(PmGraphicsGeometry *geom) {
  graphics_geometry.push_back(geom);
  }

//*============================================================*
//*==========          getGraphicsGeometry           ==========*
//*============================================================*
// get a graphics geometry object.

void
PmMolecule::getGraphicsGeometry(string name, PmGraphicsGeometry **geom)
  {
  *geom = NULL;

  for (unsigned int i = 0; i < graphics_geometry.size(); i++) {
    if (name == graphics_geometry[i]->name) {
      *geom = graphics_geometry[i];
      }
    }
  }

//*============================================================*
//*==========                display                 ==========*
//*============================================================*
// display the geometry already created for the molecule.

void
PmMolecule::display(const bool show)
  {
  fprintf (stderr, "\n>>> PmMolecule::display \n");
  }

//*============================================================*
//*==========          getAtomsGeometry              ==========*
//*============================================================*
// get a graphics atoms geometry object.

void
PmMolecule::getAtomsGeometry (string name, PmGraphicsAtoms **ageom)
  {
  *ageom = NULL; 
  PmGraphicsGeometry *geom = NULL; 
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    return;
    }

  *ageom = dynamic_cast<PmGraphicsAtoms*>(geom);
  }

//*============================================================*
//*==========              displayAtoms              ==========*
//*============================================================*
// display atoms.                     

void
PmMolecule::displayAtoms(const string gname, PmAtomFilter& filter, const bool show)
  {

  if (!pmSystem.useGraphics()) return;
  if (!num_chains) buildChains(); 
  
  int num_coords; 
  PmVector3 *coords, *colors = NULL; 
  int num_rads;
  float *rads;
  string geom_name;
  PmGraphicsAtoms *geom;

  #define ndbg_displayAtoms   
  #ifdef dbg_displayAtoms
  fprintf (stderr, "\n>>>>>> PmMolecule::displayAtoms: \n");
  fprintf (stderr, "   >>> num chains [%d] \n", num_chains);
  fprintf (stderr, "   >>> datts.color_type [%d] \n", datts.color_type);
  fprintf (stderr, "   >>> gname = %s \n", gname.c_str());
  #endif

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    string name;
    #ifdef dbg_displayAtoms
    fprintf (stderr, "---- chain %c ---- \n", chain->id);
    #endif
    buildGeomName (chain->id, pmGraphicsGeomAtoms, gname, geom_name);
    getAtomsGeometry (geom_name, &geom);

    // if there is no graphics geometry then create it //

    if (!geom) {

      // display only a selected atoms //

      if (filter.names.size()) {
        vector<PmVector3> vcoords;
        vector<float> vrads;
        vector<PmAtom> catoms;
        vector<PmVector3> vcolors;
        PmAtom atom;
        PmVector3 acolor;
        float r;
        getChainAtoms(chain->id, catoms);

        for (unsigned int j = 0; j < catoms.size(); j++) {
          atom = catoms[j];

          if (filter.hasName(atom.name)) {
            vcoords.push_back(atom.pos);
            r = atom.getRadius();
            vrads.push_back(r);
            atom.getColor (acolor);

            if (datts.color_type == PM_ATOM_COLOR_ELEMENT) {
              vcolors.push_back(acolor);
              }
            }
          }

        num_coords = vcoords.size(); 

        if (num_coords) {
          coords = new PmVector3[num_coords];
          rads = new float[num_coords];

          if (datts.color_type == PM_ATOM_COLOR_ELEMENT) {
            colors = new PmVector3[num_coords];

            for (int j = 0; j < num_coords; j++) {
              colors[j] = vcolors[j]; 
              }
            }

          for (int j = 0; j < num_coords; j++) {
            coords[j] = vcoords[j]; 
            rads[j] = vrads[j]; 
            }
          }
        }

      // display all atoms //

      else {
        getChainAtomCoords(chain->id, name, &num_coords, &coords);
        getChainAtomRadii(chain->id, name, &num_rads, &rads);
        }

      if (num_coords) {
        #ifdef dbg_displayAtoms
        fprintf (stderr, ">>>> geom [%s] \n", geom_name.c_str());
        #endif

        if (datts.color_type == PM_ATOM_COLOR_ELEMENT) {
          if (!colors) {
            getChainAtomColors (chain->id, &colors);
            }
          }

        geom = new PmGraphicsAtoms (geom_name, num_coords, coords, colors);
        addGraphicsGeometry(geom);
        geom->setRadius(num_rads, rads);
        }
      }

    PmGraphicsAttributes atts;
    atts.setColor (this->datts.color);
    atts.setShadingType(PM_GEOMETRY_SHADING_COLOR);
    geom->setDisplay(this->datts.type);
    geom->setRender(this->datts.render);
    geom->setAttributes(atts);
    geom->display();
    }
  }

//*============================================================*
//*==========              buildGeomName             ==========*
//*============================================================*

void
PmMolecule::buildGeomName (char chain, const string geom, const string gname, 
                           string& gr_name) 
  {

  string mol_name; 

  if (parent)  {
    mol_name = "domain";
    }
  else {
    mol_name = "molecule";
    }

  if (gname == "") {
    if (chain != '\0') {
      gr_name = mol_name + '[' + name + ':' + chain + ']' + geom;
      }
    else {
      gr_name = mol_name + '[' + name + ']' + geom;
      }
    }
  else {
    gr_name = mol_name + '[' + name + '_' + gname + ':' + chain + ']' + geom;
    }
  }

//*============================================================*
//*==========              displayBackbone           ==========*
//*============================================================*
// display backbone atoms.                     

void
PmMolecule::displayBackbone(const string gname, const bool show)
  {
  if (!pmSystem.useGraphics()) return;
  if (!num_chains) buildChains(); 
  
  int num_coords; 
  PmVector3 *coords, pt1, pt2, v; 
  PmGraphicsBackbone *geom;
  PmGraphicsBackboneTube *tgeom;
  string geom_name;
  int n, n1, n2, num_breaks; 
  float tol = 0.5, d;
  vector<int> ends;
  vector<PmGraphicsBackboneTube*> tgeom_list;

  #define dbg_displayBackbone
  #ifdef dbg_displayBackbone
  fprintf (stderr, "\n>>>>>> PmMolecule::displayBackbone: \n");
  fprintf (stderr, "   >>> num chains [%d] \n", num_chains);
  #endif

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    #ifdef dbg_displayBackbone
    fprintf (stderr, "---- chain %c ---- \n", chain->id);
    #endif

    buildGeomName (chain->id, pmGraphicsGeomBackbone, gname, geom_name);
    getBackboneGeometry (geom_name, &geom);

    if (!geom) {
      getChainAtomCoords (chain->id, backbone_atom, &num_coords, &coords);

      if (num_coords) {
        #ifdef dbg_displayBackbone
        fprintf (stderr, ">>> geom [%s] \n", geom_name.c_str());
        #endif

        // test for non-contiguous backbone //

        if ((type == PM_MOLECULE_NUCLEOTIDE) ||
            (type == PM_MOLECULE_NUCLEIC_ACID)) {
          tol = 1.0;
          }

        ends.clear();
        num_breaks = 0; 
        n1 = 0;

        for (int i = 0; i < num_coords-1; i++) {
          v = coords[i] - coords[i+1];
          d = v.length();
 
          if (d > tol) {
            ends.push_back(n1);
            ends.push_back(i);
            n1 = i+1;
            num_breaks += 1; 
            }
          }

        ends.push_back(n1);
        ends.push_back(num_coords-1);
        num_breaks += 1; 

        // if the backbone is broken then display //
        // using connectivity.                    //

        if (num_breaks > 1) {
          #ifdef dbg_displayBackbone
          fprintf (stderr, ">>> num breaks = %d \n", num_breaks); 
          #endif

          if (this->datts.tube) {
            fprintf (stderr, "[PmMolecule::displayBackbone] tube geom \n"); 
            tgeom = new PmGraphicsBackboneTube (geom_name, num_breaks, ends,
                                               num_coords, coords, 
                                               this->datts.tube_res, 
                                               this->datts.tube_width);
            addGraphicsGeometry (tgeom);
            }

          else {
            PmConn *conn = new PmConn(num_breaks + num_coords);
            n = 0;
            #ifdef dbg_displayBackbone
            fprintf (stderr, ">>> num breaks = %d \n", num_breaks);
            #endif

            for (int i = 0; i < num_breaks; i++) {
              n1 = ends[2*i];
              n2 = ends[2*i+1];
              (*conn)[n++] = n2 - n1 + 1;

              for (int j = n1; j <= n2; j++) {
                (*conn)[n++] = j;
                }
              }

            geom = new PmGraphicsBackbone (geom_name, num_breaks, conn, num_coords, 
                                           coords);
            addGraphicsGeometry (geom);
            }
          }
        else {
          if (this->datts.tube) {
            fprintf (stderr, "[PmMolecule::displayBackbone] tube 2 geom \n"); 
            tgeom = new PmGraphicsBackboneTube (geom_name, num_coords, coords, 
                                                this->datts.tube_res, 
                                                this->datts.tube_width);
            addGraphicsGeometry (tgeom);
            }
          else {
            geom = new PmGraphicsBackbone (geom_name, num_coords, coords);
            addGraphicsGeometry (geom);
            }
          }
        }
      else {
        fprintf (stderr, "**** WARNING[displayBackbone] no backbone %s atoms. \n", 
                 backbone_atom.c_str());
        return;
        }
      }

    PmGraphicsAttributes atts;
    atts.setLineWidth (this->datts.line_width);
    atts.setColor (this->datts.color);
    atts.setVisible (show);

    if (this->datts.tube) {
      tgeom->setAttributes(atts);
      tgeom->display();
      }
    else {
      geom->setAttributes(atts);
      geom->display();
      }
    }
  }

//*============================================================*
//*==========          getBackboneGeometry           ==========*
//*============================================================*
// get a graphics backbone geometry object.

void
PmMolecule::getBackboneGeometry (string name, PmGraphicsBackbone **bgeom)
  {
  *bgeom = NULL;
  PmGraphicsGeometry *geom = NULL;
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    return;
    }

  *bgeom = dynamic_cast<PmGraphicsBackbone*>(geom);
  }

//*============================================================*
//*==========          getBondsGeometry              ==========*
//*============================================================*
// get a graphics bonds geometry object.

void
PmMolecule::getBondsGeometry (string name, PmGraphicsBonds **bgeom)
  {
  *bgeom = NULL;
  PmGraphicsGeometry *geom = NULL;
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    return;
    }

  *bgeom = dynamic_cast<PmGraphicsBonds*>(geom);
  }

//*============================================================*
//*==========              displayBonds              ==========*
//*============================================================*
// display atom bonds.                     

void
PmMolecule::displayBonds (const string gname, const bool show)
  {
  if (!pmSystem.useGraphics()) return;
  if (!num_chains) buildChains(); 
  
  int num_coords; 
  PmVector3 *coords, *colors = NULL; 
  PmGraphicsBonds *geom;
  string geom_name, atom_name;
  vector<PmBond> bonds;

  #define ndbg_displayBonds
  #ifdef dbg_displayBonds
  fprintf (stderr, "\n>>>>>> PmMolecule::displayBonds: \n");
  fprintf (stderr, "   >>> num chains [%d] \n", num_chains);
  #endif

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    #ifdef dbg_displayBonds
    fprintf (stderr, "---- chain %c ---- \n", chain->id);
    #endif

    buildGeomName (chain->id, pmGraphicsGeomBonds, gname, geom_name);
    getBondsGeometry (geom_name, &geom);

    if (!geom) {
      getChainAtomBonds (chain->id, bonds);
      getChainAtomCoords (chain->id, atom_name, &num_coords, &coords);

      if (num_coords) {
        int num_bonds = bonds.size();
        /*
        if (adj_res && (adjacent[0] == chain->id)) {
          extra_bond = 1;
          }
        */
        #ifdef dbg_displayBonds
        fprintf (stderr, ">>>> geom [%s] \n", geom_name.c_str());
        fprintf (stderr, "   >>> num coords [%d] \n", num_coords);
        fprintf (stderr, ">>>> num bonds [%d] \n", num_bonds);
        #endif

        // need to compute the smallest index to offset the atoms  //
        // indexes that index into the atoms for the entire domain. //

        PmConn *conn = NULL; 

        if (num_bonds) {
          int min_index = bonds[i].atom1->index;
          PmAtom *a1, *a2;

          for (int i = 0; i < num_bonds; i++) {
            a1 = bonds[i].atom1;
            a2 = bonds[i].atom2;
            if (a1->index < min_index) min_index = a1->index;
            if (a2->index < min_index) min_index = a2->index;
            }

          // add the connectivity //

          conn = new PmConn(3*num_bonds);
          int n = 0;

          for (int i = 0; i < num_bonds; i++) {
            a1 = bonds[i].atom1;
            a2 = bonds[i].atom2;
            //fprintf (stderr, ">>>> bond: %s  %s \n", a1->name, a2->name);
            (*conn)[n++] = 2;
            (*conn)[n++] = a1->index - min_index;
            (*conn)[n++] = a2->index - min_index;
            /*
            fprintf (stderr, ">>>> bond: %s(%d)  %s(%d) \n", a1->name, a1->index,
                    a2->name, a2->index);
            */
            }
          }

        if (datts.bond_color_type == PM_ATOM_COLOR_ELEMENT) {
          getChainAtomColors (chain->id, &colors);
          }

        if (datts.color_type == PM_ATOM_COLOR_ELEMENT) {
          getChainAtomColors (chain->id, &colors);
          }

        geom = new PmGraphicsBonds (geom_name, num_coords, coords, num_bonds, conn,
                                    colors);
        addGraphicsGeometry (geom);
        }
      }

    if (!geom) {
      return;
      }

    if ((datts.bond_color_type == PM_ATOM_COLOR_ELEMENT) && !geom->hasColors()) {
      getChainAtomColors (chain->id, &colors);
      geom->setColors (num_coords, colors);
      }

    /*
    fprintf (stderr, ">>>> bond color = %f %f %f \n", this->datts.bond_color[0],
             this->datts.bond_color[1], this->datts.bond_color[2]); 

    fprintf (stderr, ">>>> atom color = %f %f %f \n", this->datts.atom_color[0],
             this->datts.atom_color[1], this->datts.atom_color[2]); 
    */

    geom->setAtomColor(this->datts.atom_color);
    geom->setBondColor(this->datts.bond_color);
    geom->setBondAtoms(this->datts.bond_atoms);
    geom->setAtomRenderType(datts.bond_atom_render);

    PmGraphicsAttributes atts;
    atts.setLineWidth (this->datts.line_width);
    atts.setScale(this->datts.marker_size);
    geom->setAttributes(atts);
    geom->display();
    }
  }

//*============================================================*
//*==========         displayHydrogenBonds           ==========*
//*============================================================*
// display hydrogen bonds.                     

void
PmMolecule::displayHydrogenBonds(const PmVector3& color, const float width, 
                                 const bool show)
  {
  if (!pmSystem.useGraphics()) return;
  if (!num_chains) buildChains(); 
  
  int num_verts; 
  PmVector3 *verts;
  PmGraphicsGeometry *geom;
  string gname, geom_name, atom_name;
  int num_bonds;
  vector<PmBond> bonds;
  vector<PmResidue*> rlist1, rlist2;
  PmGraphicsAttributes atts;

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain1 = &chains[i];
    getChainResidues (chain1->id, rlist1);

    for (int j = i; j < num_chains; j++) {
      PmChain *chain2 = &chains[j];
      getChainResidues (chain2->id, rlist2);
      getHydrogenBonds (chain1->id, rlist1, chain2->id, rlist2, bonds);
      num_bonds = bonds.size();

      if (num_bonds != 0) {
        gname = chain2->id;
        buildGeomName (chain1->id, pmGraphicsGeomHbonds, gname, geom_name);
        getGraphicsGeometry(geom_name, &geom);

        if (!geom) {
          verts = new PmVector3[2*num_bonds];
          num_verts = 0;

          for (int k = 0; k < num_bonds; k++) {
            verts[num_verts++] = bonds[k].atom1->pos;
            verts[num_verts++] = bonds[k].atom2->pos;
            }

          geom = new PmGraphicsLine(geom_name, num_verts, verts);
          addGraphicsGeometry (geom);
          atts.setDisjoint(true);
          atts.setColor(color); 
          atts.setLineWidth(width); 
          geom->setAttributes(atts);
          }

        geom->display();
        }
      }
    }
  }

//*============================================================*
//*==========              getHydrogenBonds          ==========*
//*============================================================*
// get the hydrogen bonds btween two residue lists.

void
PmMolecule::getHydrogenBonds (const char cid1, vector<PmResidue*> rlist1, 
                              const char cid2, vector<PmResidue*> rlist2, 
                              vector<PmBond>& bonds)
  {
  /*
  fprintf (stderr, "\n>>>>>>> PmMolecule::getHydrogenBonds \n");
  fprintf (stderr, "   >>>> cid1 %c   cid2 %c \n", cid1, cid2);
  */
  PmAtom *n, *o;
  vector<PmAtom*> oxygen1, nitrogen1, oxygen2, nitrogen2;
  int num_hbonds;
  PmVector3 v, opos, npos; 
  float d;
  float cutoff = 0.315;
  PmBond bond;

  num_hbonds = 0;
  bonds.clear();

  // if the same chain then first add backbone bonds //

  if (cid1 == cid2) {
    for (unsigned int i = 0; i < rlist1.size() - 4; i++) {
      rlist1[i]->getAtom("N", &n); 

      if (n) {
        rlist1[i+4]->getAtom("O", &o); 

        if (o) {
          v = o->pos - n->pos;
          d = v.length();

          if (d <= cutoff) {
            num_hbonds += 1;
            bond.atom1 = o;
            bond.atom2 = n;
            bonds.push_back(bond);
            }
          }
        }

      rlist1[i]->getAtom("O", &o); 

      if (o) {
        rlist1[i+4]->getAtom("N", &n); 

        if (n) {
          v = o->pos - n->pos;
          d = v.length();

          if (d <= cutoff) {
            num_hbonds += 1;
            bond.atom1 = o;
            bond.atom2 = n;
            bonds.push_back(bond);
            }
          }
        }
      }
    }

  //fprintf (stderr, ">>>>>>> num hbonds = %d \n", num_hbonds);
  }

//*============================================================*
//*==========              displayPca                ==========*
//*============================================================*
// display pca axes.   

void
PmMolecule::displayPca(PmPcaResults& pca, bool show)
  {
  if (!pmSystem.useGraphics()) return;

  PmGraphicsAxes *geom;
  string geom_name;

  buildGeomName ('\0', "axes", "", geom_name);
  getAxesGeometry (geom_name, &geom);

  if (!geom) {
    PmVector3 com, axis1, axis2, axis3;
    com = pca.com;
    axis1 = pca.s1*pca.axis1;
    axis2 = pca.s2*pca.axis2;
    axis3 = pca.s3*pca.axis3;
    geom = new PmGraphicsAxes(geom_name, com, axis1, axis2, axis3);
    addGraphicsGeometry (geom);
    }

  geom->display();
  }

//*============================================================*
//*==========          getAxesGeometry               ==========*
//*============================================================*
// get a graphics axes geometry object.

void
PmMolecule::getAxesGeometry (string str, PmGraphicsAxes **ageom)
  {
  *ageom = NULL;
  PmGraphicsGeometry *geom = NULL;
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    return;
    }

  *ageom = dynamic_cast<PmGraphicsAxes*>(geom);
  }

//*============================================================*
//*==========          getPlanesGeometry             ==========*
//*============================================================*
// get a graphics planes geometry object.

void
PmMolecule::getPlanesGeometry (string str, PmGraphicsPlanes **pgeom)
  {
  *pgeom = NULL;
  PmGraphicsGeometry *geom = NULL;
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    return;
    }

  *pgeom = dynamic_cast<PmGraphicsPlanes*>(geom);
  }

//*============================================================*
//*==========              displayBackbonePlanes     ==========*
//*============================================================*
// display backbone planes.                

void
PmMolecule::displayBackbonePlanes(const string gname, PmResidue *aux_res, 
                                  PmVector3& color, const bool show)
  {
  if (!pmSystem.useGraphics()) return;
  if (!num_chains) buildChains(); 
  
  int num_coords; 
  PmVector3 *coords, *colors = NULL; 
  int num_planes; 
  PmConn *conn;
  PmGraphicsPlanes *geom;
  string geom_name, atom_name;

  #define ndbg_displayBackbonePlanes
  #ifdef dbg_displayBackbonePlanes
  fprintf (stderr, "\n>>>>>> PmMolecule::displayBackbonePlanes: \n");
  fprintf (stderr, "   >>> num chains [%d] \n", num_chains);
  #endif

  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    #ifdef dbg_displayBackbonePlanes
    fprintf (stderr, "---- chain %c ---- \n", chain->id);
    #endif

    buildGeomName (chain->id, pmGraphicsGeomBackbonePlanes, gname, geom_name);
    getPlanesGeometry (geom_name, &geom);

    if (!geom) {
      getChainPlaneGeometry(chain->id, aux_res, num_planes, &conn, num_coords, &coords); 

      if (num_planes == 0) {
        return;
        }

      geom = new PmGraphicsPlanes(geom_name, num_coords, coords, num_planes, conn,
                                  colors);
      addGraphicsGeometry (geom);
      }

    /*
    geom->setAtomColor(this->datts.atom_color);
    geom->setPlaneColor(this->datts.bond_color);
    geom->setBondAtoms(this->datts.bond_atoms);
    geom->setAtomRenderType(datts.bond_atom_render);

    PmGraphicsAttributes atts;
    atts.setLineWidth (this->datts.line_width);
    geom->setAttributes(atts);
    */
    PmGraphicsAttributes atts;
    atts.setColor(color);
    geom->setAttributes(atts);
    geom->display();
    }
  }

//*============================================================*
//*==========              displaySurface            ==========*
//*============================================================*
// display a molecule surface.

void
PmMolecule::displaySurface(const string gname, const bool show)
  {

  if (!pmSystem.useGraphics()) return;

  if (!surface) {
    return;
    }

  surface->display(show);
  }

//*============================================================*
//*==========          getSurfaceGeometry            ==========*
//*============================================================*
// get a graphics surface geometry object.

void
PmMolecule::getSurfaceGeometry(string name, PmGraphicsSurface **sgeom)
  {
  *sgeom = NULL;
  PmGraphicsGeometry *geom = NULL;
  getGraphicsGeometry (name, &geom);

  if (!geom) {
    return;
    }

  *sgeom = dynamic_cast<PmGraphicsSurface*>(geom);
  }

//*============================================================*
//*==========              getXform                  ==========*
//*============================================================*

void
PmMolecule::getXform (PmXform& xform) {
  xform = this->xform;
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*

void 
PmMolecule::setXform (PmXform& xform) 
  {
  this->xform = xform;
  /*
  fprintf (stderr, ">>>>>> PmMolecule::setXform  [%s] \n", name.c_str());
  fprintf (stderr, ">>> xform.translation=%f %f %f\n", xform.translation[0],
           xform.translation[1], xform.translation[2]); 
  */
  //fprintf (stderr, "     num backbone geoms    [%d] \n", backbone_geom.size()); 

  for (unsigned int i = 0; i < graphics_geometry.size(); i++) {
    PmGraphicsGeometry *geom = graphics_geometry[i];
    geom->setXform(xform);
    geom->display();
    }

  if (surface) {
    surface->setXform(xform);
    }
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              addBond                   ==========*
//*============================================================*

void
PmMolecule::addBond(PmAtom *a1, PmAtom *a2, int& n, vector<PmBond>& bonds)
  {
  //fprintf (stderr, ">>>> PmMolecule::addBond  a1:%x a2:%x ] \n", a1, a2);
  if (!a1 || !a2) return;
  PmBond bond;
  bond.atom1 = a1;
  bond.atom2 = a2;
  bonds.push_back (bond);
  n = n + 1;
  }

//*============================================================*
//*==========              buildChains               ==========*
//*============================================================*
// build the chains for the molecule.

void
PmMolecule::buildChains()
  {
  if (num_chains != 0) return;
  int n, starti, alist[pmMaxAtomsPerResidue];

  #define ndbg_buildChains
  #ifdef dbg_buildChains
  fprintf (stderr, "\n>>>>>> PmMolecule::buildChains:  name [%s] \n", name.c_str());
  #endif

  //===== count chains and residues =====//

  int res_id = atoms[0].seq;
  num_residues = 1;
  char chain_id = atoms[0].chain;
  num_chains = 1;

  for (int i = 1; i < number_of_atoms; i++) {
    if (atoms[i].seq != res_id) {
      res_id = atoms[i].seq;
      num_residues += 1;
      }

    if (atoms[i].chain != chain_id) {
      chain_id = atoms[i].chain;
      num_chains += 1;
      }
    }

  #ifdef dbg_buildChains
  fprintf (stderr, "  num chains   [%d] \n", num_chains);
  fprintf (stderr, "  num residues [%d] \n", num_residues);
  #endif

  if (num_residues == 0) {
    pm_ErrorReport (PM, "domain \"%s\" has no residues", "*", this->name.c_str());

    for (int i = 0; i < number_of_atoms; i++) {
      fprintf (stderr, ">>> atom  chain=%c seq=%d \n", atoms[i].chain, atoms[i].seq);
      }
    return;
    }

  //===== create residue data =====//

  residues.resize (num_residues);

  alist[0] = 0;
  n = 1;
  num_residues = 0;
  res_id = atoms[0].seq;
  //fprintf (stderr, "---- atoms ----- \n");

  if (number_of_atoms == 1) {
    starti = 2;
    residues[0].id = res_id;
    residues[0].ctype = atoms[0].ctype;
    residues[0].mtype = atoms[0].mtype;
    residues[0].num_atoms = 1;
    residues[0].atoms = new PmAtom*[1]; 
    residues[0].atoms[0] = &atoms[0]; 
    num_residues = 1;
    }
  else {
    starti = 1;
    }

  for (int i = starti; i < number_of_atoms; i++) {
    //fprintf (stderr, "%d  seq[%d] \n", i, atoms[i].seq);

    if ((atoms[i].seq != res_id) || (i == number_of_atoms-1)) {
      residues[num_residues].id = res_id;
      residues[num_residues].ctype = atoms[i-1].ctype;
      residues[num_residues].mtype = atoms[i-1].mtype;

      if ((atoms[i].seq == res_id) && (i == number_of_atoms-1)) {
        alist[n++] = i;
        }

      residues[num_residues].num_atoms = n;
      residues[num_residues].atoms = new PmAtom*[n]; 

      for (int j = 0; j < n; j++) {
        residues[num_residues].atoms[j] = &atoms[alist[j]];
        }

      num_residues += 1;
      alist[0] = i;
      n = 1;

      if ((atoms[i].seq != res_id) && (i == number_of_atoms-1)) {
        residues[num_residues].id = atoms[i].seq;
        residues[num_residues].ctype = atoms[i].ctype;
        residues[num_residues].mtype = atoms[i].mtype;
        residues[num_residues].num_atoms = n;
        residues[num_residues].atoms = new PmAtom*[n];

        for (int j = 0; j < n; j++) {
          residues[num_residues].atoms[j] = &atoms[alist[j]];
          }

        num_residues += 1;
        }

      res_id = atoms[i].seq;
      }
    else {
      alist[n++] = i;
      }
    }


  // create chain data //

  #ifdef dbg_buildChains
  fprintf (stderr, "  >>> num residues[%d] \n", num_residues);
  #endif

  if (num_residues == 0) {
    pm_ErrorReport (PM, "domain \"%s\" has no residues", "*", this->name.c_str());
    for (int i = 0; i < number_of_atoms; i++) {
      fprintf (stderr, ">>> atom  chain=%c seq=%d name=%s\n", atoms[i].chain, 
               atoms[i].seq, atoms[i].name);
      }
    return;
    }

  //===== create residue data =====//

  chains.resize (num_chains);
  num_chains = 0;
  chain_id = residues[0].atoms[0]->chain;
  n = 1;
  chains[0].res_index = 0;
  int first_res = 0;

  for (int i = 1; i < num_residues; i++) {
    if ((residues[i].atoms[0]->chain != chain_id) || (i == num_residues-1)) { 
      if (i == num_residues-1) {
        n += 1;
        }

      chains[num_chains].id = chain_id;
      chains[num_chains].num_residues = n;
      chains[num_chains].res_index = first_res;
      first_res = i;

      #ifdef dbg_buildChains
      fprintf (stderr, "  chain [%c]  num residues [%d] \n", chain_id, n);
      #endif
      num_chains += 1;
      chain_id = residues[i].atoms[0]->chain;
      n = 1;
      }
    else {
      n += 1;
      }
    }

  if (num_residues == 1) {
    chains[num_chains].id = chain_id;
    chains[num_chains].num_residues = 1;
    chains[num_chains].res_index = 0;
    num_chains += 1;
    }

  #ifdef dbg_buildChains
  fprintf (stderr, "----- chain info ----- \n");
  for (int i = 0; i < num_chains; i++) {
    PmChain *chain = &chains[i];
    fprintf (stderr, "  chain [%c]  num residues [%d]  res_index [%d] \n", 
             chain->id, chain->num_residues, chain->res_index);
    }
  #endif

  if (!parent) {
    return;
    }

  // set atom index. its id will not change. // 
  
  for (int i = 0; i < number_of_atoms; i++) {
    atoms[i].index = i;
    }
  }

}

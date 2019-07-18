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
//* atoms:                 a t o m s                           * 
//*============================================================*

#include "atoms.h"

namespace ProteinMechanica {


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////


//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmAtoms::PmAtoms() { 
  number_of_atoms = 0;
  num_chains = 0;
  atoms.resize (100);
  }

PmAtoms::PmAtoms(int n) { 
  number_of_atoms = 0;
  num_chains = 0;
  atoms.resize (n);
  }

PmAtoms::~PmAtoms() { }

//*============================================================*
//*==========              addAtom                   ==========*
//*============================================================*
// add an atom. 

void
PmAtoms::addAtom (const int id, PmAtomName name, PmAtomElementType type, PmVector3 pos,
                  char chain, PmCompoundType ctype, PmMoleculeType mtype, int seq) 
  { 
  //fprintf (stderr, ">>>>>> PmAtoms::add:  id [%d] name [%s] \n", id, name);
  //fprintf (stderr, "   >>> seq[%d] \n", seq);

  resize ();

  atoms[number_of_atoms].setData(id, name, type, chain, seq, ctype, mtype, pos);

  if (number_of_atoms == 0) {
    extent.min = pos;
    extent.max = pos;
    }
  else {
    extent.update (pos);
    }

  atoms[number_of_atoms].index = number_of_atoms; 
  number_of_atoms += 1;
  }

//*============================================================*
//*==========              addAtom                   ==========*
//*============================================================*
// add an atom using a ptr to an atom.

void 
PmAtoms::addAtom (PmAtom *atom)
  {
  resize ();
  atoms[number_of_atoms] = *atom;

  if (number_of_atoms == 0) {
    extent.min = atom->pos;
    extent.max = atom->pos;
    }
  else {
    extent.update (atom->pos);
    }

  /*
  fprintf (stderr, ">>>>>> PmAtoms::add:  id [%d] name [%s]  pos (%g %g %g) \n", 
           atom->id, atom->name, atom->pos[0], atom->pos[1], atom->pos[2]);
  */
  number_of_atoms += 1;
  }

//*============================================================*
//*==========              getAtom                   ==========*
//*============================================================*
// get an atom using an integer id.

void 
PmAtoms::getAtom (const int id, PmAtom **atom)
  {
  *atom = NULL;

  for (int j = 0; j < number_of_atoms; j++) {
    if (id == atoms[j].id) {
      *atom = &atoms[j];
      }
    }
  }

//*============================================================*
//*==========              getNumAtoms               ==========*
//*============================================================*
// get the number of atoms. 

int 
PmAtoms::getNumAtoms () {
  return (number_of_atoms);
  }


// get the number of atoms with the given <name>. 

int 
PmAtoms::getNumAtoms (PmAtomName name) 
  {

  int num = 0;

  //fprintf (stderr, ">>>>>> PmAtoms::getNumAtoms:  name [%s] \n", name);

  for (int i = 0; i < number_of_atoms; i++) {
    if (!atoms[i].compare(name)) {
      num += 1;
      }
    }

  return (num);
  }

//*============================================================*
//*==========             getAtomCoords              ==========*
//*============================================================*
// get the coordinates for all atoms.

void
PmAtoms::getAtomCoords (int *num, PmVector3 **coords) 
  {

  *num = 0;
  *coords = NULL;

  if (number_of_atoms == 0) {
    return;
    }

  PmVector3 *atom_coords = new PmVector3[number_of_atoms];

  for (int i = 0; i < number_of_atoms; i++) {
    atom_coords[i] = atoms[i].pos; 
    }

  *num = number_of_atoms;
  *coords = atom_coords;
  }


// get the coordinates for atoms by name.

void
PmAtoms::getAtomCoords (PmAtomName name, int *num_coords, PmVector3 **coords) 
  {

  *num_coords = 0;
  *coords = NULL;

  if (number_of_atoms == 0) {
    return;
    }

  int num = getNumAtoms (name);

  if (num == 0) {
    return;
    }

  PmVector3 *vec = new PmVector3 [num];
  num = 0;

  for (int i = 0; i < number_of_atoms; i++) {
    if (!atoms[i].compare(name)) {
      vec[num] = atoms[i].pos;
      num += 1;
      }
    }

  *num_coords = num;
  *coords = vec;
  }


// get the coordinates for all atoms returned in a vector.

void 
PmAtoms::getAtomCoords (vector<PmVector3> &coords) 
  {

  if (number_of_atoms == 0) {
    coords.clear();
    return;
    }

  coords.resize(number_of_atoms);

  for (int i = 0; i < number_of_atoms; i++) {
    coords[i] = atoms[i].pos;
    }
  }


// get the coordinates for all atoms returned in a vector //

void
PmAtoms::getAtomCoords (PmAtomName name, vector<PmVector3> &coords)
  {

  int num = getNumAtoms (name);

  if (num == 0) {
    coords.clear();
    return;
    }

  coords.resize(num);
  num = 0;

  for (int i = 0; i < number_of_atoms; i++) {
    if (!atoms[i].compare(name)) {
      coords[num++] = atoms[i].pos;
      }
    }
  }

//*============================================================*
//*==========              getNumChains              ==========*
//*============================================================*
// get the number of chains.

int
PmAtoms::getNumChains () 
  {

  if (num_chains) {
    return (num_chains);
    }

  char chain_id = atoms[0].chain;
  num_chains = 1;

  for (int i = 1; i < number_of_atoms; i++) {
    if (chain_id != atoms[i].chain) {
      num_chains += 1;
      chain_id = atoms[i].chain;
      }
    }

  return (num_chains);
  }

//*============================================================*
//*==========              getAtoms                  ==========*
//*============================================================*
// get a set of atoms by atom name.

void
PmAtoms::getAtoms (PmAtomName name, vector<PmAtom*>& sel_atoms)
  {

  int n = getNumAtoms(name);

  if (n == 0) {
    return;
    }

  sel_atoms.resize(n);
  n = 0;

  for (int i = 0; i < number_of_atoms; i++) {
    if (!atoms[i].compare(name)) {
      sel_atoms[n] = &atoms[i]; 
      n += 1;
      }
    }

  }

//*============================================================*
//*==========              getAtoms                  ==========*
//*============================================================*
// get all atoms.

void
PmAtoms::getAtoms (vector<PmAtom*>& sel_atoms)
  {

  sel_atoms.resize(number_of_atoms);

  for (int i = 0; i < number_of_atoms; i++) {
    sel_atoms[i] = &atoms[i];
    }
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get the mass properties of the atoms.

void 
PmAtoms::getMassProps (PmMassProperties& props) 
  {

  float x, y, z, cx, cy, cz, cmx, cmy, cmz;
  PmAtomElementType type;
  float total_mass, mass;
  PmMatrix3x3 inertia;
  float xx, xy, xz, yy, yz, zz;
  float scale; 

  fprintf (stderr, ">>>>>> PmAtoms::getMassProps \n");
  scale = 1.0 / pmSystem.getUnitsMassScale();
  cmx = cmy = cmz = 0.0;
  cx = cy = cz = 0.0;
  total_mass = 0.0;

  //  first compute center of mass // 

  for (int i = 0; i < number_of_atoms; i++) {
    x = atoms[i].pos[0];
    y = atoms[i].pos[1];
    z = atoms[i].pos[2];
    type = atoms[i].type;
    mass = scale * atoms[i].getMass();
    total_mass += mass;
    cmx += x*mass;
    cmy += y*mass;
    cmz += z*mass;
    cx += x;
    cy += y;
    cz += z;
    }

  cmx /= total_mass; cmy /= total_mass; cmz /= total_mass;
  cx  /= number_of_atoms;  cy  /= number_of_atoms;  cz  /= number_of_atoms;

  //  compute inertia tensor // 

  xx = xy = xz = yy = yz = zz = 0.0;

  for (int i = 0; i < number_of_atoms; i++) {
    x = atoms[i].pos[0] - cmx;
    y = atoms[i].pos[1] - cmy;
    z = atoms[i].pos[2] - cmz;
    type = atoms[i].type;
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

  props.mass = total_mass;
  props.com.set(cmx, cmy, cmz);
  props.inertia = inertia;
  }

//*============================================================*
//*==========        compPrincipalComponents         ==========*
//*============================================================*
// compute principal components.

void
PmAtoms::compPrincipalComponents(PmAtomFilter& filter, PmPcaResults& pca)
  {
  //fprintf (stderr, ">>> PmAtoms::compPrincipalComponents:\n");
  int n = 0;

  for (int i = 0; i < number_of_atoms; i++) {
    if (filter.hasName(atoms[i].name)) {
      n += 1;
      }
    }

  if (!n) {
    return;
    }

  PmVector3 *coords = new PmVector3[n];
  n = 0;

  for (int i = 0; i < number_of_atoms; i++) {
    if (filter.hasName(atoms[i].name)) {
      coords[n] = atoms[i].pos;
      n += 1;
      }
    }

  //fprintf (stderr, ">>> num coords = %d \n", n);
  pm_MathPrincipalComponents(n, coords, pca);
  delete[] coords;
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              resize                    ==========*
//*============================================================*
// resize atoms if needed.

void
PmAtoms::resize ()
  {

  int size = atoms.size();

  if (number_of_atoms == size) {
    size += 100;
    atoms.resize (size);
    }
  }


}

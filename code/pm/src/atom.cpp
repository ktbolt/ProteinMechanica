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
//* atom:                 a t o m                              *
//*============================================================*

#include "pm/pm.h" 
#include "atom_prv.h" 

namespace ProteinMechanica {

////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////


//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*


//*============================================================*
//*==========              setData                   ==========*
//*============================================================*
// set the data for an atom. 

void 
PmAtom::setData (int id, PmAtomName name, PmAtomElementType type, char chain, int seq,
                 PmCompoundType ctype, PmMoleculeType mtype, PmVector3 pos) {
  this->id = id;
  strcpy (this->name, name);
  this->type = type;
  this->chain = chain;
  this->seq = seq;
  this->ctype = ctype;
  this->mtype = mtype;
  this->pos = pos;
  }

//*============================================================*
//*==========              getColor                  ==========*
//*============================================================*
// get an atom's color by element type.

void
PmAtom::getColor (PmVector3& color) {
  color.set(atom_colors[type]);
  }

//*============================================================*
//*==========              getMass                   ==========*
//*============================================================*
// get an atom's mass by element type.

float 
PmAtom::getMass () {
  return atomic_props[type].weight;
  }

//*============================================================*
//*==========              getRadius                 ==========*
//*============================================================*
// get an atom's readius by element type.

float 
PmAtom::getRadius () {
  return atomic_props[type].radius;
  }

//*============================================================*
//*==========              convType                  ==========*
//*============================================================*
// convert an atom element type into a string. 

void
PmAtom::convType(PmAtomElementType type, char str[3])
  {
  str[0] = ' ';
  str[1] = '\0';
  str[2] = '\0';

  switch (type) {
    case PM_ATOM_ELEMENT_UNKNOWN:
    break;

    case PM_ATOM_ELEMENT_CARBON:
      str[0] = 'C';
    break;

    case PM_ATOM_ELEMENT_HYDROGEN:
      str[0] = 'H';
    break;

    case PM_ATOM_ELEMENT_NITROGEN:
      str[0] = 'N';
    break;

    case PM_ATOM_ELEMENT_OXYGEN:
      str[0] = 'O';
    break;

    case PM_ATOM_ELEMENT_SULFUR:
      str[0] = 'S';
    break;

    case PM_ATOM_ELEMENT_PHOSPHORUS:
      str[0] = 'P';
    break;

    case PM_ATOM_ELEMENT_IRON:
      str[0] = 'F';
      str[1] = 'E';
    break;

    case PM_ATOM_ELEMENT_SELENIUM:
      str[0] = 'S';
      str[1] = 'E';
    break;

    default:
      str[0] = '\0';
    }
  }

//*============================================================*
//*==========              convType                  ==========*
//*============================================================*
// ---------------
// convert a str atom name into an enum.  

PmAtomElementType
PmAtom::convType(char *s) 
  {
  PmAtomElementType type;

  while (*s == ' ') s++;

  if ((*s == 'C') || (*s == 'c')) {
    type = PM_ATOM_ELEMENT_CARBON;
    }
  else if ((*s == 'N') || (*s == 'n')) {
    type = PM_ATOM_ELEMENT_NITROGEN;
    }
  else if (*s == 'O') {
    type = PM_ATOM_ELEMENT_OXYGEN;
    }
  else if (*s == 'H') {
    type = PM_ATOM_ELEMENT_HYDROGEN;
    }
  else if (*s == 'S') {
    if (*(s+1) == 'E') {
      type = PM_ATOM_ELEMENT_SELENIUM;
      }
    else {
      type = PM_ATOM_ELEMENT_SULFUR;
      }
    }
  else if (*s == 'P') {
    type = PM_ATOM_ELEMENT_PHOSPHORUS;
    }
  else {
    type = PM_ATOM_ELEMENT_UNKNOWN;
    }

  return (type);
  }

//*============================================================*
//*==========              hasName                   ==========*
//*============================================================*
// search for an atom name <str> defined in an atom filter.

bool
PmAtomFilter::hasName(const string str)
  {
  if (!names.size()) {
    return true;
    }

  bool has = false;

  for (int i = 0; i < (int)names.size(); i++) {
    if (str == names[i]) { 
      has = true;
      break;
      }
    }

  if (exclude) {
    return !has;
    }
  else {
    return has;
    }
  }

}

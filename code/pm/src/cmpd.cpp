
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
 * cmpd:             c o m p o u n d s                        *
 *============================================================*/

#include "cmpd.h" 

namespace ProteinMechanica {

bool pm_CheckAtomName(char *name) { 
  bool found = false;

  for (int i = 0; nuc_names[i] != 0; i++) {
      if (!strcmp(name,nuc_names[i])) {
          fprintf (stderr, " **** nuc atom name [%s] \n", name);
          return true;
      }
  }

  return false;
  }

void
pm_CmpdConvNucleotideType (char *s, PmCompoundType& ctype, PmMoleculeType& mtype) {

    if (*s == 'A') {
      ctype = PM_COMPOUND_ADENOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (*s == 'C') {
      ctype = PM_COMPOUND_CYTOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (*s == 'G') {
      ctype = PM_COMPOUND_GUANOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (*s == 'U') {
      ctype = PM_COMPOUND_URIDINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (*s == 'T') {
      ctype = PM_COMPOUND_THYMINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

}

//*============================================================*
//*==========              pm_CmpdConvType           ==========*
//*============================================================*
// convert a compound string name into a enum.

void
pm_CmpdConvType (char *s, PmCompoundType& ctype, PmMoleculeType& mtype)
  { 
  //fprintf (stderr, ">>>>>> pm_CmpdConvType: %s \n", s);
  int i;
  char tstr[4];
  bool found=false;

  ctype = PM_COMPOUND_UNKNOWN;
  mtype = PM_MOLECULE_UNKNOWN;

  while (*s == ' ') {
    s++;
    }

  for (unsigned i = 0; i < strlen(s); i++) {
    if (s[i] == ' ') {
      s[i] = '\0';
      break;
      }
    }

  //===== check for single character symbols =====//

  /* NOTE[Feb20,2016] comment out single char for proteins so we can process
                      single chars for rna. 
  */

  if (strlen(s) == 1) {
     pm_CmpdConvNucleotideType(s, ctype, mtype);
     return;
  }

  /*
  if (strlen(s) == 1) {
    for (i = 1; cmpd_names[i].name; i++) {
      if (*s == *cmpd_names[i].abrv1) {
        ctype = cmpd_names[i].ctype; 
        mtype = cmpd_names[i].mtype; 
        found = true;
        break;
        }

      if (tolower(*s) == *cmpd_names[i].abrv1) {
        ctype = cmpd_names[i].ctype; 
        mtype = cmpd_names[i].mtype; 
        found = true;
        break;
        }
      }
    }

  if (found) {
    return;
    }
*/

  for (i = 1; cmpd_names[i].name; i++) {
    if (!strcasecmp(s, cmpd_names[i].name) ||
        !strcasecmp(s, cmpd_names[i].abrv3) || 
        !strcasecmp(s, cmpd_names[i].abrv1)) { 
      ctype = cmpd_names[i].ctype; 
      mtype = cmpd_names[i].mtype; 
      }
    }


  // nucleic acid? //

  if (strlen(s) == 3) {
    if (s[1] == 'A') {
      ctype = PM_COMPOUND_ADENOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'C') {
      ctype = PM_COMPOUND_CYTOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'G') {
      ctype = PM_COMPOUND_GUANOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'U') {
      ctype = PM_COMPOUND_URIDINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'T') {
      ctype = PM_COMPOUND_THYMINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }
    }

  if (strlen(s) == 2) {
    if (s[1] == 'A') {
      ctype = PM_COMPOUND_ADENOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'C') {
      ctype = PM_COMPOUND_CYTOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'G') {
      ctype = PM_COMPOUND_GUANOSINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'U') {
      ctype = PM_COMPOUND_URIDINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }

    else if (s[1] == 'T') {
      ctype = PM_COMPOUND_THYMINE;
      mtype = PM_MOLECULE_NUCLEOTIDE;
      }
    }

  // check for gromacs modified terminal residue names //

  if (strlen(s) == 3) {
    tstr[0] = s[0];
    tstr[1] = s[1];
    tstr[2] = s[2];
    tstr[3] = s[3];
    
    if (!strcasecmp("lyp", s)) {
      tstr[2] = 's';
      }

    else if (!strcasecmp("hip", s)) {
      tstr[2] = 's';
      }

    else if (!strcasecmp("hid", s)) {
      tstr[2] = 's';
      }

    else if (!strcasecmp("hie", s)) {
      tstr[2] = 's';
      }

    else if (!strcasecmp("cyn", s)) {
      tstr[2] = 's';
      }

    for (i = 1; cmpd_names[i].name; i++) {
      if (!strcasecmp(tstr, cmpd_names[i].abrv3)) { 
        ctype = cmpd_names[i].ctype; 
        mtype = cmpd_names[i].mtype; 
        }
      }
    }

  if (strlen(s) == 4) {
    tstr[0] = s[1];
    tstr[1] = s[2];
    tstr[2] = s[3];
    tstr[3] = '\0';

    if (!strcasecmp("lyp", tstr)) { 
      tstr[2] = 's';
      }

    else if (!strcasecmp("hip", tstr)) { 
      tstr[2] = 's';
      }

    else if (!strcasecmp("hid", tstr)) { 
      tstr[2] = 's';
      }

    else if (!strcasecmp("hie", tstr)) { 
      tstr[2] = 's';
      }

    else if (!strcasecmp("cyn", tstr)) { 
      tstr[2] = 's';
      }

    for (i = 1; cmpd_names[i].name; i++) {
      if (!strcasecmp(tstr, cmpd_names[i].abrv3)) { 
        ctype = cmpd_names[i].ctype; 
        mtype = cmpd_names[i].mtype; 
        }
      }
    }

  if (strlen(s) == 4) {
    tstr[0] = s[0];
    tstr[1] = s[1];
    tstr[2] = s[2];
    tstr[3] = '\0';

    for (i = 1; cmpd_names[i].name; i++) {
      if (!strcasecmp(tstr, cmpd_names[i].abrv3)) { 
        ctype = cmpd_names[i].ctype; 
        mtype = cmpd_names[i].mtype; 
        }
      }
    }

  //fprintf (stderr, "   >>> CmpdType: %d \n", ctype);
  }

//*============================================================*
//*==========              pm_CmpdGetName            ==========*
//*============================================================*
// get a compound name from an enum.

void
pm_CmpdGetName (PmCompoundType type, char name[3], int n) 
  {
  for (int i = 1; cmpd_names[i].name; i++) {
    if (type == cmpd_names[i].ctype) {
      if (n == 1) {
        strcpy (name, cmpd_names[i].abrv1);
        }
      else {
        strcpy (name, cmpd_names[i].abrv3);
        }
      return;
      }
    }

  strcpy (name, cmpd_names[0].abrv3);
  }

//*============================================================*
//*==========              pm_AminoGetProp           ==========*
//*============================================================*
// get an amino acid property. 

void
pm_AminoGetProp (PmCompoundType type, PmAminoAcidPropType& prop) {
  prop = amino_props[type];
  }

//*============================================================*
//*==========              pm_CmpdConvType           ==========*
//*============================================================*
// convert a compound string name into a enum.

void
pm_CmpdMolTypeToString (PmMoleculeType mtype, string& str)
  {
  switch (mtype) {
    case PM_MOLECULE_ANY:
      str = "any";
    break;

    case PM_MOLECULE_AMINO_ACID:
      str = "protein";
    break;

    case PM_MOLECULE_ION:
      str = "ion";
    break;

    case PM_MOLECULE_NUCLEOTIDE:
      str = "nucleotide";
    break;

    case PM_MOLECULE_PROTEIN:
      str = "protein";
    break;

    case PM_MOLECULE_WATER:
      str = "water";
    break;

    default:
      str = "unknown";
    }
  }
}



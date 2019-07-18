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
 * amino:    a m i n o  a c i d s  a n d  b a s e s           *
 *============================================================*/

#include "amino.h" 

namespace ProteinMechanica {

// pm_AminoTypeConv 
// ----------------

PmAminoAcidType
pm_AminoTypeConv (char *s) 
  { 
  PmAminoAcidType type;
  int i;
  char tstr[4];

  type = PM_AMINO_ACID_UNKNOWN;

  while (*s == ' ') {
    s++;
    }

  if (strlen(s) == 1) {
    for (i = 1; amino_names[i][0]; i++) {

      if (*s == *amino_names[i][2]) {
        type = (PmAminoAcidType)i;
        return (type);
        }

      if (tolower(*s) == *amino_names[i][2]) {
        type = (PmAminoAcidType)i;
        return (type);
        }
      }
    }

  for (i = 1; amino_names[i][0]; i++) {
    if (!strcasecmp(s, amino_names[i][0]) ||
        !strcasecmp(s, amino_names[i][1]) || 
        !strcasecmp(s, amino_names[i][2])) { 
      type = (PmAminoAcidType)i;
      return (type);
      }
    }


  /*  nucleic acid?  */

  if (strlen(s) == 3) {
    if (s[1] == 'A') {
      type = PM_AMINO_ACID_ADENOSINE;
      return (type);
      }

    else if (s[1] == 'C') {
      type = PM_AMINO_ACID_CYTOSINE;
      return (type);
      }

    else if (s[1] == 'G') {
      type = PM_AMINO_ACID_GUANOSINE;
      return (type);
      }

    else if (s[1] == 'U') {
      type = PM_AMINO_ACID_URIDINE;
      return (type);
      }

    else if (s[1] == 'T') {
      type = PM_AMINO_ACID_THYMINE;
      return (type);
      }
    }


  /*  check for gromacs modified terminal residue names.  */

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

    for (i = 1; amino_names[i][0]; i++) {
      if (!strcasecmp(tstr, amino_names[i][1])) { 
        type = (PmAminoAcidType)i;
        return (type);
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

    for (i = 1; amino_names[i][0]; i++) {
      if (!strcasecmp(tstr, amino_names[i][1])) { 
        type = (PmAminoAcidType)i;
        return (type);
        }
      }
    }

  if (strlen(s) == 4) {
    tstr[0] = s[0];
    tstr[1] = s[1];
    tstr[2] = s[2];
    tstr[3] = '\0';

    for (i = 1; amino_names[i][0]; i++) {
      if (!strcasecmp(tstr, amino_names[i][1])) { 
        type = (PmAminoAcidType)i;
        return (type);
        }
      }
    }
  }


// pm_AminoNameGet  
// ---------------

void
pm_AminoNameGet (PmAminoAcidType type, int i, char *name) {
  strcpy (name, amino_names[type][i]);
  }


// pm_AminoPropGet  
// ---------------

void
pm_AminoPropGet (PmAminoAcidType type, PmAminoAcidPropType& prop) {
  prop = amino_props[type];
  }


#ifdef use

/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_AminoCharged  ****               *
 *                                                            *
 *------------------------------------------------------------*/

int
pm_AminoCharged (PmAminoAcidType type) {
  return ( (amino_props[type] == PM_AMINO_ACID_PROP_ACIDIC) ||
           (amino_props[type] == PM_AMINO_ACID_PROP_BASIC)); 
  }

int
pm_AminoBasic (PmAminoAcidType type) {
  return ( amino_props[type] == PM_AMINO_ACID_PROP_BASIC);
  }


int
pm_AminoAcidic(PmAminoAcidType type) {
  return ( amino_props[type] == PM_AMINO_ACID_PROP_ACIDIC ); 
  }

int
pm_AminoPolar(PmAminoAcidType type) {
  return ( amino_props[type] == PM_AMINO_ACID_PROP_POLAR); 
  }

int
pm_AminoHydrop(PmAminoAcidType type) {
  return ( amino_props[type] == PM_AMINO_ACID_PROP_HYDROPHOBIC); 
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_AminoPropCheck  ****             *
 *                                                            *
 *------------------------------------------------------------*/

int
pm_AminoPropCheck (PmAminoAcidType type, PmAminoAcidPropType prop)
  {

 /**************
  ***  body  ***
  **************/

  return ( (amino_props[type] == prop) );
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_AminoName2Amber  ****            *
 *                                                            *
 * convert a residue name into gromacs amber name.            *
 *------------------------------------------------------------*/

void
pm_AminoName2Amber (PmAminoAcidType type, char *name, int term, char *aname)
  {

  char last;

  int i;

  char p;

 /**************
  ***  body  ***
  **************/

  p ='\0'; 

  if (type == PM_AMINO_ACID_LYSINE) {
    strcpy (aname, "LYP");
    }
  else if (type == PM_AMINO_ACID_CYSTEINE) {
    strcpy (aname, "CYN");
    }
  else if (type == PM_AMINO_ACID_HISTIDINE) {
    strcpy (aname, "HIP");
    }
  else {
    strcpy (aname, name);
    }

  if (term == 1) {
    p = 'N'; 
    }

  else if (term == 2) {
    p = 'C'; 
    }

  if (p != '\0') {
    for (i = 4; i > 0 ; i--) {
      aname[i] = aname[i-1];
      }

    aname[0] = p;
    }
  }

#endif

};




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

#include "pm/pm.h"

namespace ProteinMechanica {

typedef struct CmpdName {
  char *name, *abrv3, *abrv1;
  PmCompoundType ctype;
  PmMoleculeType mtype;
  } CmpdName;

static char *nuc_names[] = { "OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", 
                             "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", 
                             "N2", "N3", "C4", 0 };

static CmpdName cmpd_names[] = {
  {"unknown",       "ukn",     "?", PM_COMPOUND_UNKNOWN, PM_MOLECULE_UNKNOWN},
  {"alanine",       "ala",     "a", PM_COMPOUND_ALANINE, PM_MOLECULE_AMINO_ACID},   
  {"arginine",      "arg",     "r", PM_COMPOUND_ARGININE, PM_MOLECULE_AMINO_ACID},
  {"asparagine",    "asn",     "n", PM_COMPOUND_ASPARAGINE, PM_MOLECULE_AMINO_ACID},
  {"aspartic acid", "asp",     "d", PM_COMPOUND_ASPARTIC_ACID, PM_MOLECULE_AMINO_ACID},
  {"cysteine",      "cys",     "c", PM_COMPOUND_CYSTEINE, PM_MOLECULE_AMINO_ACID},
  {"glutamic acid", "glu",     "e", PM_COMPOUND_GLUTAMIC_ACID, PM_MOLECULE_AMINO_ACID},
  {"glutamine",     "gln",     "q", PM_COMPOUND_GLUTAMINE, PM_MOLECULE_AMINO_ACID},
  {"glycine",       "gly",     "g", PM_COMPOUND_GLYCINE, PM_MOLECULE_AMINO_ACID},
  {"histidine",     "his",     "h", PM_COMPOUND_HISTIDINE, PM_MOLECULE_AMINO_ACID},
  {"isoleucine",    "ile",     "i", PM_COMPOUND_ISOLEUCINE, PM_MOLECULE_AMINO_ACID},
  {"leucine",       "leu",     "l", PM_COMPOUND_LEUCINE, PM_MOLECULE_AMINO_ACID},
  {"lysine",        "lys",     "k", PM_COMPOUND_LYSINE, PM_MOLECULE_AMINO_ACID},
  {"methionine",    "met",     "m", PM_COMPOUND_METHIONINE, PM_MOLECULE_AMINO_ACID},
  {"selenomethionine",  "mse", "m", PM_COMPOUND_SELMETHIONINE, PM_MOLECULE_AMINO_ACID},
  {"phenylalanine", "phe",     "f", PM_COMPOUND_PHENYLALANINE, PM_MOLECULE_AMINO_ACID},
  {"proline",       "pro",     "p", PM_COMPOUND_PROLINE, PM_MOLECULE_AMINO_ACID},
  {"serine",        "ser",     "s", PM_COMPOUND_SERINE, PM_MOLECULE_AMINO_ACID},
  {"threonine",     "thr",     "t", PM_COMPOUND_THREONINE, PM_MOLECULE_AMINO_ACID},
  {"tryptophan",    "trp",     "w", PM_COMPOUND_TRYPTOPHAN, PM_MOLECULE_AMINO_ACID},
  {"tyrosine",      "tyr",     "y", PM_COMPOUND_TYROSINE, PM_MOLECULE_AMINO_ACID},
  {"valine",        "val",     "v", PM_COMPOUND_VALINE, PM_MOLECULE_AMINO_ACID}, 

  {"adenosine",     "adn",     "A", PM_COMPOUND_ADENOSINE, PM_MOLECULE_NUCLEOTIDE},
  {"cytosine",      "cyt",     "C", PM_COMPOUND_CYTOSINE, PM_MOLECULE_NUCLEOTIDE},
  {"guanosine",     "gua",     "G", PM_COMPOUND_GUANOSINE, PM_MOLECULE_NUCLEOTIDE},
  {"thymine",       "thy",     "T", PM_COMPOUND_THYMINE, PM_MOLECULE_NUCLEOTIDE},
  {"uridine",       "ura",     "U", PM_COMPOUND_URIDINE, PM_MOLECULE_NUCLEOTIDE},
  {"adp",           "adp",     "h", PM_COMPOUND_ADP, PM_MOLECULE_NUCLEOTIDE},
  {"atp",           "atp",     "h", PM_COMPOUND_ATP, PM_MOLECULE_NUCLEOTIDE},

  {"adenosine",     "adn",     "DA", PM_COMPOUND_ADENOSINE, PM_MOLECULE_NUCLEOTIDE},
  {"cytosine",      "cyt",     "DC", PM_COMPOUND_CYTOSINE, PM_MOLECULE_NUCLEOTIDE},
  {"guanosine",     "gua",     "DG", PM_COMPOUND_GUANOSINE, PM_MOLECULE_NUCLEOTIDE},
  {"thymine",       "thy",     "DT", PM_COMPOUND_THYMINE, PM_MOLECULE_NUCLEOTIDE},

  {"cl",            "cli",     "h", PM_COMPOUND_CHLORINE, PM_MOLECULE_ION},
  {"na",            "cli",     "h", PM_COMPOUND_SODIUM, PM_MOLECULE_ION},

  {"hoh",           "hoh",     "h", PM_COMPOUND_WATER, PM_MOLECULE_WATER},

  {0,              "",        "",  PM_COMPOUND_UNKNOWN, PM_MOLECULE_UNKNOWN}};

 static 
 PmAminoAcidPropType amino_props[22] = {
    PM_AMINO_ACID_PROP_UNKNOWN,
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* ala  */ 
    PM_AMINO_ACID_PROP_BASIC,          /* arg  */
    PM_AMINO_ACID_PROP_POLAR,          /* asn  */
    PM_AMINO_ACID_PROP_ACIDIC,         /* asp  */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* cys  */
    PM_AMINO_ACID_PROP_ACIDIC,         /* glu  */
    PM_AMINO_ACID_PROP_POLAR,          /* gln  */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* gly  */
    PM_AMINO_ACID_PROP_BASIC,          /* his  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* ile  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* leu  */
    PM_AMINO_ACID_PROP_BASIC,          /* lys  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* met  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* mse  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* phe  */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* pro  */
    PM_AMINO_ACID_PROP_POLAR,          /* ser  */
    PM_AMINO_ACID_PROP_POLAR,          /* thr  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* trp  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* tyr  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* val  */
    };

}

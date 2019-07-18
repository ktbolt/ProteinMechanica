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

static char *amino_names[PM_AMINO_ACID_MAX_NUM+1][3] = {
  {"unknown",       "unknown", "unknown"},
  {"alanine",       "ala",     "a"},
  {"arginine",      "arg",     "r"},
  {"asparagine",    "asn",     "n"},
  {"aspartic acid", "asp",     "d"},
  {"cysteine",      "cys",     "c"},
  {"glutamic acid", "glu",     "e"},
  {"glutamine",     "gln",     "q"},
  {"glycine",       "gly",     "g"},
  {"histidine",     "his",     "h"},
  {"isoleucine",    "ile",     "i"},
  {"leucine",       "leu",     "l"},
  {"lysine",        "lys",     "k"},
  {"methionine",    "met",     "m"},
  {"phenylalanine", "phe",     "f"},
  {"proline",       "pro",     "p"},
  {"serine",        "ser",     "s"},
  {"threonine",     "thr",     "t"},
  {"tryptophan",    "trp",     "w"},
  {"tyrosine",      "tyr",     "y"},
  {"valine",        "val",     "v"},
  {"solvent",       "sol",     "h"},
  {"adp",           "adp",     "h"},
  {"cl",            "cli",     "h"},
  {"hoh",           "hoh",     "h"},
  {"adenosine",     "adn",     "A"},
  {"guanosine",     "gua",     "G"},
  {"cytosine",      "cyt",     "C"},
  {"uridine",       "ura",     "U"},
  {"thymine",       "thy",     "T"},
  {"uridine",       "h2u",     "U"},
  {"cytosine",      "omc",     "C"},
  {"guanosine",     "omg",     "G"},
  {"uridine",       "psu",     "U"},
  {"cytosine",      "5mc",     "C"},
  {"guanosine",     "7mg",     "G"},
  {"uridine",       "5mu",     "U"},
  {"adenosine",     "1ma",     "A"},
  {"guanosine",     "2mg",     "G"},
  {"guanosine",     "m2g",     "G"},
  {"",              "",        ""}};

 static PmAminoAcidPropType amino_props[PM_AMINO_ACID_MAX_NUM] = {
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
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* phe  */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* pro  */
    PM_AMINO_ACID_PROP_POLAR,          /* ser  */
    PM_AMINO_ACID_PROP_POLAR,          /* thr  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* trp  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* tyr  */
    PM_AMINO_ACID_PROP_HYDROPHOBIC,    /* val  */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* sol  */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* aden */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* guan */
    PM_AMINO_ACID_PROP_UNKNOWN,        /* cyto */
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN,        
    PM_AMINO_ACID_PROP_UNKNOWN};


}

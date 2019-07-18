
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
//* pm:             p r o t e i n   m o d e l e r              *
//*============================================================*

#ifndef _PM_H_
#define _PM_H_

//#ifndef  USE_GRAPHICS 
//#define  USE_GRAPHICS  1
//#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "pm/common.h"

using namespace std;

namespace ProteinMechanica {

#define  PM "PM"


//////////////////////////////////////////////////////////////////
//                   b a s i c    t y p e s                    //
////////////////////////////////////////////////////////////////

const int pmMaxAtomsPerResidue = 300;

typedef enum PmDbType {
  PM_DB_UNKNOWN,
  PM_DB_CERFACS,
  PM_DB_PDB,
  PM_DB_GROMACS,
  PM_DB_MRC,
  PM_DB_PM,
  PM_DB_TEXT,
  PM_DB_COORD
  } PmDbType;

typedef enum PmDbModeType {
  PM_DB_MODE_UNKNOWN,
  PM_DB_MODE_READ,
  PM_DB_MODE_WRITE,
  PM_DB_MODE_APPEND
  } PmDbModeType;

typedef enum PmMoleculeType {
  PM_MOLECULE_UNKNOWN,
  PM_MOLECULE_ANY,
  PM_MOLECULE_AMINO_ACID,
  PM_MOLECULE_ION,
  PM_MOLECULE_MISCELLANEOUS,
  PM_MOLECULE_NUCLEIC_ACID,
  PM_MOLECULE_NUCLEOTIDE,
  PM_MOLECULE_NUCLEOSIDE,
  PM_MOLECULE_PEPTIDOGLYCAN,
  PM_MOLECULE_PROTEIN,
  PM_MOLECULE_WATER
  } PmMoleculeType;

typedef enum PmObjectType {
  PM_OBJECT_UNKNOWN,
  PM_OBJECT_DOMAIN,
  PM_OBJECT_MODEL,
  PM_OBJECT_MOLECULE,
  PM_OBJECT_PARTICLE,
  PM_OBJECT_SOLID,
  PM_OBJECT_SURFACE
  } PmObjectType;

typedef enum PmAminoAcidPropType {
  PM_AMINO_ACID_PROP_UNKNOWN,
  PM_AMINO_ACID_PROP_ACIDIC,
  PM_AMINO_ACID_PROP_BASIC,
  PM_AMINO_ACID_PROP_HYDROPHOBIC,
  PM_AMINO_ACID_PROP_POLAR
  } PmAminoAcidPropType;

typedef enum PmCompoundType {
  PM_COMPOUND_UNKNOWN,
  PM_COMPOUND_ALANINE,
  PM_COMPOUND_ARGININE,
  PM_COMPOUND_ASPARAGINE,
  PM_COMPOUND_ASPARTIC_ACID,
  PM_COMPOUND_CYSTEINE,
  PM_COMPOUND_GLUTAMIC_ACID,
  PM_COMPOUND_GLUTAMINE,
  PM_COMPOUND_GLYCINE,
  PM_COMPOUND_HISTIDINE,
  PM_COMPOUND_ISOLEUCINE,
  PM_COMPOUND_LEUCINE,
  PM_COMPOUND_LYSINE,
  PM_COMPOUND_METHIONINE,
  PM_COMPOUND_SELMETHIONINE,
  PM_COMPOUND_PHENYLALANINE,
  PM_COMPOUND_PROLINE,
  PM_COMPOUND_SERINE,
  PM_COMPOUND_THREONINE,
  PM_COMPOUND_TRYPTOPHAN,
  PM_COMPOUND_TYROSINE,
  PM_COMPOUND_VALINE,
  // nucleotides
  PM_COMPOUND_ADENOSINE,
  PM_COMPOUND_CYTOSINE,
  PM_COMPOUND_GUANOSINE,
  PM_COMPOUND_THYMINE,
  PM_COMPOUND_URIDINE,
  PM_COMPOUND_ADP,
  PM_COMPOUND_ATP,
  // ions
  PM_COMPOUND_CHLORINE,
  PM_COMPOUND_SODIUM,
  // water
  PM_COMPOUND_WATER,
  PM_COMPOUND_MAX_NUM
  } PmCompoundType;


// atom definitions
#include "atom.h"

// residue definitions
#include "res.h"


// mass properties

class PM_EXPORT PmMassProperties {
   public:
     PmMassProperties() { set = false; mass = 0.0; com.set(0,0,0); };
     ~PmMassProperties() {};
     bool set;
     float mass;
     PmMatrix3x3 inertia;
     PmVector3 com;

     PmMassProperties &operator=(const PmMassProperties& rhs) {
        set = rhs.set;
        mass = rhs.mass; 
        inertia = rhs.inertia; 
        com = rhs.com; 
        return *this;
        }
   };

//////////////////////////////////////////////////////////////////
//                   u n i t s                                 //
////////////////////////////////////////////////////////////////

const double PM_UNITS_AMU  = 1.6605402e-27;    // kilograms, aka a Dalton 
const double PM_UNITS_NANO = 1.0e-09;
const double PM_UNITS_PICO = 1.0e-12;
const double PM_CONS_KB    = 1.381e-23;        // boltzmann's const (J/K)
const double PM_CONS_ETA   = 1.0e-3;      

typedef struct PmUnit {
  string name;
  double value;
  } PmUnit;

typedef enum PmUnitBaseType {
  PM_UNIT_BASE_UNKNOWN,
  PM_UNIT_BASE_LEM,                  /*  length, energy, mass  */
  PM_UNIT_BASE_LMT                   /*  length, mass, time  */
  } PmUnitBaseType;

typedef enum PmUnitType {
  PM_UNIT_UNKNOWN,
  PM_UNIT_ENERGY,
  PM_UNIT_LENGTH,
  PM_UNIT_MASS,
  PM_UNIT_TEMPERATURE,
  PM_UNIT_TIME,
  PM_UNIT_MAX_NUM
  } PmUnitType;

typedef struct PmUnitVar {
  PmUnitType type;
  string name;
  float scale;
  PmUnit unit;
  } PmUnitVar;

typedef struct PmUnitsReduced {
  float kb, kT, eta;
  } PmUnitsReduced;

typedef struct PmUnits {
  PmUnitBaseType base;
  PmUnitVar energy;
  PmUnitVar length;
  PmUnitVar mass;
  PmUnitVar temperature;
  PmUnitVar time;
  PmUnitsReduced reduced;
  } PmUnits;

//////////////////////////////////////////////////////////////////
//                   c h a i n                                 //
////////////////////////////////////////////////////////////////

class PmChain {
  public:
    char id;
    PmMoleculeType type;
    int num_residues;
    int res_index;
    int num_atoms;
    PmAtom **atoms;

    PmChain () {
      id = '\0';
      num_residues = 0;
      res_index = 0;
      num_atoms = 0;
      atoms = NULL;
      }
  };

//////////////////////////////////////////////////////////////////
//                   m o l e c u l e                           //
////////////////////////////////////////////////////////////////

typedef enum {
  PM_MOLECULE_DISPLAY_UNKNOWN,
  PM_MOLECULE_DISPLAY_POINT,
  PM_MOLECULE_DISPLAY_CROSS,
  PM_MOLECULE_DISPLAY_SPHERE,
  PM_MOLECULE_DISPLAY_SIZE
  } PmMoleculeDisplayType;

typedef enum {
  PM_MOLECULE_RENDER_UNKNOWN,
  PM_MOLECULE_RENDER_POINT,
  PM_MOLECULE_RENDER_LINE,
  PM_MOLECULE_RENDER_SOLID,
  PM_MOLECULE_RENDER_SIZE
  } PmMoleculeRenderType;

class PmMoleculeDisplayAtts {
  public:
    PmMoleculeDisplayAtts() { init(); }
    PmMoleculeDisplayType type;
    PmAtomColorType color_type, bond_color_type;
    PmMoleculeRenderType render, bond_atom_render;
    PmVector3 color, atom_color, bond_color, surf_color; 
    bool bond_atoms, base_planes, cbonly;
    float marker_size, line_width, bond_atom_radius;
    int time_step;
    float dim;
    bool tube;
    int tube_res;
    float tube_width;

    void init() {
      color_type = PM_ATOM_COLOR_RGB;
      color.set(1, 1, 1);
      atom_color.set(1, 1, 1);
      bond_color.set(1, 1, 1);
      surf_color.set(1, 1, 1);
      time_step = 0;
      type = PM_MOLECULE_DISPLAY_SPHERE;
      dim = 0.0;
      render = PM_MOLECULE_RENDER_SOLID;
      marker_size = 0.1;
      line_width = 1.0;
      bond_atoms = false;
      bond_atom_radius = 0.01;
      bond_atom_render = PM_MOLECULE_RENDER_SOLID;
      base_planes = false;
      cbonly = false;
      tube = false;
      tube_res = 4;
      tube_width = 0.01;
      }
  };

typedef struct PmBond {
  PmAtom *atom1, *atom2;
  int ai1, ai2;
  } PmBond;

typedef struct PmMolPlane {
  vector<PmAtom*> atoms;
  vector<PmVector3> vertices;
  } PmMolPlane;


// local helix properties
//-----------------------

class PmHelixProps {
  public:
     PmHelixProps() { num = 0; }
     int num;
     int start_res_num;
     int end_res_num;
     char chain_id;
     vector<int> res_id;
     vector<PmVector3> origin;
     vector<PmVector3> axis;
     vector<PmVector3> ca_pos;

     vector<float> twist_angle;
     vector<float> res_per_turn;
     vector<float> radius;
     vector<float> height;
  };

class PM_EXPORT PmQuery {
  public:
    PmQuery(){}
    string name;
    PmVector3 point;
    int entity;
  };

//////////////////////////////////////////////////////////////////
//                   s i m u l a t i o n                       //
////////////////////////////////////////////////////////////////

class PM_EXPORT PmTimeInterval {
   public:
      PmTimeInterval() { interval_set = false; }
      bool interval_set;
      vector<float>begin_intervals; 
      vector<float>end_intervals;
      bool hasTime(float t);
      void addInterval(float begin, float end);
      bool set() { return interval_set; }
   };

typedef enum {
  PM_JOINT_UNKNOWN,
  PM_JOINT_NONE,
  PM_JOINT_BALL,
  PM_JOINT_FREE,
  PM_JOINT_GIMBAL,
  PM_JOINT_HINGE,
  PM_JOINT_SLIDER,
  PM_JOINT_UNIVERSAL,
  PM_JOINT_WELD,
  PM_JOINTSIZE
  } PmJointType;

class PmMode {
  public:
    PmMode() { size = 0; } 
    int size;
    float eigen_value;
    PmVector3 *eigen_vectors;
  };


//////////////////////////////////////////////////////////////////
//                   p r o t o t y p e s                       //
////////////////////////////////////////////////////////////////

PmAminoAcidPropType
pm_AminoGetProp(PmCompoundType type);

void PM_EXPORT
pm_CmpdConvType(char *s, PmCompoundType& ctype, PmMoleculeType& mtype);

void PM_EXPORT
pm_CmpdConvNucleotideType(char *s, PmCompoundType& ctype, PmMoleculeType& mtype);

void PM_EXPORT
pm_CmpdGetName(PmCompoundType type, char name[3], int n=3);

void PM_EXPORT
pm_CmpdMolTypeToString (PmMoleculeType mtype, string& str);

PmDbType PM_EXPORT
pm_DbFormatConv (const string& format);

void PM_EXPORT
pm_ErrorReport (char *mod, char *format, ...);

void PM_EXPORT
pm_ErrorWarnReport (char *mod, char *format, ...);

int PM_EXPORT
pm_FileExists (char *name);

void PM_EXPORT
pm_PrintMsg (const char *prompt, char *format, ...);

}


// pm system

#include "pmsys.h"

#endif



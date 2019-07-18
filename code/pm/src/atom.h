
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

#ifndef _ATOM_PM_H_
#define _ATOM_PM_H_

typedef enum {
  PM_ATOM_COLOR_UNKNOWN,
  PM_ATOM_COLOR_RGB,
  PM_ATOM_COLOR_ELEMENT,
  PM_ATOM_COLOR_CHAIN,
  PM_ATOM_COLOR_FORCE,
  PM_ATOM_COLOR_GFORCE,
  PM_ATOM_COLOR_MODE_ENERGY,
  PM_ATOM_COLOR_CHARGE,
  PM_ATOM_COLOR_CHARGE_POLAR,
  PM_ATOM_COLOR_RMSD,
  PM_ATOM_COLOR_HYDROPHOBIC,
  PM_ATOM_COLOR_TEMP
  } PmAtomColorType;

typedef enum PmAtomElementType {
  PM_ATOM_ELEMENT_UNKNOWN,
  PM_ATOM_ELEMENT_CARBON,
  PM_ATOM_ELEMENT_HYDROGEN,
  PM_ATOM_ELEMENT_NITROGEN,
  PM_ATOM_ELEMENT_OXYGEN,
  PM_ATOM_ELEMENT_SULFUR,
  PM_ATOM_ELEMENT_PHOSPHORUS,
  PM_ATOM_ELEMENT_IRON,
  PM_ATOM_ELEMENT_SELENIUM,
  PM_ATOM_ELEMENT_ALL,
  PM_ATOM_ELEMENT_MAX
  } PmAtomElementType;

//===== atomic properties =====//

typedef struct PmAtomicProps {
  float radius;
  float mass;
  float weight;
  } PmAtomicProps;

//typedef std::string (PmAtomName);
typedef char (PmAtomName[5]);


// atom class //

class PM_EXPORT PmAtom {

  public:
    int id, index;
    PmAtomName name;
    int seq;
    char chain;
    PmVector3 pos;
    PmAtomElementType type;
    PmCompoundType ctype;
    PmMoleculeType mtype;

    PmAtom () {
      id = 0;
      name[0] = '\0';
      seq = 0;
      chain = '\0';
      type = PM_ATOM_ELEMENT_UNKNOWN;
      ctype = PM_COMPOUND_UNKNOWN;
      mtype = PM_MOLECULE_UNKNOWN;
      }

    PmAtom (int id, PmAtomName aname, PmAtomElementType type, char chain, int seq, 
            PmVector3 pos, PmCompoundType ctype, PmMoleculeType mtype) : 
            id(id), seq(seq), chain(chain),  pos(pos), type(type), ctype(ctype),
            mtype(mtype) {
      strcpy (name, aname);
      }

    void setData (int id, PmAtomName name, PmAtomElementType type, char chain, 
                  int seq, PmCompoundType ctype, PmMoleculeType mtype, PmVector3 pos); 

    int compare (PmAtomName name) { return (!strcmp(this->name,name)); }

    void getColor(PmVector3& color);
    float getMass();
    float getRadius();
    static PmAtomElementType convType(char *s);
    static void convType (PmAtomElementType type, char s[3]);
  };


//===== atom filter ===== //

class PM_EXPORT PmAtomFilter {
  public:
     PmAtomFilter() { mainchain = false; sidechain = false; exclude = false; 
                      peptide = false; sidechain_group = false; 
                      last_nca_atoms = false; 
                      sidechain_ca = false; sugar_phosphate = false;}
     vector<string> names;
     bool mainchain;
     bool sidechain;
     bool sidechain_group;
     bool sidechain_ca;
     bool peptide;
     bool sugar_phosphate;
     bool last_nca_atoms;
     bool exclude;
     bool hasName(const string name);
  };



#endif



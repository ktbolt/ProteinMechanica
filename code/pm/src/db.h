
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
 * db:                 d a t a b a s e                        *
 *============================================================*/

#ifndef _DB_PM_H_
#define _DB_PM_H_

#include "pm/pm.h"
#include "mol.h"
#include "grid.h"
#include "part.h"
#include "surf.h"
#include "trace.h"

namespace ProteinMechanica {

 typedef struct PmDbAtomSeq {
   char chain_id;
   int  start;
   } PmDbAtomSeq;


// PmDbInterface 
// -------------

 class PM_EXPORT PmDbInterface {

   public:
     ~PmDbInterface();
     virtual void getName(string& name) = 0;
     virtual void open (const string name, PmDbModeType otype, string data) = 0;
     virtual void read (const string mol_name) = 0;
     virtual void close(void) = 0;
     virtual PmMolecule *getMolecule(int model, PmMoleculeType type) = 0;
     virtual PmMolecule *getMolecule(string model, PmMoleculeType type) = 0;
     virtual PmGrid *getGrid() = 0;
     virtual PmSurface *getSurface() = 0;
     virtual PmTrace *getTrace() = 0;
     virtual void writeHeader() = 0;
     virtual void getModes(vector<PmMode>& modes) = 0;
     virtual void writeDomainAtoms(vector<PmMolecule*> domains, 
                                   vector<string> chains, string one_chain,
                                   bool no_hydrogen, 
                                   string start_res, string res_names, int model,
                                   string model_name) = 0;
     virtual void writeDomainSequence(vector<PmMolecule*> domains)=0;
     virtual void getModelNames(vector<string>& model_names)=0;
     static PmDbModeType convModeType(const string output);
  };
   

// PmDbPdbInterface
// ----------------
// pdb interface

 class PM_EXPORT PmDbPdbInterface: public PmDbInterface {

   public:
     PmDbPdbInterface(const string name);
     void open (const string name, PmDbModeType otype, string data);
     void getName(string& name) { name = this->name; }
     void read(const string mol_name);
     void close();
     PmMolecule *getMolecule(int model, PmMoleculeType type);
     PmMolecule *getMolecule(string model, PmMoleculeType type);
     PmSurface *getSurface() { return NULL; };
     PmTrace *getTrace() { return NULL; };
     PmGrid *getGrid() { return NULL; };
     void writeHeader();
     void writeDomainAtoms(vector<PmMolecule*> domains, vector<string> chains,
                           string one_chain, bool no_hydrogen,
                           string start_res, string res_names, int model,
                                   string model_name);
     void writeDomainSequence(vector<PmMolecule*> domains);
     void getModes(vector<PmMode>& modes) {};
     void getModelNames(vector<string>& model_names);

     void writeDynDomainAtoms(const int step, const bool init, 
                              vector<PmMolecule*> domains, vector<PmXform> xforms);

     void pm_DbPdbAtomProc (char *line);
     void pm_DbPdbAtomProc_new (char *line);
     void pm_DbPdbConectProc (char *line);
     void pm_DbPdbHelixProc(char *line);
     void pm_DbPdbHetAtomProc (char *line);
     void pm_DbPdbModelProc(char *line);
     void pm_DbPdbEndModelProc(char *line);
     void pm_DbPdbSeqresProc (char *line);
     void pm_DbPdbSheetProc(char *line);

    private:
      string name;
      FILE *fp;
      PmMoleculeType molecule_type;
      vector<PmMolecule*> protein_molecules;
      vector<PmMolecule*> water_molecules;
      vector<PmMolecule*> nucleotide_molecules;
      vector<PmMolecule*> misc_molecules;
      PmMolecule *protein, *water, *nucleotide, *misc;
      string molecule_name;
      int model_id;
      string model_name;
      vector<string> model_names;
      bool use_misc;
      vector<PmDbAtomSeq> atom_numbers;
      void getCurrentMolecule(PmMoleculeType type, PmMolecule **mol);
      void setAtomNumber(const char chain_id, const int num);
      void getAtomNumber(const char chain_id, int& num);
      void addAtoms(bool finished);
      bool classifyResAtoms();
   };


// PmDbGroInterface
// -----------------
// gromacs interface.

 class PM_EXPORT PmDbGroInterface: public PmDbInterface {

   public:
     PmDbGroInterface(const string name);
     void getName(string& name) { name = this->name; }
     void open(const string name, PmDbModeType otype, string data);
     void read(const string mol_name);
     void close();
     PmMolecule *getMolecule(int model, PmMoleculeType type);
     PmMolecule *getMolecule(string model, PmMoleculeType type);
     PmSurface *getSurface() { return NULL; };
     PmTrace *getTrace() { return NULL; };
     PmGrid *getGrid() { return NULL; };
     void writeHeader(){ };
     void getModes(vector<PmMode>& modes) {};
     void writeDomainAtoms(vector<PmMolecule*> domains, vector<string> chains,
                          string one_chain, bool no_hydrogen,
                          string start_res, string res_names, int model,
                          string model_name) {}; 
     void writeDomainSequence(vector<PmMolecule*> domains){};
     void getModelNames(vector<string>& model_names){};

    private:
      string name;
      FILE *fp;
   };

// PmDbCerfacsInterface
// --------------------
// cerfacs interface.

 class PM_EXPORT PmDbCerfacsInterface : public PmDbInterface {

   public:
     PmDbCerfacsInterface(const string name);
     void getName(string& name) { name = this->name; }
     void open(const string name, PmDbModeType otype, string data);
     void read(const string mol_name);
     void close();
     PmMolecule *getMolecule(int model, PmMoleculeType type){ return NULL; };
     PmMolecule *getMolecule(string model, PmMoleculeType type){ return NULL; };
     PmSurface *getSurface() { return NULL; };
     PmGrid *getGrid() { return NULL; };
     PmTrace *getTrace() { return NULL; };
     void writeHeader(){};
     void writeDomainAtoms(vector<PmMolecule*> domains, vector<string> chains,
                          string one_chain, bool no_hydrogen,
                          string start_res, string res_names, int model,
                          string model_name) {};
     void writeDomainSequence(vector<PmMolecule*> domains){};
     void getModes(vector<PmMode>& modes);
     void getModelNames(vector<string>& model_names){};

    private:
      string name;
      FILE *fp;

   };

// PmDbPmInterface
// ----------------
// pm interface

 class PM_EXPORT PmDbPmInterface: public PmDbInterface {

   public:
     PmDbPmInterface(const string name);
     void getName(string& name) { name = this->name; }
     void open (const string name, PmDbModeType mtype, string data);
     void read(const string name);
     void close();
     PmGrid *getGrid();
     PmMolecule *getMolecule(int model, PmMoleculeType type);
     PmMolecule *getMolecule(string model, PmMoleculeType type);
     PmParticle *getParticle();
     PmSurface *getSurface();
     PmTrace *getTrace();
     void writeHeader(){};
     void getModes(vector<PmMode>& modes) {};
     void writeDomainAtoms(vector<PmMolecule*> domains, vector<string> chains,
                           string one_chain, bool no_hydrogen,
                           string start_res, string res_names, int model,
                           string model_name) {}; 
     void writeDomainSequence(vector<PmMolecule*> domains){};

     void procGrid(char *line);
     void procParticle(char *line);
     void writeParticle(string name, vector<PmVector3> verts, vector<float> data);
     void readBinaryParticle(const string part_name, float scale);
     void procSurface (char *line);
     void writeSurface(string name, vector<PmVector3> verts, vector<int> polys);
     void readBinarySurf(const string name, float scale, bool has_charge);
     void procTrace(char *line);
     void getModelNames(vector<string>& model_names){};
     static void readEnergy(const string fname, vector<string>& names, 
                            vector<vector<float> >& values);

    private:
      string name;
      FILE *fp;
      PmMolecule *mol;
      PmGrid *grid;
      PmSurface *surf;
      PmParticle *particle;
      PmTrace *trace;
   };


// PmDbMrcInterface
// --------------------
// mrc interface.

 class PM_EXPORT PmDbMrcInterface : public PmDbInterface {

   public:
     PmDbMrcInterface(const string name);
     void getName(string& name) { name = this->name; }
     void open(const string name, PmDbModeType otype, string data);
     void read(const string mol_name);
     void close();
     PmMolecule *getMolecule(int model, PmMoleculeType type){ return NULL; };
     PmMolecule *getMolecule(string model, PmMoleculeType type){ return NULL; };
     PmSurface *getSurface() { return NULL; };
     void writeHeader(){ };
     PmGrid *getGrid();
     PmTrace *getTrace() { return NULL; };
     void writeDomainAtoms(vector<PmMolecule*> domains, vector<string> chains,
                           string one_chain, bool no_hydrogen,
                           string start_res, string res_names, int model,
                           string model_name) {};
     void writeDomainSequence(vector<PmMolecule*> domains){};
     void getModes(vector<PmMode>& modes){};
     void getModelNames(vector<string>& model_names){};

    private:
      string name;
      FILE *fp;
   };


// PmDbInterfaceSelect
// -------------------
// interface select

 class PM_EXPORT PmDbInterfaceSelect {

   public:
     PmDbInterfaceSelect(){};
     virtual PmDbInterface * create (const string name, const PmDbType type);
   };


}
#endif


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
//* mol:                    m o l e c u l e                    *
//*============================================================*

#ifndef _MOLECULE_PM_H_
#define _MOLECULE_PM_H_

#include "pm/pm.h"
#include "pobj.h"
#include "atoms.h"
#include "surf.h"
#include "graphics.h"
#include "pm/mth.h"
#include "rgn.h"

namespace ProteinMechanica {

typedef enum { 
  PM_SECONDARY_STRUCTURE_HELIX,
  PM_SECONDARY_STRUCTURE_LOOP,
  PM_SECONDARY_STRUCTURE_SHEET
  } PmSecondaryStructureType;

class PM_EXPORT PmLoop {
  public:
     int num;
     string id;
     char init_res_name[4];
     char init_chain_id;
     int init_seq_num;
     char term_res_name[4];
     char term_chain_id;
     int term_seq_num;
     int length;
  };

class PM_EXPORT PmHelix {
  public:
     int num;
     string id;
     char init_res_name[4];
     char init_chain_id;
     int init_seq_num;
     char term_res_name[4];
     char term_chain_id;
     int term_seq_num;
     int hclass;
     int length;
  };

class PM_EXPORT PmSheet {
  public:
     int num;
     string id;
     int num_strands;
     char init_res_name[4];
     char init_chain_id;
     int  init_seq_num;
     char init_icode;
     char term_res_name[4];
     char term_chain_id;
     int term_seq_num;
     char term_icode;
     int sense;
  };

typedef struct PmStructureList {
  string id;
  PmSecondaryStructureType type;
  int init_seq_num, term_seq_num;
  struct PmStructureList *next; 
  } PmStructureList;

typedef struct PmStructureListPtr {
  PmStructureList *structure;
  struct PmStructureListPtr *next;
  } PmStructureListPtr;

class PM_EXPORT PmTopology {
  public:
    char chain_id;
    int num_structures;
    PmStructureList *list;
    PmStructureListPtr **contact;
  };

typedef struct PmResidueContact {
  int num;
  PmResidue **list;
  } PmResidueContact;

class PM_EXPORT PmMolRegionParameters {
  public:
     PmMolRegionParameters() { 
       use_sidechains = false; 
       no_mainchain = false; 
       use_pca = false; 
       use_surface = false; 
       tolerance = 0.1;
       }
     bool use_sidechains;
     bool use_surface;
     bool use_pca;
     bool no_mainchain;
     float tolerance;
  };


// PmMolRegion
// -----------

class PM_EXPORT PmMolRegion : public PmRegion {
  public:
     PmMolRegion() { };
     void getCoords() { };
     void getData() { };
     void getIndices() { };
     void getRadii() { };
     string descriptor;
  };


// PmMolecule 
// ----------

 class PM_EXPORT PmMolecule : public PmPhysicalObj, public PmAtoms {

    public:

       PmMolecule(const string name);
       ~PmMolecule();

       void copy(const string name, PmMolecule **mol);

       void getAtomCoords(const string desc, PmAtomFilter& filter, 
                          vector<PmVector3>& coords);
       void getAtomCoords(const string desc, PmAtomFilter& filter, 
                          vector<PmVector3>& coords, vector<float>& rads);
       void getAtomNames(const string desc, PmAtomFilter& filter, vector<string>& names);
       void getAtoms(vector<PmAtom>& atoms);
       void getAtoms(PmAtomFilter& filter, vector<PmAtom>& atoms);
       bool getAtoms(vector<string>atoms_desc, PmAtom& atom1, PmAtom& atom2);
       void setAtomColorType(PmAtomColorType type) { this->datts.color_type = type; }
       void setAtomColor(PmVector3& color) { this->datts.atom_color = color; }
       void xformAtoms(PmXform& xform); 

       void getBackboneAtomNames(vector<string>& names);
       void addBond(PmAtom *a1, PmAtom *a2, int& n, vector<PmBond>& bonds);

       void check();

       PmChain* getChain (char chain_id); 
       void getChains (vector<PmChain*>& clist);
       void getChainAtoms(char chain_id, vector<PmAtom>& chain_atoms);
       void getChainAtomBonds(const char chain_id, vector<PmBond>& bonds);
       void getChainAtomColors (const char chain_id, PmVector3 **atom_colors);
       void getChainAtomCoords (char chain_id, string name, int *num_coords,
                                PmVector3 **coords);
       void getChainAtomRadii(char chain_id, string name, int *num, float **rads);
       void getChainPlaneGeometry(char chain_id, PmResidue *aux_res, int& num_planes, 
                                  PmConn **p_conn, int& num_coords, PmVector3 **coords);
       void getChainResidues(char chain_id, vector<PmResidue*>& rlist);
       void getChainResidues(char chain_id, char **str, vector<PmResidue*>& rlist);

       void checkClash(PmMolecule* dom, bool show);

       void addConnect(const int n, int conn[]);

       void getCoordinates(vector<PmVector3>& coords);
       void getDimensions(vector<float>& dims);
       void getRadii(vector<float>& rads);

       void getResContact();
       void checkDomainContact(vector<PmMolecule*>& list, const float tol, 
                               PmVector3& color, const bool show);
       void displayResContact(const float cutoff, const bool use_radii, 
                              PmAtomFilter& filter, PmVector3& color, bool show);

       bool checkDescriptor(const string desc);

       // domains //

       PmMolecule* createDomain (const string name, 
                                             const vector<string>& domain_names); 
       PmMolecule* createDomain (const string name, const string desc, 
                                             PmAtomFilter& filter);
       PmMolecule* createDomain (const string name);
       PmMolecule *joinDomains(const string name, const string dst_desc, 
                               PmMolecule *dst_dom, const string desc, 
                               const bool nterm, const vector<string>& chain_ids);

       // helices //

       void addHelix(PmHelix& helix);
       bool getHelix(const string id, PmHelix& helix);
       void getHelices(vector<PmHelix>& helices);
       void getHelixProps(vector<PmHelixProps>& props);
       void getHelixProps(const string desc, PmHelixProps& helix_props);
       void displayHelixProps(vector<PmHelixProps>& helix_props, 
                              PmGraphicsAttributes& atts);
       void getHelixResidues(const string dstr, vector<PmResidue*>& rlist);

       // hydrogen bonds //

       void getHydrogenBonds(const char chain_id1, vector<PmResidue*> rlist1, 
                             const char chain_id2, vector<PmResidue*> rlist2, 
                             vector<PmBond>& bonds);
       void displayHydrogenBonds(const PmVector3& color, const float width, 
                                 const bool show);

       // loops //

       void addLoop(PmLoop& loop);
       void getLoops(vector<PmLoop>& loops);
       bool getLoop(const string id, PmLoop& loop);
       void getLoopResidues(const string dstr, vector<PmResidue*>& rlist);

       // model //

       int getModelId() { return this->model_id;  }
       void setModelId(const int id) { this->model_id = id;  }
       void getModelName(string& name) { name = this->model_name;  }
       void setModelName(string name) { this->model_name = name;  }

       // modal data //

       void getMode(const int id, PmMode& mode);
       void setModes(vector<PmMode>& modes);
       void setModeColor(const PmVector3& color);
       void setModeNumber(const int number);
       void setModeScale(const float scale);

       void getMassProps(PmMassProperties& props);
       void getCenter(const string desc, PmAtomFilter& filter, PmVector3& center);
       void getName (string& name) { name = this->name;  }
       void setName (const string name) { this->name = name;  }

       void print();
       void printSecondaryStructure();

       void displayPca(PmPcaResults& pca, bool show);
       static void procQuery(PmQuery& query);
       void procQuery(PmQuery& query, char chain_id, string gtype);

       void defineRegion(const string name, const vector<int>& ids, 
                         const vector<float>& rads,
                         const vector<PmVector3>& coords, 
                         const vector<string>& atom_names=vector<string>());
       void defineRegion(const string name, const vector<int>& atom_ids);
       void defineRegion(const string name, PmAtomFilter& filter, const string desc,
                         PmMolRegionParameters& params);

       void definePcaRegion(const string name, PmAtomFilter& filter, const string desc,
                            PmMolRegionParameters& params);
       void defineSurfRegion(const string name, PmAtomFilter& filter, const string desc,
                             PmMolRegionParameters& params);

       void displayRegion(const string name, PmVector3 color, PmMoleculeRenderType rtype,
                          bool use_spheres); 
       void getRegion(const string name, PmRegion **rgn); 

       bool hasRegion(const string name);

       // residues //
       void addResBond(PmResidue *res, PmResidue *nres, vector<PmBond>& bonds);
       void addResPlane(PmResidue *res, PmResidue *nres, vector<PmMolPlane>& planes);
       int getNumResidues();
       void getResidues(const string dstr, vector<PmResidue*>& rlist);
       bool getResidueCoords (const string res, const string use_type, PmVector3& jpos);
       bool hasResidue(const char chain_id, const int id); 
       void getResidueBounds(vector<PmVector3>& center, vector<float>& radius);

       // sheets //
       void addSheet(PmSheet& sheet);
       void getSheets(vector<PmSheet>& sheets);
       bool getSheets(const string id, const int num, vector<PmSheet>& sheets);
       void getSheetResidues(const string dstr, vector<PmResidue*>& rlist);

       void getSurface (PmSurface **surf) { *surf = this->surface; }
       void setSurface (PmSurface *surf) { this->surface = surf; }

       // sidechains //
       void getSidechainGeometry(const string desc, vector<PmVector3>& coords, 
                                 vector<float>& radii);

       // topology //
       void buildTopology();
       void displayTopology(bool chain, bool contact, PmVector3& color, bool show);
       void getTopology(PmTopology **topo);

       PmMoleculeType getType();
       void setType (PmMoleculeType type);
       static void convMoleculeType(const string str, PmMoleculeType& type);

       // graphics //
       static void convMoleculeDisplayType(const string str, PmMoleculeDisplayType& type);
       static void convMoleculeRenderType(const string str, PmMoleculeRenderType& type);
       void setBondColor(PmVector3& color) { this->datts.bond_color = color; }
       void setBondColorType(PmAtomColorType type) { this->datts.bond_color_type = type;}
       void setColor(PmVector3& color) { this->datts.color = color; }
       void setDisplayBondAtoms(bool flag) { this->datts.bond_atoms = flag; }
       void setDisplayType(PmMoleculeDisplayType type) { this->datts.type = type; }
       void getExtent (PmExtent& extent ) { extent = this->extent; }
       void setLineWidth (float width) { this->datts.line_width = width; }
       void setMarkerSize (float msize) { this->datts.marker_size = msize; }
       void setRenderType(PmMoleculeRenderType type) { this->datts.render = type; }
       void setBondAtomRenderType(PmMoleculeRenderType type) { 
              this->datts.bond_atom_render = type; }
       void setSurfColor(PmVector3& color) { this->datts.surf_color = color; }

       void display(const bool show);
       void setXform (PmXform& xform);
       void getXform (PmXform& xform);
       void buildGeomName (char chain, const string geom, const string gname, 
                           string& name);

       void displayAtoms(const string gname, PmAtomFilter& filter, const bool show);
       void displayBackbone(const string name, const bool show);
       void displayBackbonePlanes(const string name, PmResidue* res, PmVector3& color,
                                  const bool show);
       void displayBonds (const string gname, const bool show);
       void displayModeVector(const string gname, PmAtomFilter& filter, const bool show);
       void displaySurface(const string gname, const bool show);

       void setTubeDisplay (const bool flag) { datts.tube = flag; }

    private:
       string name;
       int model_id;
       string model_name;
       PmMoleculeType type;
       PmMolecule *parent;
       string backbone_atom;
       vector<int> connect;
       vector<float> dimensions;

       int num_residues;
       vector<PmResidue> residues;

       int num_chains;
       vector<PmChain> chains;

       vector<PmHelix> helices;
       vector<PmLoop> loops;
       vector<PmSheet> sheets;
       PmTopology *topology;
       PmResidueContact *residue_contacts;

       PmSurface *surface;

       vector<PmMode> modes;
       int mode_number;
       PmVector3 mode_color;
       float mode_scale;

       PmXform xform;

       //vector<PmMolRegion> regions;

       vector<PmGraphicsGeometry*>  graphics_geometry;
       PmMoleculeDisplayAtts datts;

       void buildChains();
       void getChainAtoms (char chain_id, string name, vector<PmAtom*>& chain_atoms);
       void getChainAtoms (char chain_id, PmAtomFilter& filter, 
                           vector<PmAtom*>& chain_atoms);

       void addGraphicsGeometry(PmGraphicsGeometry *geom);
       void getGraphicsGeometry(const string name, PmGraphicsGeometry **geom);
       void getAtomsGeometry(const string name, PmGraphicsAtoms **geom);
       void getAxesGeometry(const string name, PmGraphicsAxes **geom);
       void getBackboneGeometry(const string name, PmGraphicsBackbone **geom);
       void getBondsGeometry(const string name, PmGraphicsBonds **geom);
       void getPlanesGeometry(const string name, PmGraphicsPlanes **geom);
       void getSurfaceGeometry(const string name, PmGraphicsSurface **geom);
   };

}

#endif


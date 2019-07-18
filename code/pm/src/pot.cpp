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
//* pot:            p o t e n t i a l                          *
//*============================================================*

#include "pot.h"
#include "sobj.h"
#include "body.h"

namespace ProteinMechanica {

#define ndbg_PmPotentialPoints_display

#define ndbg_PmPotential_cyl_use_spheres
#define ndbg_PmPotential_cyl_show_contact

#define ndbg_PmAxisPotential_compForces

#define ndbg_PmBendPotential_initConfig

#define ndbg_PmContactPotential_compForces
#define ndbg_PmContactPotential_points_points_forces
#define ndbg_PmContactPotential_sphere_sphere_forces

#define ndbg_PmEntropicSpringPotential_updateDisplay
#define ndbg_PmEntropicSpringPotential_compForces
#define ndbg_EntropicSpringForceWlc_compForce

#define ndbg_PmMolecularMechanicsPotential_compForces
#define ndbg_PmMolecularMechanicsPotential_compCoulombForces 

#define dbg_PmSpringPotential_compConnectivity
#define ndbg_PmSpringPotential_compForces
#define ndbg_PmSpringPotential_setForceConst
#define ndbg_PmSpringPotential_updateDisplay

#define ndbg_PmTorsionPotential_compForces

//////////////////////////////
// entropic spring classes //
////////////////////////////

class EntropicSpringForce {
  public:
    EntropicSpringForce();
    PmEntropicSpringFuncType force_function_type;
    bool initialized;
    virtual void initialize(const float length) = 0;
    virtual void compForce(const float dist, float& force, float& energy)=0;
    virtual void compEnergy(const float dist, float& energy)=0;
  };

// compute force from a graph //

class EntropicSpringForceGraph : public EntropicSpringForce {
  public:
    EntropicSpringForceGraph();
    int num_graph_data;
    vector<float> graph_data;
    float length;
    void compForce(const float extension, float& force, float& energy);
    void compEnergy(const float dist, float& energy);
    void initialize(const float length);
  };

// compute force using a worm-like chain formulation //

class EntropicSpringForceWlc: public EntropicSpringForce {
  public:
    EntropicSpringForceWlc();
    float persistence_length, contour_length;
    float length;
    void compEnergy(const float dist, float& energy);
    void compForce(const float dist, float& force, float& energy);
    void initialize(const float length);
  };

//////////////////////////////////////////
// parameters for Lennard-Jones energy //
////////////////////////////////////////

typedef struct LjParams {
  string atom_name;
  float sig, eps;
  } LjParams;

// LJ params: atom name, sigma (nm), eps (kJ/mol) //

LjParams LjData[] = { 
  {"C", 3.39967e-01, 3.59824e-01},
  {"N", 3.25000e-01, 7.11280e-01},
  {"O", 2.95992e-01, 8.78640e-01},
  {"P", 3.74177e-01, 8.36800e-01},
  {"S", 3.56359e-01, 1.04600e+00},
  {"",  0.0, 0.0}};

////////////////////////////////////////////////////////////////
//                    p o t e n t i al                       //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotential::PmPotential() 
  {
  this->init();
  }

PmPotential::PmPotential (const string name) 
  { 
  //fprintf (stderr, ">>> PmPotential: ctor [%s] \n", name.c_str());
  this->name = name;
  this->init();
  }

void
PmPotential::init()
  {
  geometry1 = NULL;
  geometry2 = NULL;
  sobj1 = NULL;
  sobj2 = NULL;
  verbose = false;
  current_energy = 0.0;
  report_energy = true;
  current_strain = 0.0;
  max_strain = 0.0;
  initialized = false;
  color.set(1,1,1);
  line_width = 1.0;
  activated = true;
  }

//*============================================================*
//*==========              convType                  ==========*
//*============================================================*
// convert a string potential type to a symbol.

void 
PmPotential::convType(const string str, PmPotentialType& type)
  {
  type = PM_POTENTIAL_UNKNOWN;

  if (str == "axis") {
    type = PM_POTENTIAL_AXIS;
    }

  else if (str == "bend") {
    type = PM_POTENTIAL_BEND;
    }

  else if (str == "contact") {
    type = PM_POTENTIAL_CONTACT;
    }

  else if (str == "entropic-spring") {
    type = PM_POTENTIAL_ENTROPIC_SPRING;
    }

  else if (str == "spring") {
    type = PM_POTENTIAL_SPRING;
    }

  else if (str == "torsion") {
    type = PM_POTENTIAL_TORSION;
    }

  else if (str == "molecular-mechanics") {
    type = PM_POTENTIAL_MOLECULAR_MECHANICS;
    }
  }

void
PmPotential::convType(const PmPotentialType type, string& str)
  {
  switch (type) {
    case PM_POTENTIAL_AXIS:
      str = "axis";
    break;

    case PM_POTENTIAL_BEND:
      str = "bend";
    break;

    case PM_POTENTIAL_CONTACT:
      str = "contact";
    break;

    case PM_POTENTIAL_MOLECULAR_MECHANICS:
      str = "molecular-mechanics";
    break;

    case PM_POTENTIAL_ENTROPIC_SPRING:
      str = "entropic-spring";
    break;

    case PM_POTENTIAL_SPRING:
      str = "spring";
    break;

    case PM_POTENTIAL_TORSION:
      str = "torsion";
    break;

    default:
      str = "unknown";
    }
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a potential object of a particular type.

PmPotential*
PmPotential::create(const string name, const PmPotentialType type, 
                    PmPotentialGeom *geom1,  PmPotentialGeom *geom2)
  {
  switch (type) {
    case PM_POTENTIAL_AXIS:
      {
      PmAxisPotential *axpot = new PmAxisPotential(name);
      axpot->type = type;
      axpot->geometry1 = geom1;
      axpot->geometry2 = geom2;
      geom1->getType(&axpot->geometry_type1);
      geom2->getType(&axpot->geometry_type2);
      PmPotentialLines *lgeom = dynamic_cast<PmPotentialLines*>(geom1);
      geom1->getSimulationObject(&axpot->sobj1);
      geom2->getSimulationObject(&axpot->sobj2);
      return axpot;
      }
    break;

    case PM_POTENTIAL_BEND:
      {
      PmBendPotential *bpot = new PmBendPotential(name);
      bpot->type = type;
      bpot->geometry1 = geom1;
      bpot->geometry2 = geom2;
      geom1->getType(&bpot->geometry_type1);
      geom2->getType(&bpot->geometry_type2);
      geom1->getSimulationObject(&bpot->sobj1);
      geom2->getSimulationObject(&bpot->sobj2);
      return bpot;
      }
    break;

    case PM_POTENTIAL_CONTACT:
      {
      PmContactPotential *cpot = new PmContactPotential(name);
      cpot->type = type;
      cpot->geometry1 = geom1;
      cpot->geometry2 = geom2;
      geom1->getSimulationObject(&cpot->sobj1);
      geom2->getSimulationObject(&cpot->sobj2);
      geom1->getType(&cpot->geometry_type1);
      geom2->getType(&cpot->geometry_type2);
      return cpot;
      }
    break;

    case PM_POTENTIAL_MOLECULAR_MECHANICS:
      {
      PmMolecularMechanicsPotential *mpot = new PmMolecularMechanicsPotential(name);
      mpot->type = type;
      //mpot->sobj = sobj;
      //mpot->region = rgn;
      mpot->geometry1 = geom1;
      mpot->geometry2 = geom2;
      geom1->getType(&mpot->geometry_type1);
      geom2->getType(&mpot->geometry_type2);
      geom1->getSimulationObject(&mpot->sobj1);
      geom2->getSimulationObject(&mpot->sobj2);
      return mpot;
      }
    break;

    case PM_POTENTIAL_ENTROPIC_SPRING:
      {
      PmEntropicSpringPotential *epot = new PmEntropicSpringPotential(name);
      epot->type = type;
      //epot->sobj = sobj;
      epot->geometry1 = geom1;
      epot->geometry2 = geom2;
      geom1->getSimulationObject(&epot->sobj1);
      geom2->getSimulationObject(&epot->sobj2);
      geom1->getType(&epot->geometry_type1);
      geom2->getType(&epot->geometry_type2);
      return epot;
      }
    break;

    case PM_POTENTIAL_SPRING:
      {
      PmSpringPotential *spot = new PmSpringPotential(name);
      spot->type = type;
      spot->geometry1 = geom1;
      spot->geometry2 = geom2;
      geom1->getSimulationObject(&spot->sobj1);
      geom2->getSimulationObject(&spot->sobj2);
      geom1->getType(&spot->geometry_type1);
      geom2->getType(&spot->geometry_type2);
      return spot;
      }
    break;

    case PM_POTENTIAL_TORSION:
      {
      PmTorsionPotential *tpot = new PmTorsionPotential(name);
      tpot->type = type;
      //tpot->sobj = sobj;
      tpot->geometry1 = geom1;
      tpot->geometry2 = geom2;
      geom1->getSimulationObject(&tpot->sobj1);
      geom2->getSimulationObject(&tpot->sobj2);
      geom1->getType(&tpot->geometry_type1);
      geom2->getType(&tpot->geometry_type2);
      return tpot;
      }
    break;

    default:
      return NULL;
    }
  }

//*============================================================*
//*==========              getEnergy                 ==========*
//*============================================================*

void
PmPotential::getEnergy(const float time, const bool active, float& val)
  {
  /*
  fprintf (stderr, ">>>>>> PmPotential::getEnergy active=%d report_energy=%d \n",
           active, report_energy); 
  */

  if (!active || !report_energy) {
    val = 0.0;
    }
  else {
    val = current_energy;
    }
  }

void
PmPotential::getEnergy(float& val) {
  val = current_energy;
  }

//*============================================================*
//*==========              getReportEnergy           ==========*
//*============================================================*

void
PmPotential::setReportEnergy(const bool flag) {
  report_energy = flag;
  }

void 
PmPotential::getReportEnergy(bool& flag) {
  flag = report_energy;
  }

//*============================================================*
//*==========              procQuery                 ==========*
//*============================================================*

void 
PmPotential::procQuery(PmQuery& query)
  {
  string pot_name, stype;
  unsigned int n, n1, n2, n3;
  vector<PmPotential*> pot_list;
  PmPotential *pot;
  PmPotentialGeom *pot_geom;
  float energy;
  PmPotentialType type;

  n = query.name.length();
  n1 = query.name.find("[");
  n2 = query.name.find("]");

  if ((n1 == string::npos) || (n2 == string::npos)) {
    return;
    }

  pot_name.assign(query.name, n1+1, n2-n1-1);
  //fprintf (stderr, ">>> pot_name=\"%s\" \n", pot_name.c_str());

  pmSystem.getPotential (pot_name, &pot);

  if (!pot) {
    pmSystem.getPotentialGeometry(pot_name, &pot_geom);

    if (!pot_geom) {
      return;
      }

    //fprintf (stderr, "    >>> potential geometry \n");
    pot_geom->getPotentials(pot_list);
    }
  else {
    pot_list.push_back(pot);
    }

  for (unsigned int i = 0; i < pot_list.size(); i++) {
    pot = pot_list[i];
    pot->getEnergy(energy);
    type = pot->getType();
    PmPotential::convType(type, stype);
    fprintf (stderr, "    >>> %s potential \n", stype.c_str());
    fprintf (stderr, "    >>> total energy=%g \n", energy); 
    pot->queryState(query);
    fprintf (stderr, " \n"); 
    }
  }

//*============================================================*
//*==========              getStrain                 ==========*
//*============================================================*

void
PmPotential::getStrain (const float time, const bool active, float& val)
  {
  if (!active) {
    val = 0.0;
    }
  else {
    val = current_strain;
    }
  }

//*============================================================*
//*==========              getMaxStrain              ==========*
//*============================================================*

void 
PmPotential::getMaxStrain(float& val) {
  val = max_strain;
  }

void 
PmPotential::setMaxStrain(const float val) {
  max_strain = val / 100.0;
  }

//*============================================================*
//*==========              printGeometry             ==========*
//*============================================================*

void
PmPotential::printGeometry()
  {
  if (this->geometry1) { 
    this->geometry1->print();
    }

  if (this->geometry2) { 
    this->geometry2->print();
    }
  }

//*============================================================*
//*==========              displayGeometry           ==========*
//*============================================================*
// set a flag to show potential geometry and connectivity for 
// spring potentials.

void
PmPotential::displayGeometry(const bool show)
  {

  if (!show) {
    return;
    }

  if (type == PM_POTENTIAL_ENTROPIC_SPRING) {
    PmEntropicSpringPotential *epot = dynamic_cast<PmEntropicSpringPotential*>(this);
    epot->setDisplayConn(show);
    }

  else if (type == PM_POTENTIAL_SPRING) {
    PmSpringPotential *pot = dynamic_cast<PmSpringPotential*>(this);
    pot->setDisplayConn(show);
    }

  else if (type == PM_POTENTIAL_BEND) {
    PmBendPotential *pot = dynamic_cast<PmBendPotential*>(this);
    pot->setDisplayConn(show);
    }

  else if (type == PM_POTENTIAL_TORSION) {
    PmTorsionPotential *pot = dynamic_cast<PmTorsionPotential*>(this);
    pot->setDisplayConn(show);
    }
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a potenail geometry.                

void
PmPotential::setXform(PmXform& xform)
  {
  //fprintf (stderr, ">>>>>> PmPotential::setXform\n");
  this->geometry1->setXform(xform);
  this->geometry2->setXform(xform);
  }

void
PmPotential::getXform(PmXform& xform)
  {
  this->geometry1->getXform(xform);
  }

//*============================================================*
//*==========              initSpatialDecomp         ==========*
//*============================================================*
// initialize a spatial decomposition used by some potentials.

void
PmPotential::initSpatialDecomp(PmPotentialPoints *points, 
                               PmSpatialDecomp **p_spatial_decomp, float width)
  {
  //fprintf (stderr, "\n>>>>>> PmPotential: init spatial decomp \n");
  vector<PmVector3> coords;
  PmVector3 origin, u, v, w, pt; 
  PmPcaResults pca;
  PmExtent proj, pext;
  float wd[3], du, dv, dw, d, r;
  float lu, lv, lw; 

  for (int i = 0; i < points->number; i++) {
    coords.push_back(points->coordinates[i]);
    }

  //fprintf (stderr, ">>> num pts = %d \n", points->number); 
  pm_MathPrincipalComponentsProj(coords, points->pca, proj, wd);
  pca = points->pca;

  // ------------------------------//
  // create spatial decomposition  //
  // ------------------------------//

  PmSpatialDecomp *spatial_decomp = new PmSpatialDecomp(name);

  // set origin //

  d = 0.5 + width;
  u = pca.axis1;
  v = pca.axis2;
  w = pca.axis3;

  lu = proj.max[0] - proj.min[0] + 2*d;
  lv = proj.max[1] - proj.min[1] + 2*d;
  lw = proj.max[2] - proj.min[2] + 2*d;

  origin = pca.com - 0.5*lu*u - 0.5*lv*v - 0.5*lw*w;
  spatial_decomp->setCoordSys(origin, u, v, w); 

  // set axes lenth //

  spatial_decomp->setAxesLength(lu, lv, lw);

  // set cell size //

  du = 0.5; dv = 0.5; dw = 0.5;
  spatial_decomp->setCellSize(du, dv, dw);

  // add points //

  for (unsigned int i = 0; i < coords.size(); i++) {
    pt = coords[i];
    r = points->radii[i];
    spatial_decomp->addPoint(i, pt, r);
    }

  //spatial_decomp->printInfo();
  //spatial_decomp->display();

  // test //

  #define ntest_PmMolecularMechanicsPotential_initPotential
  #ifdef test_PmMolecularMechanicsPotential_initPotential

  fprintf (stderr, "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  fprintf (stderr, "\n>>>>>>> test spatial decomp: \n");

  vector<int> ids;
  int id, n;

  //pt = origin + 0.5*lu*u + 0.5*lv*v + 0.5*lw*w;
  pt[0] = 4.64212, pt[1] = 5.15401, pt[2] = 5.41185; 
  pt[0] = 4.26566, pt[1] = 5.8598,  pt[2] = 4.45516;
  pt[0] = -2.0998, pt[1] = 18.4275, pt[2] = 24.7632;  
  r = 0.2;
  
  string geom_name;
  PmGraphicsAttributes atts, atts1;
  PmVector3 color;
  PmGraphicsSphere *sgeom1;
  PmVector3 verts[1000];
  float rads[1000];

  geom_name = "mm_show_query";
  verts[0] = pt;
  rads[0] = r;
  sgeom1 = new PmGraphicsSphere(geom_name, 1, verts, rads);
  color.set(0, 1, 0);
  atts.setColor(color);
  atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
  sgeom1->setAttributes(atts);
  sgeom1->display();

  spatial_decomp->getPoints(pt, r, ids);
  n = ids.size();
  fprintf (stderr, ">>> num points in cell = %d  | ", n);

  if (n) {
    for (int i = 0; i < n; i++) {
      id = ids[i];
      pt = coords[id];
      r = points->radii[id];
      verts[i] = pt;
      rads[i] = r;
      fprintf (stderr, " %d ", id);
      }

    fprintf (stderr, "\n");

    geom_name = "mm_show_cells_contents";
    sgeom1 = new PmGraphicsSphere(geom_name, n, verts, rads);
    color.set(1, 0, 0);
    atts.setColor(color);
    atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
    sgeom1->setAttributes(atts);
    sgeom1->display();
    }
  #endif

  *p_spatial_decomp = spatial_decomp;
  }

//*============================================================*
//*==========              setColor                  ==========*
//*============================================================*

void
PmPotential::setColor(ProteinMechanica::PmVector3& color) {
  this->color = color;
  }

//*============================================================*
//*==========              setLineWidth              ==========*
//*============================================================*

void
PmPotential::setLineWidth(float width) {
  this->line_width = width;
  }


////////////////////////////////////////////////////////////////
//              p o t e n t i al   g e o m e t r y           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotentialGeom::PmPotentialGeom()
  {
  //fprintf (stderr, ">>>>>> PmPotentialGeom::ctor\n");
  graphics_geometry = NULL;
  graphics_geometry1 = NULL;
  contact_geometry = NULL;
  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  line_width = 1.0;
  display_spheres = false;
  visible = true;
  type = PM_POTENTIAL_GEOM_UNKNOWN;
  properties.charge = 0.0;
  properties.has_charge = false;
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a potential geometry object of a particular type.
// initial geometric properties are obtained from performing 
// a principal component analysis.

PmPotentialGeom* 
PmPotentialGeom::create(const string name, const PmPotentialGeomType type, 
                        PmSimulationObj *sobj, const string rgn, 
                        PmPotentialParameters& params)
  {
  //fprintf (stderr, ">>>>>> PmPotentialGeom::create  rgn[%s] \n", rgn.c_str());
  vector<PmVector3> coords;
  vector<float> rads;
  PmPcaResults pca;
  PmExtent proj;
  float wd[3];
  float d, max_radius = 0.0;
  PmVector3 r;

  //===== get the coordinates for the region =====//

  if (sobj) {
    sobj->getRegionCoords(rgn, coords);
    sobj->getRegionRadii(rgn, rads);
    }

  //===== compute principal components =====//

  //fprintf (stderr, "   >>> num coords = %d \n", coords.size());

  if (coords.size() == 1) {
    wd[0] = rads[0];
    wd[1] = rads[0];
    wd[2] = rads[0];

    pca.s1 = 1.0;
    pca.s2 = 1.0;
    pca.s3 = 1.0;
    pca.com = coords[0];

    pca.axis1[0] = 1.0; pca.axis1[1] = 0.0; pca.axis1[2] = 0.0;
    pca.axis2[0] = 0.0; pca.axis2[1] = 1.0; pca.axis2[2] = 0.0;
    pca.axis3[0] = 0.0; pca.axis3[1] = 0.0; pca.axis3[2] = 1.0;
    }
  else if (coords.size() != 0) {
    pm_MathPrincipalComponents(coords, pca);
    pm_MathPrincipalComponentsProj(coords, pca, proj, wd);

    // add a bit more to include atom radius? //

    /*
    wd[0] += 0.2; 
    wd[1] += 0.2; 
    wd[2] += 0.2; 
    */
    }

  //===== create a geometry of a particular type =====//

  switch (type) {
    case PM_POTENTIAL_GEOM_CYLINDER:
      {
      PmPotentialCylinder *cgeom = new PmPotentialCylinder;
      cgeom->name = name;
      cgeom->rgn_name = rgn;
      cgeom->type = type;
      cgeom->simulation_obj = sobj;
      cgeom->pca = pca;
      cgeom->center = pca.com;
      cgeom->curr_center = pca.com;

      // set cylinder axis scales //

      cgeom->length = 2.0*wd[0];
      cgeom->max_rad = wd[0];
      cgeom->rads[0] = wd[0];
      cgeom->rads[1] = wd[1];
      cgeom->rads[2] = wd[2];

      // set cylinder axis //

      cgeom->u = pca.axis1;
      cgeom->curr_u = pca.axis1;
      cgeom->v = pca.axis2;
      cgeom->curr_v = pca.axis2;
      cgeom->w = pca.axis3;
      cgeom->curr_w = pca.axis3;
      return cgeom;
      }
    break;

    case PM_POTENTIAL_GEOM_ELLIPSOID:
      {
      PmPotentialEllipsoid *egeom = new PmPotentialEllipsoid;
      egeom->name = name;
      egeom->rgn_name = rgn;
      egeom->type = type;
      egeom->simulation_obj = sobj;
      egeom->pca = pca;
      egeom->center = pca.com;
      egeom->curr_center = pca.com;

      // set ellipsoid axis scales //

      egeom->max_rad = wd[0];
      egeom->rads[0] = wd[0];
      egeom->rads[1] = wd[1];
      egeom->rads[2] = wd[2];

      // set ellipsoid axis //

      egeom->u = pca.axis1;
      egeom->curr_u = pca.axis1;
      egeom->v = pca.axis2;
      egeom->curr_v = pca.axis2;
      egeom->w = pca.axis3;
      egeom->curr_w = pca.axis3;
      /*
      fprintf (stderr, "    >>> ellipsoid: axis scales (%g %g %g) \n",  wd[0], wd[1], 
               wd[2]);
      */
      return egeom;
      }
    break;

    case PM_POTENTIAL_GEOM_LINES:
      {
      PmPotentialLines *lgeom = new PmPotentialLines;
      lgeom->name = name;
      lgeom->rgn_name = rgn;
      lgeom->type = type;
      lgeom->simulation_obj = sobj;
      lgeom->origin = pca.com;
      lgeom->curr_origin = pca.com;
      lgeom->pca = pca;
      lgeom->pca_axes[0] = pca.axis1;
      lgeom->pca_axes[1] = pca.axis2;
      lgeom->pca_axes[2] = pca.axis3;
      lgeom->scales[0] = wd[0];
      lgeom->scales[1] = wd[1];
      lgeom->scales[2] = wd[2];
      lgeom->coordinates = new PmVector3[6];
      lgeom->curr_coordinates = new PmVector3[6];
      return lgeom;
      }
    break;

    case PM_POTENTIAL_GEOM_PLANE:
      {
      PmVector3 center, normal, scales;
      PmPotentialPlane *plgeom = new PmPotentialPlane;
      plgeom->name = name;
      plgeom->rgn_name = rgn;
      plgeom->type = type;
      plgeom->simulation_obj = sobj;
      plgeom->coordinates = new PmVector3[4];

      if (params.getPlaneData(center, normal)) { 
        params.getPlaneScales(scales); 
        plgeom->setPoint(center); 
        plgeom->setNormal(normal); 
        plgeom->setScale(scales[0], scales[1]); 
        }

      return plgeom;
      }
    break;

    case PM_POTENTIAL_GEOM_POINTS:
      {
      PmPotentialPoints *pgeom = new PmPotentialPoints;
      pgeom->name = name;
      pgeom->rgn_name = rgn;
      pgeom->type = type;
      pgeom->simulation_obj = sobj;
      pgeom->center = pca.com;
      pgeom->curr_center = pca.com;
      pgeom->pca = pca;
      //pgeom->radius = wd[0];
      pgeom->rads[0] = wd[0];
      pgeom->rads[1] = wd[1];
      pgeom->rads[2] = wd[2];
      pgeom->number = coords.size();
      pgeom->coordinates = new PmVector3[coords.size()];
      pgeom->curr_coordinates = new PmVector3[coords.size()];
      pgeom->radii = new float[coords.size()];
      bool use_rads = rads.size() != 0;

      for (unsigned int i = 0; i < coords.size(); i++) {
        pgeom->coordinates[i] = coords[i];
        pgeom->curr_coordinates[i] = coords[i];

        if (use_rads) {
          pgeom->radii[i] = rads[i]; 
          //fprintf (stderr, " rads[i] %g \n", rads[i]);
          }
        else {
          pgeom->radii[i] = 0.1; 
          }

        r = coords[i] - pca.com;
        d = r.length() + pgeom->radii[i];
        if (d > max_radius) max_radius = d;
        }

      pgeom->radius = max_radius;
      return pgeom;
      }
    break;

    case PM_POTENTIAL_GEOM_SPHERE:
      {
      PmPotentialSphere *sgeom = new PmPotentialSphere;
      sgeom->name = name;
      sgeom->rgn_name = rgn;
      sgeom->type = type;
      sgeom->simulation_obj = sobj;
      sgeom->pca = pca;
      sgeom->radius = wd[0];
      sgeom->center = pca.com;
      sgeom->curr_center = pca.com;
      //fprintf (stderr, "    >>> sphere: initial radius[%g] \n", wd[0]); 
      //fprintf (stderr, "        wd = %f %f %f              \n", wd[0], wd[1], wd[2]); 
      return sgeom;
      }
    break;

    default:
      return NULL;
    }
  }

//*============================================================*
//*==========              convType                  ==========*
//*============================================================*
// convert a string potential geometry type to a symbol.

void
PmPotentialGeom::convType(const string str, PmPotentialGeomType& type)
  {
  type = PM_POTENTIAL_GEOM_UNKNOWN;

  if (str == "cylinder") {
    type = PM_POTENTIAL_GEOM_CYLINDER;
    }
  else if (str == "ellipsoid") {
    type = PM_POTENTIAL_GEOM_ELLIPSOID;
    }
  else if (str == "lines") {
    type = PM_POTENTIAL_GEOM_LINES;
    }
  else if (str == "plane") {
    type = PM_POTENTIAL_GEOM_PLANE;
    }
  else if (str == "points") {
    type = PM_POTENTIAL_GEOM_POINTS;
    }
  else if (str == "sphere") {
    type = PM_POTENTIAL_GEOM_SPHERE;
    }
  }

//*============================================================*
//*==========               getType                  ==========*
//*============================================================*

void 
PmPotentialGeom::getType(PmPotentialGeomType *type) {
  *type = this->type;
  }

//*============================================================*
//*==========          getSimulationObject           ==========*
//*============================================================*

void 
PmPotentialGeom::getSimulationObject(PmSimulationObj **sobj)
  {
  *sobj = this->simulation_obj;
  }

//*============================================================*
//*==========          addPotential                  ==========*
//*============================================================*

void 
PmPotentialGeom::addPotential(PmPotential *pot) {
  potential_list.push_back(pot);
  }

//*============================================================*
//*==========          getPotentials                 ==========*
//*============================================================*

void 
PmPotentialGeom::getPotentials(vector<PmPotential*>& list)
  {
  list = potential_list;
  }

//*============================================================*
//*==========               get/setCharge            ==========*
//*============================================================*

void 
PmPotentialGeom::setCharge(const float charge)
  {
  this->properties.charge = charge;
  this->properties.has_charge = true;
  }

void 
PmPotentialGeom::getCharge(float& charge) {
  charge = this->properties.charge;
  }

//*============================================================*
//*==========          set graphics attributes       ==========*
//*============================================================*

void 
PmPotentialGeom::setColor(PmVector3& color) {
  this->color = color;
  }

void 
PmPotentialGeom::setDisplayType(PmGeometryDisplayType type) {
  this->display_type = type;
  }

void 
PmPotentialGeom::setDisplaySpheres(const bool flag) {
  this->display_spheres = flag;
  }

void 
PmPotentialGeom::setLineWidth(float width) {
  this->line_width = width;
  }

void 
PmPotentialGeom::setShadingType(PmGeometryShadingType type) {
  this->shading_type = type;
  }

//*============================================================*
//*==========                getName                 ==========*
//*============================================================*

void 
PmPotentialGeom::getName(string& name) { 
  name = this->name;  
  }

//*============================================================*
//*==========                hasName                 ==========*
//*============================================================*

bool 
PmPotentialGeom::hasName(const string& name) { 
  return name == this->name;  
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// get the xform for a potential geometry.

void
PmPotentialGeom::getXform(PmXform& xform)
  {
  xform = potential_xform;
  }

////////////////////////////////////////////////////////////////
//              s p h e r e        g e o m e t r y           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotentialSphere::PmPotentialSphere()
  {
  //fprintf (stderr, ">>>>>> PmPotentialSphere::ctor\n");
  radius = 0.0;
  mass = 0.0;
  }

void
PmPotentialSphere::print()
  {
  fprintf (stderr, "\n");
  fprintf (stderr, ">>> geometry = sphere \n");
  fprintf (stderr, ">>> center   = %f %f %f \n", center[0], center[1], center[2]);
  fprintf (stderr, ">>> radius   = %f \n", radius);
  }

//*============================================================*
//*==========              setSphereRadius           ==========*
//*============================================================*
// set the radius of a sphere potential. this will override the
// radius initially set from pca data.

void
PmPotentialSphere::setSphereRadius(const float radius) {
  this->radius = radius;
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*

void
PmPotentialSphere::display(const bool show)
  {
  //fprintf (stderr, ">>>>>> PmPotentialSphere::display \n");
  string geom_name = "potential";
  float r = radius;
  PmVector3 pos = center;
  //fprintf (stderr, "   >>> radius[%g] \n", r);
  //fprintf (stderr, "   >>> pos   (%g %g %g) \n", pos[0], pos[1], pos[2]);

  geom_name = geom_name + '[' + name + ']' + "sphere";
  int num_verts = 1;
  PmVector3 *verts = new PmVector3[1];
  verts[0] = center;
  float *rads = new float[1];
  rads[0] = r;

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setDisplayType(display_type);
  atts.setLineWidth(this->line_width); 

  //fprintf (stderr, "   >>> geom name \"%s\" \n", geom_name.c_str()); 
  PmGraphicsSphere *geom = new PmGraphicsSphere(geom_name, num_verts, verts, rads);
  geom->setAttributes(atts);
  geom->display();
  graphics_geometry = geom; 

  /*
  PmVector3 *cverts = new PmVector3[1];
  cverts[0] = center;
  PmGraphicsPoint *pgeom = new PmGraphicsPoint(geom_name, 1, cverts);
  pgeom->display();
  contact_geometry = pgeom;
  */
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of a points potential geometry.

void
PmPotentialSphere::updateDisplay()
  {
  //fprintf (stderr, ">>>>>> PmPotentialSphere::updateDisplay[%s]\n", name.c_str());
  if (!graphics_geometry) {
    return;
    }

  PmGraphicsSphere *geom = dynamic_cast<PmGraphicsSphere*>(graphics_geometry);
  geom->setXform (potential_xform);
  geom->display();
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a sphere potential geometry.

void
PmPotentialSphere::setXform(PmXform& xform)
  {
  //fprintf (stderr, ">>>>>> PmPotentialSphere::setXform  [%s]\n", name.c_str());
  if (!graphics_geometry) {
    return;
    }

  // transform center //

  PmMatrix3x3 mat = xform.matrix;
  PmVector3 xpt = mat*(center-xform.center) + xform.translation + xform.center;
  curr_center = xpt;

  // set potential geometry xform //

  potential_xform = xform;
  }

////////////////////////////////////////////////////////////////
//              p l a n e          g e o m e t r y           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotentialPlane::PmPotentialPlane()
  {
  normal.set(1,0,0);
  center.set(0,0,0);
  scale_x = scale_y = 1.0;
  color.set(1,1,1);
  }

void
PmPotentialPlane::print()
  {
  fprintf (stderr, "\n");
  fprintf (stderr, ">>> geometry = plane \n");
  fprintf (stderr, ">>> center   = %f %f %f \n", center[0], center[1], center[2]);
  }

//*============================================================*
//*==========              setNormal                 ==========*
//*============================================================*

void 
PmPotentialPlane::setNormal(const PmVector3& vec)
  {
  normal = vec;
  normal.normalize();
  }

//*============================================================*
//*==========              setPoint                  ==========*
//*============================================================*

void 
PmPotentialPlane::setPoint(const PmVector3& pt)
  {
  center = pt;
  }

//*============================================================*
//*==========              setScale                  ==========*
//*============================================================*

void
PmPotentialPlane::setScale(float sx, float sy) 
  {
  scale_x = sx;
  scale_y = sy;
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*

void
PmPotentialPlane::display(const bool show)
  {
  string geom_name = "potential";
  PmVector3 v1, v2; 
  PmVector3 pos = center;

  pm_MathBasisCompute (normal, v1, v2);

  geom_name = geom_name + '[' + name + ']' + "plane";
  int num_verts = 4;
  PmVector3 *verts = new PmVector3[4];
  coordinates[0] = verts[0] = center - scale_x*v1 - scale_y*v2;
  coordinates[1] = verts[1] = center + scale_x*v1 - scale_y*v2;
  coordinates[2] = verts[2] = center + scale_x*v1 + scale_y*v2;
  coordinates[3] = verts[3] = center - scale_x*v1 + scale_y*v2;

  PmConn *pconn = new PmConn(5);
  (*pconn)[0] = 4;
  (*pconn)[1] = 0;
  (*pconn)[2] = 1;
  (*pconn)[3] = 2;
  (*pconn)[4] = 3;

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setDisplayType(display_type);
  atts.setLineWidth(this->line_width); 

  PmGraphicsPolygon *geom = new PmGraphicsPolygon(geom_name, 1, pconn, 0, num_verts, 
                                                  verts);
  geom->setAttributes(atts);
  geom->display();
  graphics_geometry = geom; 

  PmGraphicsLine *lgeom = new PmGraphicsLine(geom_name, 4, verts);
  lgeom->setAttributes(atts);
  lgeom->display();
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of a plane potential geometry.

void
PmPotentialPlane::updateDisplay()
  {
  if (!graphics_geometry) {
    return;
    }

  /*
  PmGraphicsSphere *geom = dynamic_cast<PmGraphicsSphere*>(graphics_geometry);
  geom->setXform (potential_xform);
  geom->display();
  */
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a plane potential geometry.

void
PmPotentialPlane::setXform(PmXform& xform)
  {
  //fprintf (stderr, ">>>>>> PmPotentialPlane::setXform  [%s]\n", name.c_str());

  if (!graphics_geometry) {
    return;
    }

  // transform center //

  PmMatrix3x3 mat = xform.matrix;
  PmVector3 xpt = mat*(center-xform.center) + xform.translation + xform.center;
  curr_center = xpt;

  // set potential geometry xform //

  potential_xform = xform;
  }


////////////////////////////////////////////////////////////////
//           c y l i n d e r       g e o m e t r y           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotentialCylinder::PmPotentialCylinder()
  {
  for (int i = 0; i < 3; i++) {
    rads[i] = 0.0;
    }

  max_rad = 0.0;
  }

void
PmPotentialCylinder::print()
  {
  fprintf (stderr, ">>> geometry = cylinder \n");
  fprintf (stderr, ">>> center   = %f %f %f \n", center[0], center[1], center[2]);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*

void
PmPotentialCylinder::display(const bool show)
  {
  string geom_name = "potential";
  PmVector3 pos = center;

  geom_name = geom_name + '[' + name + ']' + "cylinder";
  int num_verts = 1;
  PmVector3 *verts = new PmVector3[1];
  verts[0] = center;

  // set ellipsoid axis data //

  PmEllipsoidData *data = new PmEllipsoidData[1];
  data[0].radius1 = rads[0]; data[0].axis1 = u;
  data[0].radius2 = rads[1]; data[0].axis2 = v;
  data[0].radius3 = rads[2]; data[0].axis3 = w;

  PmGraphicsEllipsoid *geom = new PmGraphicsEllipsoid(geom_name, num_verts, verts, data);
  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setDisplayType(display_type);
  geom->setAttributes(atts);
  geom->display();
  graphics_geometry = geom;

  #ifdef dbg_PmPotential_cyl_show_contact
  PmVector3 *cverts = new PmVector3[1];
  PmVector3 color;
  cverts[0] = center;
  PmGraphicsPoint *pgeom = new PmGraphicsPoint(geom_name, 1, cverts);
  PmGraphicsAttributes atts;
  color.set(1,0,1);
  atts.setColor(color);
  atts.setScale(0.1);
  pgeom->setAttributes(atts);
  pgeom->display();
  contact_geometry = pgeom;
  #endif

  #ifdef dbg_PmPotential_cyl_use_spheres
  PmVector3 *sverts = new PmVector3[2];
  sverts[0] = center - (rads[0]-rads[1])*u;
  sverts[1] = center + (rads[0]-rads[1])*u;
  float srads[2];
  srads[0] = rads[1];
  srads[1] = rads[1];
  PmGraphicsSphere *sgeom = new PmGraphicsSphere(geom_name, 2, sverts, srads);
  sgeom->display();
  graphics_geometry1 = sgeom;
  #else
  graphics_geometry1 = NULL;
  #endif
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of a cylinder potential geometry.

void
PmPotentialCylinder::updateDisplay()
  {
  if (!graphics_geometry) {
    return;
    }

  PmGraphicsCylinder *geom = dynamic_cast<PmGraphicsCylinder*>(graphics_geometry);

  if (geom) {
    geom->setXform (potential_xform);
    geom->display();
    }
  }

//*============================================================*
//*==========              getGeometry               ==========*
//*============================================================*
// get the current geometry for a cylinder.

void
PmPotentialCylinder::getGeometry(PmVector3& pos, PmVector3& gu, PmVector3& gv,
                                 PmVector3& gw)
  {
  pos = curr_center;
  gu = curr_u;
  gv = curr_v;
  gw = curr_w;
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a cylinder potential geometry.

void
PmPotentialCylinder::setXform(PmXform& xform)
  {
  if (!graphics_geometry) {
    return;
    }

  PmMatrix3x3 mat = xform.matrix;

  // transform center //

  curr_center = mat*(center-xform.center) + xform.translation + xform.center;

  // transform axes //

  curr_u = mat*u;
  curr_v = mat*v;
  curr_w = mat*w;

  // transform potential geometry //

  /*
  PmGraphicsEllipsoid *geom = dynamic_cast<PmGraphicsEllipsoid*>(graphics_geometry);
  geom->setXform (xform);
  geom->display();

  if (graphics_geometry1) {
    PmGraphicsSphere *sgeom = dynamic_cast<PmGraphicsSphere*>(graphics_geometry1);
    sgeom->setXform (xform);
    sgeom->display();
    }
  */

  potential_xform = xform;
  }

////////////////////////////////////////////////////////////////
//           e l l i p s o i d     g e o m e t r y           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotentialEllipsoid::PmPotentialEllipsoid()
  {
  for (int i = 0; i < 3; i++) {
    rads[i] = 0.0;
    }
  max_rad = 0.0;
  }

void
PmPotentialEllipsoid::print()
  {
  fprintf (stderr, ">>> geometry = ellipsoid \n");
  fprintf (stderr, ">>> radii    = %f %f %f \n", rads[0], rads[1], rads[2]);
  fprintf (stderr, ">>> center   = %f %f %f \n", center[0], center[1], center[2]);
  }

//*============================================================*
//*==========          setEllipsoidScale             ==========*
//*============================================================*
// scale the ellipsoid axes by the given amount.

void
PmPotentialEllipsoid::setEllipsoidScale (const float scale) {
  this->rads[0] *= scale;
  this->rads[1] *= scale;
  this->rads[2] *= scale;
  this->max_rad = this->rads[0];
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*

void
PmPotentialEllipsoid::display(const bool show)
  {
  string geom_name = "potential";
  PmVector3 pos = center;

  geom_name = geom_name + '[' + name + ']' + "ellipsoid";
  int num_verts = 1;
  PmVector3 *verts = new PmVector3[1];
  verts[0] = center;

  // set ellipsoid axis data //

  PmEllipsoidData *data = new PmEllipsoidData[1];
  data[0].radius1 = rads[0]; data[0].axis1 = u; 
  data[0].radius2 = rads[1]; data[0].axis2 = v; 
  data[0].radius3 = rads[2]; data[0].axis3 = w; 

  PmGraphicsEllipsoid *geom = new PmGraphicsEllipsoid(geom_name, num_verts, verts, data);
  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setDisplayType(display_type);
  atts.setLineWidth(this->line_width); 
  geom->setAttributes(atts);
  geom->display();
  graphics_geometry = geom;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of a points potential geometry.

void
PmPotentialEllipsoid::updateDisplay()
  {
  if (!graphics_geometry) {
    return;
    }

  PmGraphicsEllipsoid *geom = dynamic_cast<PmGraphicsEllipsoid*>(graphics_geometry);
  geom->setXform (potential_xform);
  geom->display();
  }

//*============================================================*
//*==========              getGeometry               ==========*
//*============================================================*
// get the current geometry for an ellipsoid.    

void
PmPotentialEllipsoid::getGeometry(PmVector3& pos, PmVector3& gu, PmVector3& gv, 
                                  PmVector3& gw)
  {
  pos = curr_center;
  gu = curr_u; 
  gv = curr_v; 
  gw = curr_w;
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a ellipsoid potential geometry.

void
PmPotentialEllipsoid::setXform(PmXform& xform)
  {
  if (!graphics_geometry) {
    return;
    }

  /*
  fprintf (stderr, ">>>>>> PmPotentialEllipsoid::setXform  [%s]\n", name.c_str());
  fprintf (stderr, "   >>> center (%g %g %g) \n", xform.center[0], xform.center[1], 
           xform.center[2]);
  fprintf (stderr, "   >>> angles (%g %g %g) \n", xform.angles[0], xform.angles[1], 
           xform.angles[2]);
  fprintf (stderr, "   >>> translation (%g %g %g) \n", xform.translation[0], 
           xform.translation[1], xform.translation[2]);
  */
  PmMatrix3x3 mat = xform.matrix;

  // transform center //

  curr_center = mat*(center-xform.center) + xform.translation + xform.center;

  // transform axes //

  curr_u = mat*u;
  curr_v = mat*v;
  curr_w = mat*w;

  potential_xform = xform;
  }

////////////////////////////////////////////////////////////////
//              l i n e s          g e o m e t r y           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPotentialLines::PmPotentialLines()
  {
  number = 0; 
  num_axes = 0;
  origin.set(0,0,0);
  color.set(1,1,1);
  scale = 1.0;
  scales[0] = 1.0;
  scales[1] = 1.0;
  scales[2] = 1.0;
  }

void
PmPotentialLines::print()
  {
  fprintf (stderr, "\n");
  fprintf (stderr, ">>> geometry = lines \n");
  }

//*============================================================*
//*==========              setAxes                   ==========*
//*============================================================*
// set axes.

void
PmPotentialLines::setAxes(const vector<int>& ids)
  {
  int id;

  for (int i = 0; i < ids.size(); i++) {
    id = ids[i];

    if ((id >= 1) && (id <= 3)) {
      this->axes[i] = this->pca_axes[id-1];
      this->axes_ids[i] = id-1;
      this->num_axes += 1;
      }
    }
  }

void
PmPotentialLines::setAxes(const int id, const int aid)
  {
  if ((id < 1) || (id > 3) || (aid < 1) || (aid > 3)) {
    return;
    }

  this->axes[id-1] = this->pca_axes[aid-1];

  if (this->num_axes < 3) {
    this->num_axes += 1;
    }
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*

void
PmPotentialLines::display(const bool show)
  {
  /*
  fprintf (stderr, ">>>>>> PmPotentialLines::display name = %s\n", name.c_str());
  fprintf (stderr, "   >>> origin = %f %f %f \n", origin[0], origin[1], origin[2]);   
  fprintf (stderr, "   >>> scales = %f %f %f \n", scales[0], scales[1], scales[2]);   
  */
  string geom_name = "potential";
  PmVector3 v1, v2; 
  PmVector3 pos = this->origin;
  float s1, s2, s3;

  geom_name = geom_name + '[' + name + ']' + "lines";
  int num_verts = 4;
  PmVector3 *verts = new PmVector3[6];
  s1 = this->scales[0]; 
  s2 = this->scales[1]; 
  s3 = this->scales[2]; 
  number = 6;
  coordinates[0] = verts[0] = this->origin;
  coordinates[1] = verts[1] = this->origin + this->scale*s1*this->pca.axis1;
  coordinates[2] = verts[2] = this->origin;
  coordinates[3] = verts[3] = this->origin + this->scale*s2*this->pca.axis2;
  coordinates[4] = verts[4] = this->origin;
  coordinates[5] = verts[5] = this->origin + this->scale*s3*this->pca.axis3;

  PmGraphicsAttributes atts;
  atts.setColor(color);
  atts.setDisplayType(display_type);
  atts.setDisjoint(true);
  atts.setLineWidth(this->line_width); 

  PmGraphicsLine *geom = new PmGraphicsLine(geom_name, 6, verts);
  geom->setAttributes(atts);
  geom->display();
  graphics_geometry = geom; 
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of a lines potential geometry.

void
PmPotentialLines::updateDisplay()
  {
  if (!graphics_geometry) {
    return;
    }

  //PmGraphicsLine *geom = dynamic_cast<PmGraphicsLine*>(graphics_geometry);
  //geom->setXform (potential_xform);
  graphics_geometry->setXform (potential_xform);
  graphics_geometry->display();
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a lines potential geometry.

void
PmPotentialLines::setXform(PmXform& xform)
  {
  //fprintf (stderr, ">>>>>> PmPotentialLines::setXform  [%s]\n", name.c_str());

  if (!graphics_geometry) {
    return;
    }

  // transform center //

  PmMatrix3x3 mat = xform.matrix;
  PmVector3 xpt = mat*(origin-xform.center) + xform.translation + xform.center;
  curr_origin = xpt;

  /*
  fprintf (stderr, "   >>> mat(0,*) = %f %f %f \n", mat(0,0), mat(0,1), mat(0,2)); 
  fprintf (stderr, "   >>> mat(1,*) = %f %f %f \n", mat(1,0), mat(1,1), mat(1,2)); 
  fprintf (stderr, "   >>> mat(2,*) = %f %f %f \n", mat(2,0), mat(2,1), mat(2,2)); 
  */

  // transform point coordinates //

  for (int i = 0; i < number; i++) {
    PmVector3 pos = coordinates[i];
    curr_coordinates[i] = mat*(pos-xform.center) + xform.translation + xform.center;
    }

  // set potential geometry xform //

  potential_xform = xform;
  }

////////////////////////////////////////////////////////////////
//           p o i n t s           g e o m e t r y           //
//////////////////////////////////////////////////////////////

void
PmPotentialPoints::print()
  {
  fprintf (stderr, ">>> geometry = points \n");
  fprintf (stderr, ">>> center   = %f %f %f \n", center[0], center[1], center[2]);
  }

//*============================================================*
//*==========              setSphereRadius           ==========*
//*============================================================*

void
PmPotentialPoints::setSphereRadius(const float val)
  {
  //fprintf (stderr, ">>>>>> PmPotentialPoints::setSphereRadius [%s] \n", name.c_str()); 
  //fprintf (stderr, ">>> val=%f \n", val); 

  for (unsigned int i = 0; i < this->number; i++) {
    this->radii[i] = val;
    }
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*

void
PmPotentialPoints::display(const bool show)
  {
  #ifdef dbg_PmPotentialPoints_display
  fprintf (stderr, ">>>>>> PmPotentialPoints::display  show[%d] \n", show);
  fprintf (stderr, ">>> number=%d \n", number);
  fprintf (stderr, ">>> display_spheres=%d \n", display_spheres);
  fprintf (stderr, ">>> show=%d \n", show);
  #endif 

  if (!visible) return;

  if (!number) {
    pm_ErrorWarnReport (PM, "potential \"%s\" geometry has no points.", "*", 
                        name.c_str());
    return;
    }

  PmGraphicsAttributes atts;
  string geom_name = "potential";
  PmVector3 pos = center;

  geom_name = geom_name + '[' + name + ']' + "points";
  int num_verts = number;
  PmVector3 *verts = new PmVector3[num_verts];
  float *srads = new float[number];

  for (int i = 0; i < num_verts; i++) {
    verts[i] = this->coordinates[i];
    srads[i] = this->radii[i];
    }

  atts.setColor(color);

  if (display_spheres) {
    //fprintf (stderr, ">>> spheres \n");
    PmGraphicsSphere *geom = new PmGraphicsSphere(geom_name, num_verts, verts, srads);
    atts.setDisplayType(display_type);
    atts.setColor(color);
    atts.setVisible(show);
    geom->setAttributes(atts);
    geom->display();
    graphics_geometry = geom;
    }
  else {
    PmGraphicsAttributes atts;
    atts.setMarker(false);
    atts.setScale(1.0);
    atts.setColor(color);
    atts.setVisible(show);
    PmGraphicsPoint *geom = new PmGraphicsPoint(geom_name, num_verts, verts);
    geom->setAttributes(atts);
    geom->display();
    graphics_geometry = geom;
    }
  }

//*============================================================*
//*==========              getCurrentAxes            ==========*
//*============================================================*
// get the current xformed axes. 

void 
PmPotentialPoints::getCurrentAxes(PmVector3& pos, PmVector3& u, PmVector3& v, 
                                  PmVector3& w) 
  {
  PmMatrix3x3 mat = potential_xform.matrix;
  pos = mat*(center-potential_xform.center) + potential_xform.translation + 
        potential_xform.center;
  u = mat*pca.axis1;
  v = mat*pca.axis2;
  w = mat*pca.axis3;
  }

//*============================================================*
//*==========              getPcaScales              ==========*
//*============================================================*
// get the pca scales.              

void 
PmPotentialPoints:: getPcaScales(float s[3])
  {
  s[0] = rads[0];
  s[1] = rads[1];
  s[2] = rads[2];
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the xform for a points potential geometry.

void
PmPotentialPoints::setXform(PmXform& xform)
  {
  #ifdef dbg_PmPotentialPoints_setXform
  fprintf (stderr, ">>>>>> PmPotentialPoints::setXform  [%s]\n", name.c_str());
  fprintf (stderr, ">>> center=%f %f %f \n", center[0],  center[1], center[2]);  
  #endif

  PmMatrix3x3 mat = xform.matrix;
  /*
  fprintf (stderr, "   >>> mat(0,*) = %f %f %f \n", mat(0,0), mat(0,1), mat(0,2)); 
  fprintf (stderr, "   >>> mat(1,*) = %f %f %f \n", mat(1,0), mat(1,1), mat(1,2)); 
  fprintf (stderr, "   >>> mat(2,*) = %f %f %f \n", mat(2,0), mat(2,1), mat(2,2)); 
  */

  // transform center //

  curr_center = mat*(center-xform.center) + xform.translation + xform.center;

  // transform point coordinates //

  for (int i = 0; i < number; i++) {
    PmVector3 pos = coordinates[i];
    curr_coordinates[i] = mat*(pos-xform.center) + xform.translation + xform.center;
    }

  // set potential geometry xform //
  
  potential_xform = xform;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the display of a points potential geometry.

void
PmPotentialPoints::updateDisplay()
  {
  if (!graphics_geometry) {
    return;
    }

  /*
  PmMatrix3x3 mat = potential_xform.matrix;
  fprintf (stderr, ">>>>>> PmPotentialPoints::updateDisplay [%s]\n", name.c_str());
  fprintf (stderr, "   >>> mat(0,*) = %f %f %f \n", mat(0,0), mat(0,1), mat(0,2)); 
  fprintf (stderr, "   >>> mat(1,*) = %f %f %f \n", mat(1,0), mat(1,1), mat(1,2)); 
  fprintf (stderr, "   >>> mat(2,*) = %f %f %f \n", mat(2,0), mat(2,1), mat(2,2)); 
  */
  PmGraphicsSphere *geom = dynamic_cast<PmGraphicsSphere*>(graphics_geometry);

  if (!geom) {
    PmGraphicsPoint *geom = dynamic_cast<PmGraphicsPoint*>(graphics_geometry);
    geom->setXform (potential_xform);
    geom->display();
    }
  else {
    geom->setXform (potential_xform);
    geom->display();
    }
  }

////////////////////////////////////////////////////////////////
//           c o n t a c t   p o t e n t i a l               //
//////////////////////////////////////////////////////////////

typedef void (PmContactPotential::*CompForcesFptr)(PmContactPotential *pot, 
              const bool cmpf);

static CompForcesFptr 
contact_force_funcs[PM_POTENTIAL_GEOM_SIZE][PM_POTENTIAL_GEOM_SIZE] = { 

  { &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces },

    // cylinder //

  { &PmContactPotential::no_forces, 
    &PmContactPotential::cylinder_cylinder_forces, 
    &PmContactPotential::cylinder_ellipsoid_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::cylinder_plane_forces, 
    &PmContactPotential::cylinder_points_forces, 
    &PmContactPotential::cylinder_sphere_forces },

    // ellipsoid //

  { &PmContactPotential::no_forces, 
    &PmContactPotential::ellipsoid_cylinder_forces, 
    &PmContactPotential::ellipsoid_ellipsoid_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::ellipsoid_plane_forces, 
    &PmContactPotential::ellipsoid_points_forces, 
    &PmContactPotential::ellipsoid_sphere_forces },

    // lines //

  { &PmContactPotential::no_forces,
    &PmContactPotential::no_forces,
    &PmContactPotential::no_forces,
    &PmContactPotential::no_forces,
    &PmContactPotential::no_forces,
    &PmContactPotential::no_forces,
    &PmContactPotential::no_forces },

    // plane //

  { &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::plane_ellipsoid_forces,
    &PmContactPotential::no_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::plane_points_forces,
    &PmContactPotential::plane_sphere_forces},

    // points //

  { &PmContactPotential::no_forces, 
    &PmContactPotential::points_cylinder_forces, 
    &PmContactPotential::points_ellipsoid_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::points_plane_forces, 
    &PmContactPotential::points_points_forces, 
    &PmContactPotential::points_sphere_forces },

    // sphere //

  { &PmContactPotential::no_forces, 
    &PmContactPotential::sphere_cylinder_forces, 
    &PmContactPotential::sphere_ellipsoid_forces, 
    &PmContactPotential::no_forces, 
    &PmContactPotential::sphere_plane_forces, 
    &PmContactPotential::sphere_points_forces, 
    &PmContactPotential::sphere_sphere_forces } };

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmContactPotential::PmContactPotential(const string name)
  {
  this->name = name;
  strength = 1.0;
  current_energy = 0.0; 
  contact = false;
  use_impulse = false;
  restitution = 1.0;
  maximum_force = 100.0;
  spatial_decomp = NULL;
  }

//*============================================================*
//*==========              initPotential             ==========*
//*============================================================*
// initialize potential.

void
PmContactPotential::initPotential() {
  initialized = true;
  }

//*============================================================*
//*==========              checkContact              ==========*
//*============================================================*
// check for contact.                                         

void
PmContactPotential::checkContact(bool& flag) {
  flag = contact;
  }

//*============================================================*
//*==========              setImpulse                ==========*
//*============================================================*
// set impulse flag for computing contact forces using impulses.

void 
PmContactPotential::setImpulse(const bool val)
  {
  use_impulse = val;
  }

//*============================================================*
//*==========              setContactStrength        ==========*
//*============================================================*

void
PmContactPotential::setContactStrength(float val) {
  this->strength = val;
  }

//*============================================================*
//*==========          setRestitutionCoef            ==========*
//*============================================================*
// set coefficient of restitution used in computing impulse
// forces.

void 
PmContactPotential::setRestitutionCoef(const float val)
  {
  if ((val >= 0.0) && (val <= 1.0)) {
    restitution = val;
    }
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces between a set of contact potentials.

void
PmContactPotential::compForces(const bool upd, const float time, const bool active,
                               const bool verbose)
  {
  #ifdef dbg_PmContactPotential_compForces
  fprintf (stderr, ">>>>>> PmContactPotential::compForces [%s] \n", name.c_str());
  #endif
  this->verbose = verbose;

  if (!active) {
    return;
    }

  #ifdef dbg_PmContactPotential_compForces
  fprintf (stderr, ">>> contact type %s[%d]-[%d] \n", name.c_str(), geometry_type1,
           geometry_type2);
  #endif

  (this->*contact_force_funcs[geometry_type1][geometry_type2])(NULL, true);
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the geometry for a contact potential.

void
PmContactPotential::updateDisplay()
  {
  //fprintf (stderr, ">>>>>> PmContactPotential::updateDisplay[%s] \n", name.c_str());
  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              compImpulse               ==========*
//*============================================================*

void 
PmContactPotential::compImpulse (PmContactPotential *potB, PmVector3& cpt,
                                 PmVector3& normal, float *impulse)
  {

  #define ndbg_PmContactPotential_compImpulse 
  #ifdef dbg_PmContactPotential_compImpulse 
  fprintf (stderr, "\n>>>>>> compContactImpulse \n");
  #endif

  PmContactPotential *potA;
  PmMassProperties propsA, propsB;
  float massA, massB;
  PmMatrix3x3 inertiaA, inertiaB, inertiaAInv, inertiaBInv, Ia, Ib; 
  PmMatrix3x3 matrixA, matrixB, matrixATransp, matrixBTransp;
  PmXform xformA, xformB;
  PmBody *bodyA, *bodyB;
  bool singular;
  PmVector3 posA, posB, velAP, velBP, velAB, rAP, rBP;
  PmVector3 v1, v2, v3, va, vb; 
  float c1, c2, c3;
  string nameA, nameB;

  potA = this;
  potA->getName(nameA);
  potB->getName(nameB);

  #ifdef dbg_PmContactPotential_compImpulse 
  fprintf (stderr, ">>> potA  \"%s\"  \n", nameA.c_str());
  fprintf (stderr, ">>> potB  \"%s\"  \n", nameB.c_str());
  #endif

  // get information for body A //

  bodyA = (PmBody*)(potA->sobj1);
  bodyA->getMassProps(propsA);
  bodyA->getVelocity(cpt, velAP);
  bodyA->getPosition(posA);
  massA = propsA.mass;
  inertiaA = propsA.inertia;
  potA->getXform(xformA);
  matrixA = xformA.matrix;
  matrixATransp = xformA.matrix;
  matrixATransp.transpose(); 
  pm_MathMatrix3x3Inverse (inertiaA, inertiaAInv, &singular);
  Ia = matrixA * inertiaAInv * matrixATransp;
  rAP = cpt - posA;

  // get information for body B //

  bodyB = (PmBody*)(potB->sobj1);
  bodyB->getMassProps(propsB);
  bodyB->getVelocity(cpt, velBP);
  bodyB->getPosition(posB);
  massB = propsB.mass;
  inertiaB = propsB.inertia;
  potB->getXform(xformB);
  matrixB = xformB.matrix;
  matrixBTransp = xformB.matrix;
  matrixBTransp.transpose(); 
  pm_MathMatrix3x3Inverse (inertiaB, inertiaBInv, &singular);
  Ib = matrixB * inertiaBInv * matrixBTransp;
  rBP = cpt - posB;

  // compute impulse //

  velAB = velAP - velBP;
  v1 = Ia * rAP.cross(normal);
  va = v1.cross(rAP);
  v1 = Ib * rBP.cross(normal);
  vb = v1.cross(rBP);

  v3 = va + vb;
  c1 = normal*v3;
  c2 = normal*normal * (1.0 / massA + 1.0 / massB);
  c3 = normal*velAB;
  *impulse = -(1.0 + restitution)*c3 / (c1 + c2); 

  #ifdef dbg_PmContactPotential_compImpulse 
  fprintf (stderr, "       massA = %f \n", propsA.mass); 
  fprintf (stderr, "       velA  = %f %f %f \n", velAP[0], velAP[1], velAP[2]); 
  fprintf (stderr, "       massB = %f \n", propsB.mass); 
  fprintf (stderr, "       velB  = %f %f %f \n", velBP[0], velBP[1], velBP[2]); 
  fprintf (stderr, "       -----------------\n");
  fprintf (stderr, "       velAB = %f %f %f \n", velAB[0], velAB[1], velAB[2]); 
  fprintf (stderr, "       norm  = %f %f %f \n", normal[0], normal[1], normal[2]); 
  fprintf (stderr, "       c1    = %f \n", c1); 
  fprintf (stderr, "       c2    = %f \n", c2); 
  fprintf (stderr, "       c3    = %f \n", c3); 
  fprintf (stderr, "       -----------------\n");
  fprintf (stderr, "       impulse = %f \n", *impulse); 
  #endif
  }

//*============================================================*
//*==========      cylinder_cylinder_forces          ==========*
//*============================================================*
// Compute the forces between two cylinder contact potentials

void
PmContactPotential::cylinder_cylinder_forces(PmContactPotential *opot, const bool cmpf)
  {
  fprintf (stderr, "\n>>>>>> PmContactPotential::cylinder_cylinder_forces \n"); 
  PmPotentialCylinder *cyl1, *cyl2;
  PmVector3 dpos, pos1, pos2, u1, v1, w1, u2, v2, w2;
  PmVector3 force, R, npt1, npt2, cpt, dir, vec, pt;
  float *r1, *r2, dist, s, lambda, cfunc;
  float max_r1, max_r2;
  PmMatrix3x3 A, B;

  cyl1= dynamic_cast<PmPotentialCylinder*>(geometry1);
  cyl1->getGeometry(pos1, u1, v1, w1);
  r1 = cyl1->rads; 
  max_r1 = cyl1->max_rad;

  cyl2 = dynamic_cast<PmPotentialCylinder*>(geometry2);
  cyl2->getGeometry(pos2, u2, v2, w2);
  r2 = cyl2->rads; 
  max_r2 = cyl2->max_rad;

  // first do a simple extent check //

  cyl1->contact = false;
  cyl2->contact = false;

  R = pos2 - pos1;
  dist = R.length();

  if (dist > max_r1 + max_r2) {
    return;
    }

  // compute quadratic forms (A and B) for ellipsoids //

  pm_MathEcfQuadForm (r1, u1, v1, w1, A);
  pm_MathEcfQuadForm (r2, u2, v2, w2, B);

  // compute contact //

  pm_MathEcfContactFuncSolve (A, B, R, &lambda, &cfunc);

  if (cfunc > 1.0) {
    return;
    }

  // compute contact points //

  pm_MathEcfContactPtComp (A, B, pos1, pos2, lambda, cfunc, npt1, npt2, cpt);

  #ifdef dbg_PmPotential_cyl_show_contact
  cyl1->contact_geometry->update(1, &npt1);
  cyl2->contact_geometry->update(1, &npt2);
  //cyl1->contact_geometry->update(1, &cpt);
  //cyl2->contact_geometry->update(1, &cpt);
  fprintf (stderr, "   >>> cpt (%g %g %g)  \n", cpt[0], cpt[1], cpt[2]); 
  fprintf (stderr, "       npt1 (%g %g %g)  \n", npt1[0], npt1[1], npt1[2]); 
  fprintf (stderr, "       npt2 (%g %g %g)  \n", npt2[0], npt2[1], npt2[2]); 
  #endif

  cyl1->contact = true;
  cyl1->contact_pt = cpt;
  cyl2->contact = true;
  cyl2->contact_pt = cpt;

  // set contact forces //

  s = this->strength;
  dir = npt2 - npt1;
  dir.normalize();
  //mag = dir.length();
  force = s*dir;
  this->sobj1->addForce(npt1, force);
  //fprintf (stderr, "   >>> npt1 force (%g %g %g)  \n", force[0], force[1], force[2]); 
  //fprintf (stderr, "       strength = %g  \n", s); 
  //fprintf (stderr, "       mag = %g  \n", mag); 

  force = -s*dir;
  this->sobj2->addForce(npt2, force);
  }

//*============================================================*
//*==========      cylinder_ellipsoid_forces         ==========*
//*============================================================*
// Compute the forces between two cylinder contact potentials

void
PmContactPotential::cylinder_ellipsoid_forces(PmContactPotential *pot, const bool cmpf)
  {
  }

//*============================================================*
//*==========      cylinder_plane_forces             ==========*
//*============================================================*
// Compute the forces between cylinder and plane contact potentials

void
PmContactPotential::cylinder_plane_forces(PmContactPotential *pot, const bool cmpf)
  {
  }

//*============================================================*
//*==========      cylinder_points_forces            ==========*
//*============================================================*
// Compute the forces between two cylinder contact potentials

void
PmContactPotential::cylinder_points_forces(PmContactPotential *pot, const bool cmpf)
  {
  }

//*============================================================*
//*==========      cylinder_sphere_forces            ==========*
//*============================================================*
// Compute the forces between two cylinder contact potentials

void
PmContactPotential::cylinder_sphere_forces(PmContactPotential *pot, const bool cmpf)
  {
  }

//*============================================================*
//*==========      ellipsoid_cylinder_forces         ==========*
//*============================================================*

void
PmContactPotential::ellipsoid_cylinder_forces(PmContactPotential *pot, const bool cmpf)
  {
  }

//*============================================================*
//*==========      ellipsoid_ellipsoid_forces        ==========*
//*============================================================*
// Compute the forces between two ellipsoid contact potentials
// using a "hard ellipsoid" appoarach (Perram) which computes 
// the distance of closest approach of the two ellipsoids. 

void
PmContactPotential::ellipsoid_ellipsoid_forces(PmContactPotential *pot, const bool cmpf)
  {
  #define ndbg_PmContactPotential_ellipsoid_ellipsoid_forces
  #ifdef dbg_PmContactPotential_ellipsoid_ellipsoid_forces
  fprintf (stderr, ">>>>>> PmContactPotential::ellipsoid_ellipsoid_forces[%s]\n", 
           name.c_str());
  #endif

  PmPotentialEllipsoid *ellipsoid1, *ellipsoid2;
  PmVector3 dpos, pos1, pos2, u1, v1, w1, u2, v2, w2;
  PmVector3 force, R, npt1, npt2, cpt, dir, v;
  float *r1, *r2, dist, s, lambda, cfunc;
  float max_r1, max_r2, d;
  PmMatrix3x3 A, B;
  //PmMathContactEllipsoids ellipsoids;

  ellipsoid1 = dynamic_cast<PmPotentialEllipsoid*>(geometry1);
  ellipsoid1->getGeometry(pos1, u1, v1, w1);
  r1 = ellipsoid1->rads; 
  max_r1 = ellipsoid1->max_rad;

  ellipsoid2 = dynamic_cast<PmPotentialEllipsoid*>(geometry2);
  ellipsoid2->getGeometry(pos2, u2, v2, w2);
  r2 = ellipsoid2->rads; 
  max_r2 = ellipsoid2->max_rad;

  // first do a simple extent check //

  ellipsoid1->contact = false;
  ellipsoid2->contact = false;

  //pm_MathEcfCompContact (ellipsoids, contact);

  R = pos2 - pos1;
  dist = R.length();

  if (dist > max_r1 + max_r2) {
    return;
    }

  #ifdef dbg_PmContactPotential_ellipsoid_ellipsoid_forces
  fprintf (stderr, "   >>> pos1 (%g %g %g) \n", pos1[0], pos1[1], pos1[2]); 
  fprintf (stderr, "   >>> u1   (%g %g %g) \n", u1[0], u1[1], u1[2]); 
  fprintf (stderr, "   >>> v1   (%g %g %g) \n", v1[0], v1[1], v1[2]); 
  fprintf (stderr, "   >>> w1   (%g %g %g) \n", w1[0], w1[1], w1[2]); 
  fprintf (stderr, "\n");
  fprintf (stderr, "   >>> pos2 (%g %g %g) \n", pos2[0], pos2[1], pos2[2]); 
  fprintf (stderr, "   >>> u2   (%g %g %g) \n", u2[0], u2[1], u2[2]); 
  fprintf (stderr, "   >>> v2   (%g %g %g) \n", v2[0], v2[1], v2[2]); 
  fprintf (stderr, "   >>> w2   (%g %g %g) \n", w2[0], w2[1], w2[2]); 
  #endif

  // compute quadratic forms (A and B) for ellipsoids //

  pm_MathEcfQuadForm (r1, u1, v1, w1, A);
  pm_MathEcfQuadForm (r2, u2, v2, w2, B);

  // compute contact //

  pm_MathEcfContactFuncSolve (A, B, R, &lambda, &cfunc);

  if (cfunc > 1.0) {
    return;
    }

  #ifdef dbg_PmContactPotential_ellipsoid_ellipsoid_forces
  fprintf (stderr, "   >>> contact.  cfunc [%g] \n", cfunc); 
  #endif

  pm_MathEcfContactPtComp (A, B, pos1, pos2, lambda, cfunc, npt1, npt2, cpt);
  ellipsoid1->contact = true;
  ellipsoid1->contact_pt = cpt;
  ellipsoid2->contact = true;
  ellipsoid2->contact_pt = cpt;
  contact = true;

  // compute contact force using impulses //

  if (use_impulse) {
    PmVector3 normal;
    float impulse;
    normal = 2.0*A*(pos1 - cpt);
    normal.normalize();

    // compute impulse //

    this->compImpulse (pot, cpt, normal, &impulse);

    s = impulse*this->strength;

    if (s > maximum_force) {
      fprintf (stderr, "**** WARNING: potential \"%s\" force is large.  s1 = %g \n", 
               name.c_str(), s);
      fprintf (stderr, "     impulse = %g \n", impulse);
      }

    dir = normal;
    current_energy = 0.0; 

    #ifdef dbg_PmContactPotential_compImpulse 
    static PmGraphicsLine *lineAP = NULL;
    static PmGraphicsLine *lineNormal = NULL;
    string geom_name;
    PmVector3 verts[2];
    PmGraphicsAttributes atts;
    PmVector3 color;
  
    color.set(1,1,1);
    atts.setColor(color);

    verts[0] = pos1;
    verts[1] = cpt;

    if (!lineAP) {
      geom_name = "contact_rAP";
      lineAP = new PmGraphicsLine(geom_name, 2, verts);
      lineAP->setAttributes(atts);
      lineAP->display();
      }
    else {
      lineAP->update(2, verts);
      }

    verts[0] = cpt;
    verts[1] = cpt + normal;

    if (!lineNormal) {
      geom_name = "contact_Normal";
      lineNormal = new PmGraphicsLine(geom_name, 2, verts);
      lineNormal->setAttributes(atts);
      lineNormal->display();
      }
    else {
      lineNormal->update(2, verts);
      }
    #endif
    }

  // compute contact force as a constant //

  else {
    dir = npt2 - npt1;
    d = dir.length();
    dir.normalize();
    s = d*this->strength;
    current_energy = this->strength;
    }

  // set contact forces //

  if (cmpf) {
    force = s*dir;
    this->sobj1->addForce(npt1, force);
    force = -s*dir;
    this->sobj2->addForce(npt2, force);
    }

  if (verbose) {
    string pname;
    fprintf (stderr, ">>> potential \"%s\" in contact. \n", name.c_str());
    }
  }

//*============================================================*
//*==========      ellipsoid_plane_forces            ==========*
//*============================================================*

void
PmContactPotential::ellipsoid_plane_forces(PmContactPotential *pot, const bool cmpf)
  {
  plane_ellipsoid_forces(this, cmpf);
  }

//*============================================================*
//*==========      ellipsoid_sphere_forces           ==========*
//*============================================================*
// Compute the forces between ellipsoid and sphere contact 
// potentials using a "hard ellipsoid" appoarach (Perram).

void
PmContactPotential::ellipsoid_sphere_forces(PmContactPotential *opot, const bool cmpf)
  {
  #define ndbg_PmContactPotential_ellipsoid_sphere_forces
  #ifdef dbg_PmContactPotential_ellipsoid_sphere_forces
  fprintf (stderr, ">>>>>> PmContactPotential::ellipsoid_sphere_forces[%s]\n", 
           name.c_str());
  #endif

  PmPotentialEllipsoid *ellipsoid;
  PmPotentialSphere *sphere;
  PmVector3 dpos, pos1, pos2, u1, v1, w1, u2, v2, w2;
  PmVector3 force, R, npt1, npt2, cpt, dir;
  float *r1, r2[3], dist, s, lambda, cfunc;
  float max_r1, max_r2;
  PmMatrix3x3 A, B;
  PmSimulationObj *sobj1, *sobj2;

  ellipsoid = dynamic_cast<PmPotentialEllipsoid*>(geometry1);

  if (ellipsoid == NULL) {
    ellipsoid = dynamic_cast<PmPotentialEllipsoid*>(geometry2);
    sphere = dynamic_cast<PmPotentialSphere*>(geometry1);
    sobj2 = this->sobj1;
    sobj1 = this->sobj2;
    }
  else {
    sphere = dynamic_cast<PmPotentialSphere*>(geometry2);
    sobj1 = this->sobj1;
    sobj2 = this->sobj2;
    }
 
  ellipsoid->getGeometry(pos1, u1, v1, w1);
  r1 = ellipsoid->rads; 
  max_r1 = ellipsoid->max_rad;

  r2[0] = sphere->radius;
  r2[1] = sphere->radius;
  r2[2] = sphere->radius;
  max_r2 = sphere->radius;
  pos2 = sphere->curr_center;
  u2.set(1.0, 0.0, 0.0);
  v2.set(0.0, 1.0, 0.0);
  w2.set(0.0, 0.0, 1.0);

  // first do a simple extent check //

  ellipsoid->contact = false;
  sphere->contact = false;

  R = pos2 - pos1;
  dist = R.length();

  if (dist > max_r1 + max_r2) {
    return;
    }

  // build quadratic forms (A and B) for ellipsoids //

  pm_MathEcfQuadForm (r1, u1, v1, w1, A);
  pm_MathEcfQuadForm (r2, u2, v2, w2, B);

  // compute contact //

  pm_MathEcfContactFuncSolve (A, B, R, &lambda, &cfunc);
  //fprintf (stderr, "   >>> cfunc  [%g] \n", cfunc); 
  //fprintf (stderr, "   >>> lambda [%g] \n", lambda); 

  if (cfunc > 1.0) {
    return;
    }

  #ifdef dbg_PmContactPotential_ellipsoid_ellipsoid_forces
  fprintf (stderr, "   >>> contact.  cfunc [%g] \n", cfunc); 
  #endif

  pm_MathEcfContactPtComp (A, B, pos1, pos2, lambda, cfunc, npt1, npt2, cpt);
  ellipsoid->contact = true;
  ellipsoid->contact_pt = cpt;
  sphere->contact = true;
  //sphere->contact_pt = cpt;

  // set contact forces //

  s = this->strength;
  dir = npt2 - npt1;
  dir.normalize();
  force = s*dir;
  sobj1->addForce(npt1, force);

  force = -s*dir;
  sobj2->addForce(npt2, force);
  #ifdef dbg_PmContactPotential_ellipsoid_sphere_forces
  fprintf (stderr, "   >>> fmag = %g \n", s); 
  #endif
  }

//*============================================================*
//*==========      ellipsoid_points_forces           ==========*
//*============================================================*
// Compute the forces between ellipsoid and points contact 
// potentials using a "hard ellipsoid" appoarach (Perram).

void
PmContactPotential::ellipsoid_points_forces(PmContactPotential *opot, const bool cmpf)
  {
  #define ndbg_PmContactPotential_ellipsoid_points_forces
  #ifdef dbg_PmContactPotential_ellipsoid_points_forces
  fprintf (stderr, "\n>>>>>> PmContactPotential::ellipsoid_points_forces[%s]\n", 
           name.c_str());
  #endif

  PmPotentialEllipsoid *ellipsoid;
  PmPotentialPoints *points;
  PmVector3 dpos, pos1, pos2, u1, v1, w1, u2, v2, w2;
  PmVector3 force, R, npt1, npt2, cpt, dir, vec;
  float *r1, r2[3], dist, s, lambda, cfunc;
  float max_r1, max_r2;
  PmMatrix3x3 A, B;
  PmSimulationObj *sobj1, *sobj2;

  ellipsoid = dynamic_cast<PmPotentialEllipsoid*>(geometry1);

  if (ellipsoid == NULL) {
    ellipsoid = dynamic_cast<PmPotentialEllipsoid*>(geometry2);
    points = dynamic_cast<PmPotentialPoints*>(geometry1);
    sobj2 = this->sobj1;
    sobj1 = this->sobj2;
    }
  else {
    points = dynamic_cast<PmPotentialPoints*>(geometry2);
    sobj1 = this->sobj1;
    sobj2 = this->sobj2;
    }

  ellipsoid->getGeometry(pos1, u1, v1, w1);
  r1 = ellipsoid->rads; 
  max_r1 = ellipsoid->max_rad;

  u2.set(1.0, 0.0, 0.0);
  v2.set(0.0, 1.0, 0.0);
  w2.set(0.0, 0.0, 1.0);
  max_r2 = points->radius;
  pos2 = points->curr_center;

  // first do a simple extent check //

  ellipsoid->contact = false;
  points->contact = false;

  R = pos2 - pos1;
  dist = R.length();

  if (dist > max_r1 + max_r2) {
    return;
    }

  // build quadratic form for ellipsoid //

  pm_MathEcfQuadForm (r1, u1, v1, w1, A);

  // check for contact with each point and the ellipsoid  //
 
  for (int k = 0; k < points->number; k++) {
    pos2 = points->curr_coordinates[k];
    r2[0] = points->radii[k];
    r2[1] = points->radii[k];
    r2[2] = points->radii[k];
    R = pos2 - pos1;
    dist = R.length();

    // quick extent check //
    if (dist > max_r1 + r2[0]) {
      continue;
      }

    // bounding box check //
    if ((fabs(R*u1) >= r1[0]+r2[0]) ||
        (fabs(R*v1) >= r1[1]+r2[0]) ||
        (fabs(R*w1) >= r1[2]+r2[0])) {
      continue;
      }

    // build quadratic form for sphere //

    pm_MathEcfQuadForm (r2, u2, v2, w2, B);

    // compute contact //

    pm_MathEcfContactFuncSolve (A, B, R, &lambda, &cfunc);

    if (cfunc <= 1.0) {
      #ifdef dbg_PmContactPotential_ellipsoid_ellipsoid_forces
      fprintf (stderr, "   >>> contact[%d]  cfunc [%g] \n", k, cfunc); 
      #endif

      pm_MathEcfContactPtComp (A, B, pos1, pos2, lambda, cfunc, npt1, npt2, cpt);
      ellipsoid->contact = true;
      ellipsoid->contact_pt = cpt;
      points->contact = true;

      // set contact forces //

      s = this->strength;
      dir = npt2 - npt1;
      dir.normalize();
      force = s*dir;
      sobj1->addForce(npt1, force);

      force = -s*dir;
      sobj2->addForce(pos2, force);
      }
    }
  }

//*============================================================*
//*==========            plane_points_forces         ==========*
//*============================================================*

void
PmContactPotential::plane_points_forces(PmContactPotential *pot, const bool cmpf)
  {
  //fprintf (stderr, "\n>>>>>> PmPotential::plane_points_forces[%s] \n", name.c_str());
  float d, s, r, max_d;
  PmPotentialPlane *plane;
  PmPotentialPoints *points;
  PmVector3 normal, pos;
  PmVector3 force; 
  PmSimulationObj *sobj1, *sobj2;
  int num_contact;

  // check which geometry is the plane //

  points = dynamic_cast<PmPotentialPoints*>(geometry1);

  if (!points) {
    points = dynamic_cast<PmPotentialPoints*>(geometry2);
    plane = dynamic_cast<PmPotentialPlane*>(geometry1);
    sobj2 = this->sobj1;
    sobj1 = this->sobj2;
    }
  else {
    plane = dynamic_cast<PmPotentialPlane*>(geometry2);
    sobj1 = this->sobj1;
    sobj2 = this->sobj2;
    }

  points->contact = false;
  s = this->strength;
  normal = plane->normal;
  num_contact = 0;
  max_d = 0.0;

  for (int i = 0; i < points->number; i++) {
    r = points->radii[i];
    pos = points->curr_coordinates[i];
    d = normal*(pos - plane->center) - r;

    if (d < 0.0) {
      if (fabs(d) > 1.0) {
        force = fabs(d)*s*normal;
        }
      else {
        force = s*normal;
        }

      if (fabs(d) > max_d) max_d = fabs(d);
      this->sobj2->addForce(pos, force);
      points->contact = true;
      num_contact += 1;
      }
    }

  /*
  if (num_contact) {
    fprintf (stderr, ">>>>>> PmPotential::plane_points_forces num contact=%d maxd=%f\n", 
             num_contact, max_d);
    }
  */
  }

//*============================================================*
//*==========      plane_sphere_forces               ==========*
//*============================================================*

void
PmContactPotential::plane_ellipsoid_forces(PmContactPotential *pot, const bool cmpf)
  {
  //fprintf (stderr, "\n>>>>>> PmPotential::plane_ellipsoid_forces[%s] \n", name.c_str());
  float d, r, max_rad, s, *rads;
  PmPotentialPlane *plane;
  PmPotentialEllipsoid *ellipsoid;
  PmVector3 normal, pos, u, v, w, spts[40], pt; 
  PmVector3 force, fpt; 
  int n;

  ellipsoid = dynamic_cast<PmPotentialEllipsoid*>(geometry2);
  plane = dynamic_cast<PmPotentialPlane*>(geometry1);
  rads = ellipsoid->rads;
  max_rad = ellipsoid->max_rad;
  ellipsoid->getGeometry(pos, u, v, w);

  normal = plane->normal;
  ellipsoid->contact = false;
  d = normal*(pos - plane->center) - max_rad;

  if (d >= 0.0) {
    return;
    }

  n = 0;
  spts[n++] = pos + rads[0]*u;
  spts[n++] = pos - rads[0]*u;
  spts[n++] = pos + rads[1]*v;
  spts[n++] = pos - rads[1]*v;
  spts[n++] = pos + rads[2]*w;
  spts[n++] = pos - rads[2]*w;

  for (int i = 0; i < n; i++) {
    pt = spts[i];
    d = normal*(pt - plane->center);

    if (d < 0.0) {
      ellipsoid->contact = true;
      fpt = pt;
      break;
      }
    }

  if (!ellipsoid->contact) {
    return; 
    } 

  s = this->strength;
  force = s*normal;
  this->sobj2->addForce(fpt, force);
  }

//*============================================================*
//*==========      plane_sphere_forces               ==========*
//*============================================================*

void
PmContactPotential::plane_sphere_forces(PmContactPotential *pot, const bool cmpf)
  {
  //fprintf (stderr, "\n>>>>>> PmPotential::plane_sphere_forces[%s] \n", name.c_str());

  float d, rad, s;
  PmPotentialPlane *plane;
  PmPotentialSphere *sphere;
  PmVector3 normal, fpt, pos, force;

  sphere = dynamic_cast<PmPotentialSphere*>(geometry2);
  plane = dynamic_cast<PmPotentialPlane*>(geometry1);
  rad = sphere->radius;
  pos = sphere->curr_center;
  normal = plane->normal;
  d = normal*(pos - plane->center) - rad;

  if (d >= 0.0) {
    sphere->contact = false;
    return;
    }

  s = this->strength;
  force = s*normal;
  fpt = pos - normal*rad;
  this->sobj2->addForce(fpt, force);
  sphere->contact = true;

  /*
  fprintf (stderr, "  >>> contact  d = %g \n", d);
  fprintf (stderr, "  >>> force = %g %g %g \n", force[0], force[1], force[2] );
  */
  }

//*============================================================*
//*==========      points_cylinder_forces            ==========*
//*============================================================*

void
PmContactPotential::points_cylinder_forces(PmContactPotential *pot, const bool cmpf)
  {
  cylinder_points_forces(this, cmpf);
  }

//*============================================================*
//*==========      points_ellipsoid_forces           ==========*
//*============================================================*

void
PmContactPotential::points_ellipsoid_forces(PmContactPotential *pot, const bool cmpf) {
  ellipsoid_points_forces(this, cmpf);
  }

//*============================================================*
//*==========      points_plane_forces               ==========*
//*============================================================*

void
PmContactPotential::points_plane_forces(PmContactPotential *pot, const bool cmpf) 
  {
  plane_points_forces(this, cmpf);
  }

//*============================================================*
//*==========      points_points_forces              ==========*
//*============================================================*
// compute contact forces between two points potentials. 

void
PmContactPotential::points_points_forces(PmContactPotential *opot, const bool cmpf)
  {
  #ifdef dbg_PmContactPotential_points_points_forces
  fprintf (stderr, "\n>>>>>> PmPotential::points_points_forces[%s]\n", name.c_str());
  #endif

  PmPotentialPoints *points1, *points2;
  float r1, r2, dist, s;
  PmVector3 pos1, pos2, dpos, force;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points1->contact = false;
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  points2->contact = false;
  s = this->strength;

  // do a an extent overlap test                   //
  // NOTE: this is only valid for rigid xforms.    //

  dpos = points1->curr_center - points2->curr_center;
  dist = dpos.length();

  if (dist > points1->radius + points2->radius) {
    return;
    }

  // check for overlap of ellipsoid extents //

  PmMathContactEllipsoids ellipsoids;
  bool debug, contact;
  debug = true;

  if (points1->number > 10) {
    points1->getCurrentAxes(ellipsoids.pos1, ellipsoids.u1, ellipsoids.v1, ellipsoids.w1);
    points1->getPcaScales(ellipsoids.r1);
    points2->getCurrentAxes(ellipsoids.pos2, ellipsoids.u2, ellipsoids.v2, ellipsoids.w2);
    points2->getPcaScales(ellipsoids.r2);

    for (int i = 0; i < 3; i++) {
      ellipsoids.r1[i] += 0.4;
      ellipsoids.r2[i] += 0.4;
      }

    pm_MathEcfCompContact (ellipsoids, contact);

    if (!contact) {
      //fprintf (stderr, "     no ellipsoid overlap \n");
      return;
      }
    }

#ifdef brute_PmContactPotential_points_points_forces

  // check all points //

  for (int i = 0; i < points1->number; i++) {
    r1 = points1->radii[i];
    pos1 = points1->curr_coordinates[i];

    for (int j = 0; j < points2->number; j++) {
      r2 = points2->radii[j];
      pos2 = points2->curr_coordinates[j];
      dpos = pos2 - pos1; 
      dist = dpos.length();

      if (dist < r1 + r2) {
        //fprintf (stderr, "   >>> contact  %d %d \n", i, j);
        points1->contact = true;
        points2->contact = true;
        force = s*dpos / (dist*dist);
        this->sobj2->addForce(pos2, force);
        force = -force;
        this->sobj1->addForce(pos1, force);
        }
      }
    }
#endif

  // use spatial decomposition //

  if (!spatial_decomp) {
    this->initSpatialDecomp(points1, &spatial_decomp, 0.0);
    }

  // note: is this the right geometry? //
  spatial_decomp->setXform (this->geometry1->potential_xform);
  vector<int> ids;
  int num_sd_contact = 0;
  int id;
  float min_dist = 1e6;
  //fprintf (stderr, "\n---- spatial decomp ---------- \n");

  for (int i = 0; i < points2->number; i++) {
    r2 = points2->radii[i];
    pos2 = points2->curr_coordinates[i];
    spatial_decomp->getPoints(pos2, r2, ids);
    //fprintf (stderr, "%d: num ids = %d \n", i, ids.size());

    if (ids.size() != 0) {
      for (unsigned int j = 0; j < ids.size(); j++) {
        id = ids[j];
        r1 = points1->radii[id];
        pos1 = points1->curr_coordinates[id];
        dpos = pos2 - pos1;
        dist = dpos.length();

        if (dist < r1 + r2) {
          points1->contact = true;
          points2->contact = true;
          #ifdef dbg_PmContactPotential_points_points_forces
          fprintf (stderr, ">>> contact \n");
          #endif

          if (dist  < min_dist) { 
            min_dist = dist; 
            }

          dpos.normalize();
          force = -s*dpos / (dist*dist);
          this->sobj1->addForce(pos1, force);

          #ifdef dbg_PmContactPotential_points_points_forces
          fprintf (stderr, ">>> force on sobj1 = %f %f %f\n", force[0],force[1],force[2]); 
          #endif

          force = -force;
          this->sobj2->addForce(pos2, force);
          num_sd_contact += 1;
          }
        }
      }
    }

  if (num_sd_contact) {
    points1->contact = true;
    points2->contact = true;
    //fprintf (stderr, ">>> num_sd_contact = %d \n", num_sd_contact); 
    }
  }

//*============================================================*
//*==========      points_sphere_forces              ==========*
//*============================================================*
// compute contact forces between a points and a sphere 
// potential. 

void
PmContactPotential::points_sphere_forces(PmContactPotential *opot, const bool cmpf)
  {
  //fprintf (stderr, ">>>>>> PmPotential::points_sphere_forces[%s]\n", name.c_str());
  PmPotentialSphere *sphere;
  PmPotentialPoints *points;
  float sphere_r, r, dist, s;
  PmVector3 sphere_pos, pos, dpos, force;
  PmSimulationObj *sobj1, *sobj2;

  points = dynamic_cast<PmPotentialPoints*>(geometry1);

  if (!points) {
    points = dynamic_cast<PmPotentialPoints*>(geometry2);
    sphere = dynamic_cast<PmPotentialSphere*>(geometry1);
    sobj2 = this->sobj1;
    sobj1 = this->sobj2;
    }
  else {
    sphere = dynamic_cast<PmPotentialSphere*>(geometry2);
    sobj1 = this->sobj1;
    sobj2 = this->sobj2;
    }

  points->contact = false;

  sphere_r = sphere->radius;
  sphere_pos = sphere->curr_center;
  sphere->contact = false;
  s = this->strength;

  for (int i = 0; i < points->number; i++) {
    r = points->radii[i];
    pos = points->curr_coordinates[i];
    dpos = sphere_pos - pos; 
    dist = dpos.length();

    if (dist < r + sphere_r) {
      //fprintf (stderr, "   >>> contact. \n");
      sphere->contact = true;
      points->contact = true;
      force = s*dpos;
      sobj2->addForce(sphere_pos, force);
      force = -force;
      sobj1->addForce(pos, force);
      }
    }
  }

//*============================================================*
//*==========      sphere_cylinder_forces            ==========*
//*============================================================*

void
PmContactPotential::sphere_cylinder_forces(PmContactPotential *pot, const bool cmpf)
  {
  cylinder_sphere_forces(this, cmpf);
  }

//*============================================================*
//*==========      sphere_ellipsoid_forces           ==========*
//*============================================================*

void
PmContactPotential::sphere_ellipsoid_forces(PmContactPotential *pot, const bool cmpf)
  {
  ellipsoid_sphere_forces(this, cmpf);
  }

//*============================================================*
//*==========      sphere_plane_forces               ==========*
//*============================================================*

void
PmContactPotential::sphere_plane_forces(PmContactPotential *pot, const bool cmpf)
  {
  plane_sphere_forces(this, cmpf);
  }

//*============================================================*
//*==========      sphere_points_forces              ==========*
//*============================================================*

void
PmContactPotential::sphere_points_forces(PmContactPotential *pot, const bool cmpf)
  {
  points_sphere_forces(this, cmpf);
  }

//*============================================================*
//*==========         sphere_sphere_forces           ==========*
//*============================================================*
// compute the forces between two sphere contact potentials.

void 
PmContactPotential::sphere_sphere_forces(PmContactPotential *pot, const bool cmpf)
  {
  #ifdef dbg_PmContactPotential_sphere_sphere_forces
  fprintf (stderr, ">>>>>> PmPotential::sphere_sphere_forces[%s]\n", name.c_str());
  fprintf (stderr, ">>> strength=%f \n", strength); 
  #endif
  PmPotentialSphere *sphere1, *sphere2;
  PmVector3 dir, pos1, pos2;
  PmVector3 force, cpt; 
  float r1, r2, dist, s;

  sphere1 = dynamic_cast<PmPotentialSphere*>(geometry1);
  r1 = sphere1->radius;
  pos1 = sphere1->curr_center;

  sphere2 = dynamic_cast<PmPotentialSphere*>(geometry2);
  r2 = sphere2->radius;
  pos2 = sphere2->curr_center;

  #ifdef dbg_PmContactPotential_sphere_sphere_forces
  fprintf (stderr, ">>> r1=%f pos1 (%g %g %g) \n", r1, pos1[0], pos1[1], pos1[2]); 
  fprintf (stderr, ">>> r2=%f pos2 (%g %g %g) \n", r2, pos2[0], pos2[1], pos2[2]); 
  #endif

  dir = pos2 - pos1; 
  dist = dir.length();

  #ifdef dbg_PmContactPotential_sphere_sphere_forces
  fprintf (stderr, ">>> dist=%f \n", dist); 
  #endif

  if (dist < r1 + r2) {
    sphere1->contact = true;
    sphere2->contact = true;

    //===== compute contact force using impulses =====//

    if (use_impulse) {
      PmVector3 normal;
      float impulse;
      normal = dir; 
      normal.normalize();
      cpt = pos1 + r1*normal;
  
      // compute impulse //

      this->compImpulse (pot, cpt, normal, &impulse);

      s = -impulse*this->strength;
      current_energy = 0.0; 

      static PmGraphicsLine *lineNormal = NULL;
      string geom_name;
      PmVector3 verts[2];
      PmGraphicsAttributes atts;
      PmVector3 color;
  
      color.set(1,1,1);
      atts.setColor(color);

      verts[0] = pos1;
      verts[1] = cpt;

      /*
      if (!lineAP) {
        geom_name = "contact_rAP";
        lineAP = new PmGraphicsLine(geom_name, 2, verts);
        lineAP->setAttributes(atts);
        lineAP->display();
        }
      else {
        lineAP->update(2, verts);
        }
      */

      verts[0] = cpt;
      verts[1] = cpt + normal;
      color.set(1,1,0);
      atts.setColor(color);

      if (!lineNormal) {
        geom_name = "contact_Normal";
        lineNormal = new PmGraphicsLine(geom_name, 2, verts);
        lineNormal->setAttributes(atts);
        lineNormal->display();
        }
      else {
        lineNormal->update(2, verts);
        }
      }

    //===== contact using strength =====//
    else {
      s = this->strength / (dist*dist);
      }

    contact = true;

    if (verbose) {
      fprintf (stderr, ">>> potential \"%s\" in contact. \n", name.c_str());
      }

    if (cmpf) { 
      force = -s*dir;
      this->sobj1->addForce(pos1, force);
      #ifdef dbg_PmContactPotential_sphere_sphere_forces
      fprintf (stderr, ">>> contact. \n");
      fprintf (stderr, ">>> force=%f \n", s);
      #endif
      force = s*dir;
      this->sobj2->addForce(pos2, force);
      }

    current_energy = this->strength / dist;
    }
  else {
    sphere1->contact = false;
    sphere2->contact = false;
    current_energy = 0.0; 
    contact = false;
    }
  }

////////////////////////////////////////////////////////////////
//   e n t r o p i c     s p r i n g     p o t e n t i a l   //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmEntropicSpringPotential::PmEntropicSpringPotential(const string name)
  {
  //fprintf (stderr, ">>> PmEntropicSpringPotential: ctor [%s] \n", name.c_str());
  this->name = name;
  this->force_const = 1.0;
  this->graphics_conn_geometry = NULL;
  this->show_conn = false;
  this->conn_color.set(1,1,1);
  this->conn_line_width = 1.0;
  this->force_function_type = PM_ENTROPIC_SPRING_FUNC_UNKNOWN;
  this->force_function = NULL;
  }

//*============================================================*
//*==========              convFunctionType          ==========*
//*============================================================*

void 
PmEntropicSpringPotential::convFunctionType(const string str, 
                                            PmEntropicSpringFuncType& type) 
  {
  type = PM_ENTROPIC_SPRING_FUNC_UNKNOWN;

  if (str == "graph") {
    type = PM_ENTROPIC_SPRING_FUNC_GRAPH;
    }

  else if (str == "wlc") {
    type = PM_ENTROPIC_SPRING_FUNC_WLC;
    }
  }

//*============================================================*
//*==========          setForceFunctionType          ==========*
//*============================================================*

void 
PmEntropicSpringPotential::setForceFunctionType(PmEntropicSpringFuncType type)
  {
  #ifdef dbg_PmEntropicSpringPotential_setForceFunctionType
  fprintf (stderr, ">>>>>> PmEntropicSpringPotential::setForceFunctionType: \n"); 
  fprintf (stderr, ">>> type=%d \n", type); 
  #endif
  force_function_type = type;

  switch (type) {
    case PM_ENTROPIC_SPRING_FUNC_GRAPH:
      {
      EntropicSpringForceGraph *graph = new EntropicSpringForceGraph;
      force_function = graph;
      //fprintf (stderr, ">>> setForceFunction graph \n");
      }
    break;

    case PM_ENTROPIC_SPRING_FUNC_WLC:
      {
      EntropicSpringForceWlc *wlc = new EntropicSpringForceWlc;
      force_function = wlc;
      //fprintf (stderr, ">>> setForceFunction wlc \n");
      }
    break;

    default:
      fprintf (stderr, ">>> setForceFunction ? \n");
    break;
    }
  }

//*============================================================*
//*==========          setForceGraphData             ==========*
//*============================================================*

void
PmEntropicSpringPotential::setForceGraphData(const vector<float>& data)
  {

  if (force_function_type != PM_ENTROPIC_SPRING_FUNC_GRAPH) {
    return;
    }

  EntropicSpringForceGraph *func = (EntropicSpringForceGraph*)force_function; 
  func->graph_data.push_back(0.0);
  func->graph_data.push_back(0.0);

  for (unsigned int i = 0; i < data.size(); i++) {
    func->graph_data.push_back(data[i]);
    }

  func->num_graph_data = func->graph_data.size() / 2;
  }

//*============================================================*
//*==========          setWlcParams                  ==========*
//*============================================================*

void 
PmEntropicSpringPotential::setWlcParams(const float lp, const float lc)
  {
  EntropicSpringForceWlc *func = (EntropicSpringForceWlc*)force_function; 
  func->persistence_length = lp;
  func->contour_length = lc;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the geometry for a spring potential.

void
PmEntropicSpringPotential::updateDisplay()
  {
  #ifdef dbg_PmEntropicSpringPotential_updateDisplay
  fprintf (stderr, ">>> PmSpringPotential:: updateDisplay[%s] \n", name.c_str());
  fprintf (stderr, "    show conn = %d \n", show_conn); 
  #endif

  if (this->show_conn) {
    displayConnectivity();
    }

  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces between two spring potentials.

void 
PmEntropicSpringPotential::compForces(const bool upd, const float time, const bool active,
                                      const bool verbose)
  {
  #ifdef dbg_PmEntropicSpringPotential_compForces
  fprintf (stderr, ">>>>>> PmSpringPotential::compForces name=%s  time = %f \n", 
           name.c_str(), time);
  #endif

  if (!active) {
    return;
    }

  this->verbose = verbose;

  // compute forces //

  PmPotentialPoints *points1, *points2;
  PmVector3 pos1, pos2, dir, force;
  float dist, s, mag, fmag;
  float energy;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  s = this->force_const;

  EntropicSpringForce *func = (EntropicSpringForce*)force_function; 

  if (!func->initialized) {
    dir = points2->coordinates[0] - points1->coordinates[0];
    mag = dir.length();
    func->initialize(mag);
    }

  pos1 = points1->curr_coordinates[0];
  pos2 = points2->curr_coordinates[0];
  dir = pos2 - pos1;
  dist = dir.length();
  func->compForce(dist, fmag, energy);
  force = fmag*dir;
  this->sobj1->addForce(pos1, force);
  force = -force;
  this->sobj2->addForce(pos2, force);
  this->current_energy = energy; 

  #ifdef dbg_PmEntropicSpringPotential_compForces
  fprintf (stderr, ">>> entropic dist = %f  fmag = %f \n", dist, fmag); 
  fprintf (stderr, ">>> force=%f %f %f \n", force[0], force[1], force[2]); 
  #endif
  }

//*============================================================*
//*==========              setForceConst             ==========*
//*============================================================*
// set the force constant. 

void 
PmEntropicSpringPotential::setForceConst(const float val) {
  this->force_const = val;
  }

//*============================================================*
//*==========              setDisplayConn            ==========*
//*============================================================*
// set show the spring connectivity.

void
PmEntropicSpringPotential::setDisplayConn(const bool flag) {
  //fprintf (stderr, ">>>>>> PmSpringPotential::setDisplayConn  flag=%d \n", flag);
  show_conn = flag;
  }

//*============================================================*
//*==========              displayConn               ==========*
//*============================================================*
// show the spring connectivity.

void
PmEntropicSpringPotential::displayConnectivity()
  {
  //fprintf (stderr, ">>>>>> PmSpringPotential::displayConnectivity \n");
  PmPotentialPoints *points1, *points2;
  PmVector3 pos1, pos2;  

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);

  if (!graphics_conn_geometry) {
    connectivity_verts = new PmVector3[2];
    }

  connectivity_verts[0] = points1->curr_coordinates[0];
  connectivity_verts[1] = points2->curr_coordinates[0];

  if (!graphics_conn_geometry) { 
    string geom_name = "entropic spring potential";
    geom_name = geom_name + '[' + name + ']' + "spring_conn";

    PmGraphicsAttributes atts;
    atts.setColor(this->geometry1->color); 
    atts.setLineWidth(this->geometry1->line_width); 

    PmGraphicsLine *geom = new PmGraphicsLine(geom_name, 2, connectivity_verts);
    geom->setAttributes(atts);
    geom->display();
    graphics_conn_geometry = geom;
    }

  else {
    PmGraphicsLine *geom = dynamic_cast<PmGraphicsLine*>(graphics_conn_geometry);
    geom->update(2, connectivity_verts);
    }
  }

////////////////////////////////////////////////////////////////
//            entropic spring force functions                //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

EntropicSpringForce::EntropicSpringForce()
  {
  initialized = false;
  }

EntropicSpringForceGraph::EntropicSpringForceGraph()
  {
  }

EntropicSpringForceWlc::EntropicSpringForceWlc()
  {
  persistence_length = 0.4;
  contour_length = 0.34;
  }

//*============================================================*
//*==========             compEnergy                 ==========*
//*============================================================*

void
EntropicSpringForceGraph::compEnergy(const float dist, float& energy)
  {
  energy = 0;
  }

void 
EntropicSpringForceWlc::compEnergy(const float dist, float& energy)
  {
  //fprintf (stderr, ">>>>>> EntropicSpringForceWlc::compEnergy \n");
  float kt = 4.1;
  float s, c;
  float extension;

  extension = dist - length;

  if (extension < 0) {
    energy = 0;
    return;
    }

  s = kt / persistence_length;
  c = extension / contour_length;
  energy = s * ( 0.5*extension*c - 0.25*extension + 0.25*contour_length / (1.0-c) 
                 - 0.25*contour_length);
  }

//*============================================================*
//*==========              compForce                 ==========*
//*============================================================*
// compute the force for the given distance.

void 
EntropicSpringForceGraph::compForce(const float dist, float& force, float& energy)
  {
  float s, d1, d2, f1, f2;
  float extension;
  force = 0.0;
  extension = length - dist;

  if (extension <= graph_data[0]) {
    return;
    }

  for (int i = 0; i < num_graph_data-1; i++) {
    d1 = graph_data[2*i];
    f1 = graph_data[2*i+1];
    d2 = graph_data[2*(i+1)];
    f2 = graph_data[2*(i+1)+1];

    if (extension <= d2) {
      break;
      }
    }

  s = (f2 - f1) / (d2 - d1);
  force = f1 + (extension - d1) * s; 
  }

// worm-like chain //

void
EntropicSpringForceWlc::compForce(const float dist, float& force, float& energy)
  {
  #ifdef dbg_EntropicSpringForceWlc_compForce
  fprintf (stderr, ">>>>>> EntropicSpringForceWlc::compForce \n");
  fprintf (stderr, ">>> contour_length=%f \n", contour_length); 
  fprintf (stderr, ">>> persistence_length=%f \n", persistence_length); 
  #endif

  if ((persistence_length == 0.0) || (contour_length == 0.0)) {
    pm_ErrorWarnReport (PM, "wlc entropic spring parametres not set.", "*");
    return;
    }

  float kt = 4.1;
  float s, c, c2;
  float extension;

  extension = dist - this->length;
  #ifdef dbg_EntropicSpringForceWlc_compForce
  fprintf (stderr, ">>> extension=%f \n", extension); 
  #endif

  if (extension < 0) {
    force = 0;
    energy = 0;
    return;
    }

  s = kt / persistence_length;
  c = extension / contour_length;
  c2 = (1.0 - c);
  force = s * (0.25 / c2*c2 - 0.25 + c);
  energy = s * ( 0.5*extension*c - 0.25*extension + 0.25*contour_length / (1.0-c) 
                 - 0.25*contour_length);

  #ifdef dbg_EntropicSpringForceWlc_compForce
  fprintf (stderr, ">>> extension=%f  force=%f \n", extension, force);
  fprintf (stderr, ">>> length=%f  dist=%f \n", length, dist);
  fprintf (stderr, ">>> s=%f  c=%f \n", s, c);
  fprintf (stderr, ">>> force=%f \n", force);
  fprintf (stderr, ">>> energy=%f \n", energy);
  #endif
  }

//*============================================================*
//*==========              initialize                ==========*
//*============================================================*

void 
EntropicSpringForceGraph::initialize(const float length)
  {
  initialized = true;
  graph_data[0] = length;
  this->length = length;
  }

void
EntropicSpringForceWlc::initialize(const float length)
  {
  initialized = true;
  this->length = length;
  }

////////////////////////////////////////////////////////////////
//           s p r i n g     p o t e n t i a l               //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSpringPotential::PmSpringPotential(const string name)
  {
  //fprintf (stderr, ">>> PmSpringPotential: ctor [%s] \n", name.c_str());
  this->name = name;
  this->force_const = 1.0;
  this->cutoff = 0.8;
  this->number_springs = 0;
  this->no_connections_warning = false;
  this->connectivity_computed = false;
  this->connectivity = NULL;
  this->connectivity_verts = NULL;
  this->connectivity_colors = NULL;
  this->graphics_conn_geometry = NULL;
  this->show_conn = false;
  this->paired_conn = false;
  this->map_strain = false;
  this->strains = NULL;
  this->map_energy = false;
  this->energy = NULL;
  this->map_scale_set = false;
  this->map_scale_min = 0.0;
  this->map_scale_max = 0.0;
  this->max_distance = 0.0;
  this->ljspring = false;
  }

//*============================================================*
//*==========            getNumSprings               ==========*
//*============================================================*
// get number of springs.                                  

void
PmSpringPotential::getNumSprings(int& num) {
  num = number_springs;
  }

//*============================================================*
//*==========            getNumInteractions          ==========*
//*============================================================*
// get number of interactions (springs).

int 
PmSpringPotential::getNumInteractions()
  {
  if (!this->connectivity) {
    compConnectivity();
    }

  return number_springs;
  }

//*============================================================*
//*==========      setPairedConnectivity             ==========*
//*============================================================*
// set the paired connectivity for spring potentials.

void 
PmSpringPotential::setPairedConnectivity(const bool flag) {
  paired_conn = flag;
  }

//*============================================================*
//*==========           setLJSpring                  ==========*
//*============================================================*
// set Lennard-Jones spring option.                    

void 
PmSpringPotential::setLJSpring(const bool val)
  {
  ljspring = val;
  }

void
PmSpringPotential::getLJSpring(bool& val)
  {
  val = ljspring;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// update the geometry for a spring potential.

void
PmSpringPotential::updateDisplay()
  {
  #ifdef dbg_PmSpringPotential_updateDisplay
  fprintf (stderr, ">>> PmSpringPotential:: updateDisplay[%s] \n", name.c_str());
  fprintf (stderr, ">>> show conn= %d \n", this->show_conn); 
  #endif

  if (this->show_conn) {
    displayConnectivity();
    }

  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces between two spring potentials.

void 
PmSpringPotential::compForces(const bool upd, const float time, const bool active,
                              const bool verbose)
  {
  #define ndbg_PmSpringPotential_compForces
  #ifdef dbg_PmSpringPotential_compForces
  fprintf (stderr, ">>>>>> PmSpringPotential: comp forces name=%s time=%f \n", 
           name.c_str(), time);
  fprintf (stderr, ">>> active=%d \n", active); 
  fprintf (stderr, ">>> force const=%f \n", this->force_const); 
  #endif

  if (!active) {
    this->current_energy = 0.0;
    this->current_strain = 0.0;
    return;
    }

  this->verbose = verbose;

  if (!this->connectivity_computed) {
    compConnectivity();
    }

  if (number_springs == 0) {
    return;
    }

  //===== compute forces =====//

  int i, j;
  PmPotentialPoints *points1, *points2;
  PmVector3 pos1, pos2, dir, force;
  float df, dist, s, mag, total_energy, total_strain, strain;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  s = this->force_const;
  total_energy = 0.0;
  total_strain = 0.0;

  for (int k = 0; k < number_springs; k++) {
    i = connectivity[2*k];
    j = connectivity[2*k+1];
    pos1 = points1->curr_coordinates[i];
    pos2 = points2->curr_coordinates[j];
    dir = pos2 - pos1;
    dist = dir.length();
    current_distances[k] = dist;

    if (this->ljspring) {
      s = force_consts[k];
      }

    if (distances[k] == 0.0) {
      mag = 0.0;
      df = 0.0;
      strain = 0.0; 
      }
    else {
      df = dist - distances[k];
      mag = s*df; 
      strain = df / distances[k];
      }

    if (dist == 0.0) {
      force.set(0,0,0);
      }
    else {
      force = mag*dir/dist;
      }

    if ((max_strain != 0.0) && (strain > max_strain)) {
      this->strains[k] = strain;
      spring_active[k] = false;
      continue;
      }

    if ((max_distance != 0.0) && (dist > max_distance)) {
      spring_active[k] = false;
      continue;
      }

    spring_active[k] = true;
    this->sobj1->addForce(pos1, force);
    force = -force;
    this->sobj2->addForce(pos2, force);
    this->energy[k] = 0.5*s*df*df;
    total_energy += this->energy[k];

    if (distances[k] == 0.0) {
      this->strains[k] = 0.0;
      }
    else {
      this->strains[k] = strain; 
      }

    total_strain += this->strains[k];
    }

  this->current_energy = total_energy;
  this->current_strain = total_strain;
  }

//*============================================================*
//*==========              setForceConst             ==========*
//*============================================================*
// set the force constant. 

void 
PmSpringPotential::setForceConst(const float val) {
  this->force_const = val;
  #ifdef dbg_PmSpringPotential_setForceConst
  fprintf (stderr, ">>>>>> PmSpringPotential::setForceConst %s \n", name.c_str());
  fprintf (stderr, ">>> force const=%f \n", this->force_const); 
  #endif
  }

void
PmSpringPotential::getForceConst(float& val) 
  {
  if (force_consts.size()) {
    val = 0.0;

    for (unsigned int i = 0; i < force_consts.size(); i++) {
      val += force_consts[i];
      }

    val = val / force_consts.size();
    }
  else {
    val = this->force_const;
    }
  }

//*============================================================*
//*==========             queryState                 ==========*
//*============================================================*

void
PmSpringPotential::queryState(PmQuery& query)
  {
  float strain, dist;

  int i = query.entity-1;

  if (i >= number_springs) {
    return;
    }

  fprintf (stderr, "    >>> strain=%f \n", strains[i]*100.0); 
  fprintf (stderr, "    >>> energy=%f \n", energy[i]); 
  fprintf (stderr, "    >>> init distance=%f %\n", distances[i]); 
  fprintf (stderr, "    >>> current distance=%f %\n", current_distances[i]); 
  }

//*============================================================*
//*==========              getStrain                 ==========*
//*============================================================*
// get the strain at a particular spring.

void
PmSpringPotential::getStrain(const int i, float& strain) 
  {
  if (i >= number_springs) {
    strain = 0.0;
    return;
    }

  strain = strains[i];
  }

//*============================================================*
//*==========              getDistance               ==========*
//*============================================================*
// get the distance for a particular spring.

void
PmSpringPotential::getDistance(const int i, float& dist)
  {
  if (i >= number_springs) {
    dist = 0.0;
    return;
    }

  dist = distances[i];
  }

//*============================================================*
//*==========              setMaxDistance            ==========*
//*============================================================*

void 
PmSpringPotential::setMaxDistance(const float val)
  {
  max_distance = val;
  }

//*============================================================*
//*==========              setCutoff                 ==========*
//*============================================================*
// set the cutoff for spring connectivity.

void
PmSpringPotential::setCutoff(const float val) {
  this->cutoff = val;
  }

//*============================================================*
//*==========              compConnectivity          ==========*
//*============================================================*
// compute the spring connectivity.

void
PmSpringPotential::compConnectivity()
  {
  #define ndbg_PmSpringPotential_compConnectivity
  #ifdef dbg_PmSpringPotential_compConnectivity
  fprintf (stderr, "\n>>>>>> PmSpringPotential::compConnectivity \n");
  fprintf (stderr, ">>> cutoff = %f \n", cutoff);
  fprintf (stderr, ">>> paired conn = %d \n", paired_conn);
  #endif

  PmPotentialPoints *points1, *points2;
  float dist, r1, r2;
  PmVector3 pos1, pos2, dir;
  vector<int> conn;
  vector<float> dist0;
  bool verbose = pmSystem.getCmdVerbose();

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  number_springs = 0;
  connectivity_computed = true;

  #ifdef dbg_PmSpringPotential_compConnectivity
  fprintf (stderr, ">>> num pts1 %d \n", points1->number); 
  fprintf (stderr, ">>> num pts2 %d \n", points2->number); 
  #endif

  // create springs only between pairs of points. //
  // this is for hydrogen and ionic bonds.        //

  if (paired_conn) {
    if (points1->number != points2->number) {
      pm_ErrorWarnReport (PM, "paired-potential \"%s\" doesn't have equal points.", 
                          "*", name.c_str());
      return;
      }

    for (int i = 0; i < points1->number; i++) {
      pos1 = points1->coordinates[i];
      pos2 = points2->coordinates[i];
      dir = pos2 - pos1;
      dist = dir.length();
      dist0.push_back(dist);
      conn.push_back(i);
      conn.push_back(i);
      number_springs += 1;
      }
    }

  //===== connect all points within cutoff =====//

  else {
    PmPhysicalObj *pobj1, *pobj2;
    PmRegion *rgn1, *rgn2;

    if (this->ljspring) {
      points1->simulation_obj->getPhysicalObj(&pobj1);
      points2->simulation_obj->getPhysicalObj(&pobj2);
      pobj1->getRegion(points1->rgn_name, &rgn1);
      pobj2->getRegion(points2->rgn_name, &rgn2);
      //fprintf (stderr, ">>> rgn1 atoms size=%d \n", rgn1->atom_names.size()); 
      //fprintf (stderr, ">>> rgn2 atoms size=%d \n", rgn2->atom_names.size()); 

      if (rgn1->atom_names.empty() || rgn2->atom_names.empty()) {
        pm_ErrorWarnReport (PM, "no atom names for LJ springs.", "*");
        return;
        }
      }

    for (int i = 0; i < points1->number; i++) {
      pos1 = points1->coordinates[i];
      r1 = points1->radii[i];

      for (int j = 0; j < points2->number; j++) {
        pos2 = points2->coordinates[j];
        r2 = points2->radii[j];
        dir = pos2 - pos1;
        dist = dir.length();

        if (dist - r1 - r2 <= cutoff) {

          if (this->ljspring) {
            float k, r0;
            getLJSpringConst (rgn1->atom_names[i], rgn2->atom_names[j], dist, k, r0);

            if ((fabs(dist-r0) / dist) < 0.3) {
              this->force_consts.push_back(k);
              conn.push_back(i);
              conn.push_back(j);
              dist0.push_back(dist);
              number_springs += 1;

              //fprintf (stderr, ">>> %d:%d  r1=%f r2=%f dist=%f \n", i, j, r1, r2, dist);
              if (verbose) {
                fprintf (stderr, ">>> add spring=%d atoms=%s-%s dist=%g k=%f r0=%f \n", 
                         number_springs, rgn1->atom_names[i].c_str(), 
                         rgn2->atom_names[j].c_str(), dist, k, r0); 
                }
              }
            }

          else {
            /*
            fprintf (stderr, ">>> add spring=%d atoms=%s-%s dist=%g \n", 
                     number_springs, rgn1->atom_names[i].c_str(), 
                     rgn2->atom_names[j].c_str(), dist); 
            */
            //fprintf (stderr, ">>> add spring=%d dist=%g \n", number_springs, dist); 
            conn.push_back(i);
            conn.push_back(j);
            dist0.push_back(dist);
            number_springs += 1;
            }
          }
        }
      }
    }

  #ifdef dbg_PmSpringPotential_compConnectivity
  fprintf (stderr, ">>> num springs = %d \n", number_springs); 
  #endif

  if (number_springs == 0) {
    if (!no_connections_warning) {
      pm_ErrorWarnReport (PM, "spring potential \"%s\" has no connections within cutoff.", 
                          "*", name.c_str());
      no_connections_warning = true;
      }
    return;
    }

  connectivity = new int[2*number_springs];
  distances = new float[number_springs];
  current_distances = new float[number_springs];
  strains = new float[number_springs];
  energy = new float[number_springs];
  spring_active = new bool[number_springs];

  for (unsigned int i = 0; i < conn.size(); i++) {
    connectivity[i] = conn[i];
    }

  for (unsigned int i = 0; i < dist0.size(); i++) {
    distances[i] = dist0[i];
    }

  for (int i = 0; i < number_springs; i++) {
    strains[i] = 0.0;
    energy[i] = 0.0;
    spring_active[i] = true;
    }
  }

//*============================================================*
//*==========              getLJSpringConst          ==========*
//*============================================================*
// get a spring constant approximating a Lennard-Jones potential
// for a pair of atoms.

void 
PmSpringPotential::getLJSpringConst (const string atom_name1, const string atom_name2, 
                                     const float dist, float& k, float& r0)
  {
  float sig, eps, sig1, eps1, sig2, eps2;
  float v;
  int n=0;
  sig1 = sig2 = 0.0;
  eps1 = eps2 = 0.0;

  while (LjData[n].atom_name != "") {
    if (atom_name1[0] == LjData[n].atom_name[0]) {
      sig1 = LjData[n].sig;
      eps1 = LjData[n].eps;
      }

    if (atom_name2[0] == LjData[n].atom_name[0]) {
      sig2 = LjData[n].sig;
      eps2 = LjData[n].eps;
      }

    n += 1;
    }

  sig = (sig1 + sig2) / 2.0;
  eps = sqrt(eps1*eps2);

  // convert kJ / mol to pN.nm //

  eps = eps / 1.654;

  // compute minimum energy distance //

  r0 = pow(2.0, 1.0/6.0) * sig;

  // compute spring constant //

  k = 36.0*eps /(pow(2.0, 2.0/3.0) * sig*sig);

  //v = 4*eps*(pow(sig/dist, 12.0) - pow(sig/dist,6.0));
  //fprintf (stderr, ">>> sig=%f eps=%f energy=%g \n", sig, eps, v);
  }

//*============================================================*
//*==========              setDisplayConn            ==========*
//*============================================================*
// set show the spring connectivity.

void
PmSpringPotential::setDisplayConn(const bool flag) {
  show_conn = flag;
  updateDisplay();
  }

//*============================================================*
//*==========              setMapEnergy              ==========*
//*============================================================*
// set the flag to map energy into colors.

void
PmSpringPotential::setMapEnergy(const bool flag) {
  map_energy = flag;
  }

//*============================================================*
//*==========              setMapScale               ==========*
//*============================================================*
// set the map scale.

void 
PmSpringPotential::setMapScale(const float min, const float max)
  {
  this->map_scale_set = true;
  this->map_scale_min = min;
  this->map_scale_max = max;
  }

//*============================================================*
//*==========              setMapStrain              ==========*
//*============================================================*
// set the flag to map strain into colors.

void
PmSpringPotential::setMapStrain(const bool flag) {
  map_strain = flag;
  }

//*============================================================*
//*==========              displayConn               ==========*
//*============================================================*
// show the spring connectivity.

void
PmSpringPotential::displayConnectivity()
  {
  //fprintf (stderr, ">>>>>> PmSpringPotential::displayConnectivity \n");
  //fprintf (stderr, ">>> number_springs=%d\n", number_springs);
  int i, j;
  PmPotentialPoints *points1, *points2;
  PmVector3 pos1, pos2;  

  if (!this->connectivity_computed) {
    compConnectivity();
    }

  if (number_springs == 0) {
    return;
    }

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);

  if (!graphics_conn_geometry) { 
    connectivity_verts = new PmVector3[2*number_springs];
    connectivity_colors = new PmVector3[number_springs];
    }

  if (this->map_strain) {
    float vmin, vmax;
    vmin = vmax = strains[0];

    if (map_scale_set) {
      vmin = map_scale_min;
      vmax = map_scale_max;
      }
    else {
      for (int i = 0; i < number_springs; i++) {
        if (strains[i] < vmin) vmin = strains[i];
        if (strains[i] > vmax) vmax = strains[i];
        }
      }

    PmGraphicsInterface::mapDataToColors(vmin, vmax, number_springs, strains, 
                                         connectivity_colors);

    /*
    vmin = vmax = strains[0];
    for (int i = 0; i < number_springs; i++) {
      if (strains[i] < vmin) vmin = strains[i];
      if (strains[i] > vmax) vmax = strains[i];
      }
    fprintf (stderr, ">>> vmin=%f  vmax=%f \n", vmin,vmax); 
    */
    }

  if (this->map_energy) {
    float vmin, vmax;

    if (map_scale_set) {
      vmin = map_scale_min;
      vmax = map_scale_max;
      }
    else {
      vmin = vmax = energy[0];

      for (int i = 0; i < number_springs; i++) {
        if (energy[i] < vmin) vmin = energy[i];
        if (energy[i] > vmax) vmax = energy[i];
        }
      }

    PmGraphicsInterface::mapDataToColors(vmin, vmax, number_springs, energy,
                                         connectivity_colors);
    }

  for (int k = 0; k < number_springs; k++) {
    i = connectivity[2*k];
    j = connectivity[2*k+1];
    connectivity_verts[2*k] = points1->curr_coordinates[i];
    connectivity_verts[2*k+1] = points2->curr_coordinates[j];

    if (!spring_active[k]) {
      connectivity_colors[k].set(0,0,0);
      }
    }

  if (!graphics_conn_geometry) { 
    string geom_name = "potential";
    geom_name = geom_name + '[' + name + ']' + "spring_conn";
    //fprintf (stderr, ">>> geom_name=%s\n", geom_name.c_str());

    PmGraphicsAttributes atts;
    atts.setDisjoint(true);
    atts.setColor(this->color); 
    atts.setLineWidth(this->line_width); 

    PmGraphicsLine *geom = new PmGraphicsLine(geom_name, 2*number_springs, 
                                              connectivity_verts);
    if (this->map_strain || this->map_energy) {
      geom->setColors(number_springs, connectivity_colors);
      }

    geom->setAttributes(atts);
    geom->display();
    graphics_conn_geometry = geom;
    }

  else {
    PmGraphicsLine *geom = dynamic_cast<PmGraphicsLine*>(graphics_conn_geometry);
    geom->update(2*number_springs, connectivity_verts);

    if (this->map_strain || this->map_energy) {
      geom->updateColors(number_springs, connectivity_colors);
      }
    }
  }

////////////////////////////////////////////////////////////////
//           b e n d         p o t e n t i a l               //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmBendPotential::PmBendPotential(const string name)
  {
  //fprintf (stderr, ">>> PmBendPotential: ctor [%s] \n", name.c_str());
  this->has_origin = false;
  this->name = name;
  this->init = false;
  this->force_const = 1.0;
  //this->connectivity = NULL;
  this->graphics_conn_geometry = NULL;
  this->show_conn = false;
  this->conn_color.set(1,1,1);
  this->conn_line_width = 1.0;
  }

//*============================================================*
//*==========               setPoint                 ==========*
//*============================================================*
// set spring point.                           

void 
PmBendPotential::setPoint(const PmVector3& point) 
  {
/*
  this->point = point;
  PmPotentialPoints *pgeom = dynamic_cast<PmPotentialPoints*>(this->geometry1);

  if (pgeom->number == 1) {
    pgeom->number = 2;
    pgeom->coordinates[1] = point; 
    pgeom->curr_coordinates[1] = point; 
    }
  else {
    pgeom->number = 1;
    pgeom->center = point;
    pgeom->curr_center = point;
    pgeom->coordinates[0] = point; 
    pgeom->curr_coordinates[0] = point; 
    }
*/
  }

//*============================================================*
//*==========               setOrigin                ==========*
//*============================================================*
// set spring origin.                           

void 
PmBendPotential::setOrigin(const PmVector3& origin) 
  {
/*
  this->origin = origin;
  has_origin = true;
  PmPotentialPoints *pgeom = dynamic_cast<PmPotentialPoints*>(this->geometry1);

  if (pgeom->number == 1) {
    pgeom->number = 2;
    pgeom->coordinates[1] = pgeom->coordinates[0];
    pgeom->curr_coordinates[1] = pgeom->curr_coordinates[0];
    }
  else {
    pgeom->number = 1;
    }

  pgeom->center = origin;
  pgeom->curr_center = origin;
  pgeom->coordinates[0] = origin;
  pgeom->curr_coordinates[0] = origin;
*/
  }

//*============================================================*
//*==========            getNumInteractions          ==========*
//*============================================================*
// get number of interactions (springs).

int 
PmBendPotential::getNumInteractions() {
  return 1;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// compute the forces between two spring potentials.

void
PmBendPotential::updateDisplay()
  {
  //fprintf (stderr, ">>> PmBendPotential:: updateDisplay[%s] \n", name.c_str());

  if (this->show_conn) {
    displayConnectivity();
    }

  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces between two spring potentials.

void 
PmBendPotential::compForces(const bool upd, const float time, const bool active,
                              const bool verbose)
  {
  /*
  fprintf (stderr, ">>> PmBendPotential: comp forces [%s] \n", name.c_str());
  fprintf (stderr, "   >>> num pot[%d] \n", potentials.size());
  */
  this->verbose = verbose;

  if (!active) {
    return;
    }

  //===== compute forces =====//

  PmPotentialPoints *points1, *points2;
  PmVector3 pt1, pt2, origin, dir, force;
  PmVector3 v1, v2; 
  float ang, s, dang; 

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  s = this->force_const;

  if (!init) {
    initConfig();
    }

  if (points1->number == 2) {
    pt1 = points1->curr_coordinates[0];
    origin = points1->curr_coordinates[1];
    pt2 = points2->curr_coordinates[0];
    }
  else {
    pt1 = points2->curr_coordinates[0];
    origin = points2->curr_coordinates[1];
    pt2 = points1->curr_coordinates[0];
    }

  v1 = pt1 - origin;
  v2 = pt2 - origin;

  v1.normalize();
  v2.normalize();
  dir = v1.cross(v2);
  ang = asin(dir.length());
  dang = ang - init_angle;
  force = -s*dang*dir;
  //fprintf (stderr, ">>> dang=%f \n", dang);
  //fprintf (stderr, ">>> s=%f \n", s); 
  //fprintf (stderr, ">>> force=%f %f %f\n", force[0], force[1], force[2]); 

  this->sobj1->addTorque(this->origin, force);
  force = -force;
  this->sobj2->addTorque(this->origin, force);

  current_energy = 0.5*s*dang*dang;
  }

//*============================================================*
//*==========              setForceConst             ==========*
//*============================================================*
// set the force constant. 

void 
PmBendPotential::setForceConst(const float val) {
  this->force_const = val;
  }

//*============================================================*
//*==========              setDisplayConn            ==========*
//*============================================================*
// set show the spring connectivity.

void
PmBendPotential::setDisplayConn(const bool flag) {
  show_conn = flag;
  }

//*============================================================*
//*==========              displayConn               ==========*
//*============================================================*
// show the spring connectivity.

void
PmBendPotential::displayConnectivity()
  {
  //fprintf (stderr, ">>>>>> PmBendPotential::displayConnectivity \n");
  int n;
  PmPotentialPoints *points1, *points2;
  PmVector3 pos1, pos2;  

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  n = points1->number + points2->number;
  //fprintf (stderr, ">>> n=%d \n", n);

  if (!graphics_conn_geometry) { 
    connectivity_verts = new PmVector3[n];
    }

  connectivity_verts[0] = points1->curr_coordinates[0];
  connectivity_verts[1] = points1->curr_coordinates[1];
  connectivity_verts[2] = points2->curr_coordinates[0];

  if (!graphics_conn_geometry) { 
    string geom_name = "potential";
    geom_name = geom_name + '[' + name + ']' + "angle_conn";

    PmGraphicsAttributes atts;
    atts.setColor(this->geometry1->color); 
    atts.setLineWidth(this->geometry1->line_width); 

    PmGraphicsLine *geom = new PmGraphicsLine(geom_name, 3, connectivity_verts);
    geom->setAttributes(atts);
    geom->display();
    graphics_conn_geometry = geom;
    }

  else {
    PmGraphicsLine *geom = dynamic_cast<PmGraphicsLine*>(graphics_conn_geometry);
    geom->update(3, connectivity_verts);
    }
  }

//*============================================================*
//*==========              compConnectivity          ==========*
//*============================================================*
// compute the spring connectivity.

void
PmBendPotential::compConnectivity(PmBendPotential *pot)
  {
  //fprintf (stderr, "\n>>>>>> PmSpringPotential::compConnectivity \n");
  //PmPotentialPoints *points1, *points2;
  //float dist, r1, r2;
  }

//*============================================================*
//*==========              initConfig                ==========*
//*============================================================*
// initialize the spring configuration.                       

void
PmBendPotential::initConfig()
  {
  #ifdef dbg_PmBendPotential_initConfig
  fprintf (stderr, "\n>>>>>> PmBendPotential::initConfig \n");
  #endif
 
  PmPotentialPoints *points1, *points2;
  PmVector3 pt1, pt2, origin; 
  PmVector3 dir, v1, v2;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);

  #ifdef dbg_PmBendPotential_initConfig
  fprintf (stderr, ">>> points1 number=%d \n", points1->number); 
  fprintf (stderr, ">>> points2 number=%d \n", points2->number); 
  #endif

  if (points1->number == 2) {
    this->point1 = points1->coordinates[0];
    this->origin = points1->coordinates[1];
    this->point2 = points2->coordinates[0];
    }
  else {
    this->point1 = points2->coordinates[0];
    this->origin = points2->coordinates[1];
    this->point2 = points1->coordinates[0];
    }

  v1 = this->point1 - this->origin;
  v2 = this->point2 - this->origin;
  v1.normalize();
  v2.normalize();
  dir = v1.cross(v2);
  init_angle = asin(dir.length());
  init = true;
  }

////////////////////////////////////////////////////////////////
//           a x i s         p o t e n t i a l               //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmAxisPotential::PmAxisPotential(const string name)
  {
  //fprintf (stderr, ">>> PmAxisPotential: ctor [%s] \n", name.c_str());
  this->name = name;
  this->init = false;
  this->restraint_cutoff = 1e6;
  this->restraint_strength = 0.0;
  this->axes_strength[0] = 0.0;
  this->axes_strength[1] = 0.0;
  this->point_strength = 0.0;

  this->graphics_geometry = NULL;
  this->show = false;
  this->color.set(1,1,1);
  this->line_width = 1.0;
  }

//*============================================================*
//*==========              setPointStrength          ==========*
//*============================================================*
// set point strength.

void 
PmAxisPotential::setPointStrength(const float val) {
  point_strength = val;
  }

//*============================================================*
//*==========              setRestStrength           ==========*
//*============================================================*
// set restraint strength.

void 
PmAxisPotential::setRestStrength(const float val)
  {
  this->restraint_strength = val;
  }

//*============================================================*
//*==========              setRestCutoff             ==========*
//*============================================================*
// set restraint cutoff.

void 
PmAxisPotential::setRestCutoff(const float val)
  {
  this->restraint_cutoff = val;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// compute the forces between two spring potentials.

void
PmAxisPotential::updateDisplay()
  {
  //fprintf (stderr, ">>> PmAxisPotential:: updateDisplay[%s] \n", name.c_str());

  if (this->show) {
    displayConnectivity();
    }

  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces for axis potentials.

void 
PmAxisPotential::compForces(const bool upd, const float time, const bool active,
                            const bool verbose)
  {
  #ifdef dbg_PmAxisPotential_compForces
  fprintf (stderr, ">>>>>> PmAxisPotential: comp forces [%s] \n", name.c_str());
  fprintf (stderr, ">>> restraint_cutoff=%f \n", this->restraint_cutoff); 
  fprintf (stderr, ">>> restraint_strength=%f \n", this->restraint_strength); 
  #endif
  this->verbose = verbose;

  if (!active) {
    return;
    }

  PmPotentialLines *lines1 = dynamic_cast<PmPotentialLines*>(geometry1);
  PmPotentialLines *lines2 = dynamic_cast<PmPotentialLines*>(geometry2);

  // compute restraint force between com //

  if (restraint_strength != 0.0) {
    PmVector3 pos1, pos2, dir, force;
    float dist, dist0, s, fmag;

    pos1 = lines1->origin;
    pos2 = lines2->origin;
    dir = pos2 - pos1; 
    dist0 = dir.length();

    pos1 = lines1->curr_origin;
    pos2 = lines2->curr_origin;
    dir = pos2 - pos1; 
    dist = dir.length();

    if (dist > this->restraint_cutoff) {
      s = this->restraint_strength;

      if (dist0 == 0.0) {
        fmag = 0.0;
        }
      else {
        fmag = s*dist / dist0;
        }

      //fprintf (stderr, ">>> dist = %f  fmag = %f \n", dist, fmag);
  
      if (dist == 0.0) {
        force.set(0,0,0);
        }
      else {
        force = fmag*dir/dist;
        }

      this->sobj1->addForce(pos1, force);
      force = -force;
      this->sobj2->addForce(pos2, force);
      }
    }

  // compute forces between axes //

  int id1, id2, j;
  PmVector3 dir, v, v1, v2, force;
  float s, ang;

  #ifdef dbg_PmAxisPotential_compForces
  fprintf (stderr, ">>> lines1->num_axes=%d \n", lines1->num_axes); 
  fprintf (stderr, ">>> lines2->num_axes=%d \n", lines2->num_axes); 
  #endif

  for (int i = 0; i < lines1->num_axes; i++) {
    s = this->axes_strength[i];
    id1 = lines1->axes_ids[i];
    v1 = lines1->curr_coordinates[2*id1+1] - lines1->curr_coordinates[2*id1]; 
    id2 = lines2->axes_ids[i];
    v2 = lines2->curr_coordinates[2*id2+1] - lines2->curr_coordinates[2*id2]; 

    v1.normalize();
    v2.normalize();
    //dir = v2.cross(v1);
    dir = v1.cross(v2);
    ang = asin(dir.length());

    #ifdef dbg_PmAxisPotential_compForces
    fprintf (stderr, ">>> v1 = %f %f %f \n", v1[0], v1[1], v1[2]);
    fprintf (stderr, ">>> v2 = %f %f %f \n", v2[0], v2[1], v2[2]);
    fprintf (stderr, ">>> s = %f  \n", s);
    fprintf (stderr, ">>> ang = %f  \n", ang);
    #endif

    dir.normalize();
    force = s*ang*dir;
    this->sobj1->addTorque(lines1->origin, force);
    force = -force;
    this->sobj2->addTorque(lines2->origin, force);
    }

  // compute forces between points //

  if (this->point_strength != 0.0) {
    int id1, id2, j;
    PmVector3 dir, v, v1, v2, force;
    PmVector3 pos1, pos2; 
    float dist, fmag, s;
    s = this->point_strength;

    for (int i = 0; i < 4; i++) {
      if (i == 3) {
        pos1 = lines1->curr_origin;
        pos2 = lines2->curr_origin;
        }
      else {
        pos1 = lines1->curr_coordinates[2*i+1];
        pos2 = lines2->curr_coordinates[2*i+1];
        }

      dir = pos2 - pos1; 
      dist = dir.length();
      dir.normalize();
      fmag = s;
      //fprintf (stderr, ">>> fmag = %f  \n", fmag);

      force = fmag*dir;
      this->sobj1->addForce(pos1, force);
      force = -force;
      this->sobj2->addForce(pos2, force);
      }
    } 
  }

//*============================================================*
//*==========              setAxesStrength           ==========*
//*============================================================*
// set the axes potential strength.

void 
PmAxisPotential::setAxesStrength(const int id, const float val) 
  {
  this->axes_strength[id-1] = val;
  }

void 
PmAxisPotential::setAxesStrength(const vector<float>& vals)
  {
  for (unsigned int i = 0; i < vals.size(); i++) {
    axes_strength[i] = vals[i];
    }
  }

//*============================================================*
//*==========              setDisplayConn            ==========*
//*============================================================*
// set show the potentail connectivity.

void
PmAxisPotential::setDisplayConn (const bool flag) {
  show = flag;
  }

//*============================================================*
//*==========              displayConnectivity       ==========*
//*============================================================*

void
PmAxisPotential::displayConnectivity()
  {
  }

////////////////////////////////////////////////////////////////
//  m o l e c u l a r  m e c h a n i c s   p o t e n t i a l //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmMolecularMechanicsPotential::PmMolecularMechanicsPotential(const string name)
  {
  //fprintf (stderr, ">>> PmMolecularMechanicsPotential: ctor [%s] \n", name.c_str());
  this->name = name;
  this->initialized = false;
  this->maximum_force = 100.0;
  this->strength = 1.0;
  this->graphics_geometry = NULL;
  this->current_energy = 0.0;
  this->num_contact = 0;
  }

//*============================================================*
//*==========              setCutoff                 ==========*
//*============================================================*
// set the cutoff for calculating energy   

void
PmMolecularMechanicsPotential::setCutoff(const float val) {
  this->cutoff = val;
  }

//*============================================================*
//*==========              setStrength               ==========*
//*============================================================*
// set stength.

void
PmMolecularMechanicsPotential::setStrength(const float val)
  {
  this->strength = val;
  }

//*============================================================*
//*==========              setTerms                  ==========*
//*============================================================*
// set terms.                            

void
PmMolecularMechanicsPotential::setTerms(const PmPotentialMolMechTerms& terms)
  {
  this->terms = terms;  
  }

//*============================================================*
//*==========              getParameters             ==========*
//*============================================================*
// get parameters.

void
PmMolecularMechanicsPotential::getParameters(PmPotentialMolMechParameters& params)
  {
  params = this->parameters;
  }

//*============================================================*
//*==========              setParameters             ==========*
//*============================================================*
// set parameters.

void 
PmMolecularMechanicsPotential::setParameters(const PmPotentialMolMechParameters& params)
  {
  this->parameters = params;
  }

//*============================================================*
//*==========              setMaxForce               ==========*
//*============================================================*

void 
PmMolecularMechanicsPotential::setMaxForce(const float val)
  {
  maximum_force = val;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*

void
PmMolecularMechanicsPotential::updateDisplay()
  {
  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              queryState                ==========*
//*============================================================*

void
PmMolecularMechanicsPotential::queryState(PmQuery& query)
  {
  fprintf (stderr, "    >>> num contact=%d %\n", num_contact); 
  }

//*============================================================*
//*==========              compCoulombForces         ==========*
//*============================================================*
// compute the Coulomb forces between potentials.

void
PmMolecularMechanicsPotential::compCoulombForces(const bool upd, const float time,
                                                 const bool active, const bool verbose)
  {
  #ifdef dbg_PmMolecularMechanicsPotential_compCoulombForces 
  fprintf (stderr,"\n>>>>>> PmMolecularMechanicsPotential::compCoulombForces \n"); 
  fprintf (stderr, ">>> name=%s \n", name.c_str());
  fprintf (stderr, ">>> cutoff=%f \n", this->cutoff);
  #endif

  if (!active) {
    return;
    }

  PmPotentialPoints *points1, *points2;
  float r1, r2, dist, q1, q2, fmag, total_force, charge_scale;
  float charge1, charge2;
  PmVector3 pos1, pos2, dir, force;
  bool debug = false;
  int num_int;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points1->contact = false;
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  points2->contact = false;

  //===== do a quick extent overlap test =====//

  dir = points1->curr_center - points2->curr_center;
  dist = dir.length();

  if (dist > this->cutoff + points2->radius) {
    #ifdef dbg_PmMolecularMechanicsPotential_compCoulombForces 
    fprintf (stderr,">>> *** no overlap *** \n"); 
    #endif
    return;
    }

  //===== initiaize charge data =====//

  if (charge_data1.empty()) {
    this->initChargeData();
    }

  //===== compute forces =====//

  num_int = 0;
  total_force = 0.0;
  charge_scale = this->parameters.charge_scale;

  for (int i = 0; i < points1->number; i++) {
    r1 = points1->radii[i];
    pos1 = points1->curr_coordinates[i];
    q1 = charge_scale*this->charge_data1[i];

    for (int j = 0; j < points2->number; j++) {
      r2 = points2->radii[j];
      pos2 = points2->curr_coordinates[j];
      dir = pos1 - pos2;
      dist = dir.length();

      if (dist - r1 - r2 <= this->cutoff) {
        q2 = this->charge_data2[j];
        fmag = -q1*q2 / (dist*dist);
        num_int += 1;
        dir.normalize();
        //fprintf (stderr, " q1 %f  q2 %f \n", q1, q2); 

        if (fabs(fmag) > this->maximum_force) {
          fmag = this->maximum_force * pm_MathSign(fmag);
          }

        //fprintf (stderr, ">>> fmag = %g  \n", fmag); 
        force = fmag*dir;
        this->sobj2->addForce(pos2, force);

        force = -force;
        this->sobj1->addForce(pos1, force);
        total_force += fmag;
        }
      }
    }
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces between potentials.

void
PmMolecularMechanicsPotential::compForces(const bool upd, const float time,
                                          const bool active,
                                          const bool verbose)
  {
  //#define dbg_PmMolecularMechanicsPotential_compForces
  #ifdef dbg_PmMolecularMechanicsPotential_compForces
  fprintf (stderr, "\n>>>>>> PmMolecularMechanicsPotential::compForces \n");
  fprintf (stderr, ">>> name=%s \n", name.c_str()); 
  fprintf (stderr, ">>> max force=%f \n", maximum_force); 
  fprintf (stderr, ">>> terms.lj_repulsive=%d \n", terms.lj_repulsive); 
  #endif

  this->verbose = verbose;

  //===== initialize potential data =====//

  if (!initialized) {
    initPotential();
    }

  if (!active) {
    return;
    }

  //===== calculate coulomb forces =====//

  if (terms.coulomb) {
    this->compCoulombForces(upd, time, active, verbose);
    }

  if (!terms.lj_repulsive) {
    return;
    }

  //===== compute repulsive forces =====//

  PmPotentialPoints *points1, *points2;
  float r1, r2, dist, s, sigma, rp, fmag;
  PmVector3 pos1, pos2, dir, force;
  int id;
  PmMathContactEllipsoids ellipsoids;
  bool contact;
  bool debug = false;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points1->contact = false;
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  points2->contact = false;
  s = this->strength; 
  sigma = 0.34;
  this->current_energy = 0.0;

  #ifdef dbg_PmMolecularMechanicsPotential_compForces
  fprintf (stderr, ">>> lj_rp=%f \n", parameters.lj_repulsive_power); 
  fprintf (stderr, ">>> sigma=%f \n", sigma); 
  fprintf (stderr, ">>> s=%f \n", s); 
  #endif

  //===== do a quick extent overlap test  =====//

  dir = points1->curr_center - points2->curr_center;
  dist = dir.length();

  if (dist > points1->radius + points2->radius) {
    #ifdef dbg_PmMolecularMechanicsPotential_compForces
    fprintf (stderr, ">>> *** no sphere overlap *** \n");
    #endif
    return;
    }

  //===== check for overlap of ellipsoid extents =====//

  if (points1->number > 10) { 
    points1->getCurrentAxes(ellipsoids.pos1, ellipsoids.u1, ellipsoids.v1, ellipsoids.w1);
    points1->getPcaScales(ellipsoids.r1); 
    points2->getCurrentAxes(ellipsoids.pos2, ellipsoids.u2, ellipsoids.v2, ellipsoids.w2);
    points2->getPcaScales(ellipsoids.r2); 

    for (int i = 0; i < 3; i++) {
      ellipsoids.r1[i] += 0.4; 
      ellipsoids.r2[i] += 0.4; 
      }

    #ifdef dbg_PmMolecularMechanicsPotential_compForces
    fprintf (stderr, ">>> check ellipsoid overlap \n");
    fprintf (stderr, "    pos1 = %f %f %f         \n", ellipsoids.pos1[0], 
             ellipsoids.pos1[1], ellipsoids.pos1[2]);
    fprintf (stderr, "    u1 = %f %f %f         \n", ellipsoids.u1[0], 
             ellipsoids.u1[1], ellipsoids.u1[2]);
    fprintf (stderr, "    v1 = %f %f %f         \n", ellipsoids.v1[0], 
             ellipsoids.v1[1], ellipsoids.v1[2]);
    fprintf (stderr, "    w1 = %f %f %f         \n", ellipsoids.w1[0], 
             ellipsoids.w1[1], ellipsoids.w1[2]);
    #endif

    pm_MathEcfCompContact (ellipsoids, contact);

    if (!contact) {
      if (debug) fprintf (stderr, "     no ellipsoid overlap \n");
      #ifdef dbg_PmMolecularMechanicsPotential_compForces
      fprintf (stderr, "    **** no overlap \n");
      #endif
      return;
      }
    }

  //===== brute force check of all points =====//

  #define nuse_brute_PmMolecularMechanicsPotential_compForces
  #ifdef use_brute_PmMolecularMechanicsPotential_compForces
  //bool test_brute = true;
  bool test_brute = false;
  bool brute_contact;
  int num_brute_contact = 0;
  brute_contact = false;

  for (int i = 0; i < points1->number; i++) {
    r1 = points1->radii[i];
    pos1 = points1->curr_coordinates[i];

    for (int j = 0; j < points2->number; j++) {
      r2 = points2->radii[j];
      pos2 = points2->curr_coordinates[j];
      dir = pos2 - pos1; 
      dist = dir.length();

      if (dist < r1 + r2) {
        num_brute_contact += 1;
        brute_contact = true;
        }
      }
    }
  #endif

  //===== use spatial decomposition =====//

  spatial_decomp->setXform (this->geometry1->potential_xform);
  vector<int> ids;
  int num_sd_contact = 0;
  float min_dist = 1e6;
  float lj_rp = parameters.lj_repulsive_power;
  this->current_energy = 0.0;

  for (int i = 0; i < points2->number; i++) {
    r2 = points2->radii[i];
    pos2 = points2->curr_coordinates[i];
    spatial_decomp->getPoints(pos2, r2, ids);

    if (ids.size() != 0) {
      for (unsigned int j = 0; j < ids.size(); j++) {
        id = ids[j];
        r1 = points1->radii[id];
        pos1 = points1->curr_coordinates[id];
        dir = pos2 - pos1;
        dist = dir.length();

        if (dist < r1 + r2) {
          points1->contact = true;
          points2->contact = true;
          rp = lj_rp * pow(sigma/dist, lj_rp) / dist;
          fmag = (s*rp) / dist;
          //fprintf (stderr, " %d %d  fmag = %g \n", i, id, fmag);
          //fprintf (stderr, "#### energy = %g  s=%g  rp=%g \n", s*rp, s, rp); 
          this->current_energy += s*rp;

          if (dist  < min_dist) { 
            min_dist = dist; 
            }

          if (fmag > this->maximum_force) { 
            /*
            fprintf (stderr, "**** WARNING: potential \"%s\" force is large. \n", 
                     name.c_str());
            fprintf (stderr, "     fmag = %g \n", fmag);
            */
            fmag = this->maximum_force;
            dir.normalize();
            }

          force = fmag * dir;
          this->sobj2->addForce(pos2, force);

          force = -force;
          this->sobj1->addForce(pos1, force);
          num_sd_contact += 1;
          //fprintf (stderr, ">>> %s c12 = %f \n", name.c_str(), c12);
          }
        }
      }
    }

  /*
  fprintf (stderr, ">>> num_brute_contact=%d   num_sd_contact=%d \n", num_brute_contact,
           num_sd_contact); 
  */

  #ifdef use_brute_PmMolecularMechanicsPotential_compForces
  if (num_sd_contact != num_brute_contact) {
    fprintf (stderr, "**** WARNING: num_brute_contact = %d but num_sd_contact = %d \n", 
             num_brute_contact, num_sd_contact);
    }
  #endif

  #ifdef dbg_PmMolecularMechanicsPotential_compForces
  fprintf (stderr, ">>> num_sd_contact = %d \n", num_sd_contact);
  #endif

  /*
  if (num_sd_contact) {
    fprintf (stderr,"\n>>> MM: comp forces %s -> %s ", name.c_str(), pname.c_str());
    fprintf (stderr," #contact = %d ", num_sd_contact);
    fprintf (stderr," min dist = %f \n", min_dist);
    }
  */

  if ((time == 0.0) && (num_sd_contact != 0)) {
    if (num_sd_contact) {
      pm_ErrorWarnReport (PM, "molecular mechanics potential %s has %d atoms in initial contact.", "*", name.c_str(), num_sd_contact); 
      }
    }

  this->num_contact = num_sd_contact;
  }

//*============================================================*
//*==========              initPotential             ==========*
//*============================================================*
// initialize some data structures used by the potential.

void
PmMolecularMechanicsPotential::initPotential()
  {
  PmPotentialPoints *points = dynamic_cast<PmPotentialPoints*>(geometry1);

  // create spatial decomposition  //

  this->initSpatialDecomp(points, &spatial_decomp, 0.0);
  initialized = true;
  }

//*============================================================*
//*==========              initChargeData            ==========*
//*============================================================*
// initialize charge data.                                         

void 
PmMolecularMechanicsPotential::initChargeData()
  {
  PmPotentialPoints *points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  PmPotentialPoints *points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  float q1 = geometry1->properties.charge; 
  float q2 = geometry2->properties.charge; 

  if (this->parameters.surface_charge) { 
/*
    vector<int> ids;
    this->sobj1->getRegionIndices(region, ids);
    this->sobj1->getRegionData(region, charge_data);

    if (charge_data.empty()) {
      pm_ErrorWarnReport (PM, "potential \"%s\" does not have surface charge data", "*", 
                          name.c_str());
      for (int i = 0; i < points->number; i++) {
        charge_data.push_back(0.0);
        }
      }

*/
    }
  else {
    for (int i = 0; i < points1->number; i++) {
      charge_data1.push_back(q1);
      }

    for (int i = 0; i < points2->number; i++) {
      charge_data2.push_back(q2);
      }
    }
  }

////////////////////////////////////////////////////////////////
//            molecular mechanics terms functions            //
//////////////////////////////////////////////////////////////

bool 
PmPotentialMolMechTerms::setTerm(const string term)
  {
  if (term == "lj-repulsive") {
    this->lj_repulsive = true;
    return true;
    }

  if (term == "lj-attractive") {
    this->lj_attractive = true;
    return true;
    }

  if (term == "coulomb") {
    this->coulomb = true;
    return true;
    }

  return false;
  }

////////////////////////////////////////////////////////////////
//        molecular mechanics parameter functions            //
//////////////////////////////////////////////////////////////

PmPotentialMolMechParameters::PmPotentialMolMechParameters()
  {
  lj_repulsive_power = 10.0;
  lj_attractive_power = 6.0;
  charge = 0.0;
  surface_charge = false;
  charge_scale = 1.0;
  }

bool 
PmPotentialMolMechParameters::setParameter(const string param)
  {
  string pname, value;
  int epos;
  epos = param.find('=');

  if (epos == -1) {
    return false;
    }

  // parse name=value expression //

  for (int i = 0; i < epos; i++) {
    pname.push_back(param[i]);
    }

  for (unsigned int i = epos+1; i < param.size(); i++) {
    value.push_back(param[i]);
    }

  // lj-power=<int>-<int> //

  if (pname == "lj-power") {
    char val[100];
    int n;
    float rp, ap;
    epos = value.find('-');

    if (epos == -1) {
      return false;
      }

    n = value.size();

    for (int i = 0; i < n; i++) {
      val[i] = value[i];
      }

    val[n] = '\n';
    sscanf (val, "%f-%f", &rp, &ap);
    //fprintf (stderr, ">>> rp = %f  ap = %f \n", rp, ap);
    this->lj_repulsive_power = rp;
    this->lj_attractive_power = ap;
    }

  // charge=<float> //

  else if (pname == "charge") {
    char val[100];
    int n;
    n = value.size();

    if (value == "surface") {
      this->surface_charge = true;
      }
    else {
      for (int i = 0; i < n; i++) {
        val[i] = value[i];
        }

      val[n] = '\n';
      sscanf (val, "%f", &this->charge);
      }

    //fprintf (stderr, ">>> charge = %f \n", charge);
    }

  else if (pname == "charge-scale") {
    char val[100];
    int n;
    n = value.size();

    for (int i = 0; i < n; i++) {
      val[i] = value[i];
      }

    val[n] = '\n';
    sscanf (val, "%f", &this->charge_scale);
    }

  return true;
  }

////////////////////////////////////////////////////////////////
//          t o r s i o n    p o t e n t i a l               //
//////////////////////////////////////////////////////////////
//
//  define points as:
//
//           2'-----1'    geometry2
//           .
//           .
//           .
//           2------1    geometry1

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmTorsionPotential::PmTorsionPotential(const string name)
  {
  //fprintf (stderr, ">>> PmTorsionPotential: ctor [%s] \n", name.c_str());
  this->has_origin = false;
  this->name = name;
  this->init = false;
  this->force_const = 1.0;
  this->graphics_conn_geometry = NULL;
  this->show_conn = false;
  this->conn_color.set(1,1,1);
  this->conn_line_width = 1.0;
  }

//*============================================================*
//*==========               setPoints                ==========*
//*============================================================*
// set spring points.                           

void 
PmTorsionPotential::setPoints(const PmVector3& point1, const PmVector3& point2)
  {

/*
  this->point1 = point1;
  this->point2 = point2;
  PmPotentialPoints *pgeom = dynamic_cast<PmPotentialPoints*>(this->geometry1);

  if (pgeom->number == 1) {
    pgeom->number = 3;
    pgeom->coordinates[1] = point1; 
    pgeom->curr_coordinates[1] = point1; 
    pgeom->coordinates[2] = point2; 
    pgeom->curr_coordinates[2] = point2; 
    }
  else {
    pgeom->number = 2;
    pgeom->center = point1;
    pgeom->curr_center = point1;
    pgeom->coordinates[0] = point1; 
    pgeom->curr_coordinates[0] = point1; 
    pgeom->coordinates[1] = point2; 
    pgeom->curr_coordinates[1] = point2; 
    }
*/
  }

//*============================================================*
//*==========               setOrigin                ==========*
//*============================================================*
// set spring origin.                           

void 
PmTorsionPotential::setOrigin(const PmVector3& origin) 
  {
/*
  this->origin = origin;
  has_origin = true;

  PmPotentialPoints *pgeom = dynamic_cast<PmPotentialPoints*>(this->geometry1);

  if (pgeom->number == 2) {
    pgeom->number = 3;
    pgeom->coordinates[2] = pgeom->coordinates[1];
    pgeom->curr_coordinates[2] = pgeom->curr_coordinates[1];
    pgeom->coordinates[1] = pgeom->coordinates[0];
    pgeom->curr_coordinates[1] = pgeom->curr_coordinates[0];
    }
  else {
    pgeom->number = 1;
    }

  pgeom->center = origin;
  pgeom->curr_center = origin;
  pgeom->coordinates[0] = origin;
  pgeom->curr_coordinates[0] = origin;
*/
  }

//*============================================================*
//*==========            getNumInteractions          ==========*
//*============================================================*
// get number of interactions (springs).

int 
PmTorsionPotential::getNumInteractions() {
  return 1;
  }

//*============================================================*
//*==========              updateDisplay             ==========*
//*============================================================*
// compute the forces between two spring potentials.

void
PmTorsionPotential::updateDisplay()
  {
  //fprintf (stderr, ">>> PmTorsionPotential:: updateDisplay[%s] \n", name.c_str());

  if (this->show_conn) {
    displayConnectivity();
    }

  this->geometry1->updateDisplay();
  this->geometry2->updateDisplay();
  }

//*============================================================*
//*==========              compForces                ==========*
//*============================================================*
// compute the forces between two spring potentials.

void 
PmTorsionPotential::compForces(const bool upd, const float time, const bool active,
                               const bool verbose)
  {
  #ifdef dbg_PmTorsionPotential_compForces
  fprintf (stderr, ">>> PmTorsionPotential: comp forces [%s] \n", name.c_str());
  #endif
  this->verbose = verbose;

  if (!active) {
    return;
    }

  //===== compute forces =====//

  PmPotentialPoints *points1, *points2;
  PmVector3 pt1, pt2, origin, dir, force;
  PmVector3 v1, v2, v3, v12, v23; 
  float dang, ang, s, a1, a2; 

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  s = this->force_const;

  if (!init) {
    initConfig();
    }

  //  define points as:
  //
  //           2'-----1'    geometry2
  //           .
  //           .
  //           2------1    geometry1

  v1 = points1->curr_coordinates[1] - points1->curr_coordinates[0];
  v2 = points2->curr_coordinates[1] - points1->curr_coordinates[1];
  v3 = points2->curr_coordinates[0] - points2->curr_coordinates[1];
  origin = points1->curr_coordinates[2]; 

  v12 = v1.cross(v2);
  v23 = v2.cross(v3);

  a1 = (sqrt(v2*v2)*v1) * v23;
  a2 = v12 * v23;
  ang = atan2(a1,a2);
  dang = (init_angle - ang);
  v2.normalize();

  force = s*dang*v2;
  this->sobj1->addTorque(origin, force);
  force = -force;
  this->sobj2->addTorque(origin, force);

  #ifdef dbg_PmTorsionPotential_compForces
  fprintf (stderr, ">>> dang=%f   s=%f  \n", dang, s);
  #endif

#ifdef test_PmTorsionPotential_compForces
  PmRigidBodyState state;
  float angs1[3], angs2[3], rot[3][3];
  PmBody *body;
  body = (PmBody*)(this->sobj);
  body->getState(state);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rot[i][j] = state.rotation_matrix(i,j);
      }
    }

  pm_MathExtractRotations(rot, angs1);
  //fprintf (stderr, "\n>>> this angs = %f %f %f  \n", angs1[0], angs1[1], angs1[2]);
  body = (PmBody*)(pot->sobj);
  body->getState(state);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rot[i][j] = state.rotation_matrix(i,j);
      }
    }

  pm_MathExtractRotations(rot, angs2);
  //fprintf (stderr, "    pot angs = %f %f %f  \n", angs2[0], angs2[1], angs2[2]);

  for (int i = 0; i < 3; i++) {
    if (i == 0) dir.set(1,0,0);
    if (i == 1) dir.set(0,1,0);
    if (i == 2) dir.set(0,0,1);

    ang = 0.5*(angs1[i] - angs2[i]);

    force = s*ang*dir;
    this->sobj->addTorque(origin, force);

    force = -force;
    pot->sobj->addTorque(origin, force);
    }
#endif
  }

//*============================================================*
//*==========              setForceConst             ==========*
//*============================================================*
// set the force constant. 

void 
PmTorsionPotential::setForceConst(const float val) {
  this->force_const = val;
  }

//*============================================================*
//*==========              setDisplayConn            ==========*
//*============================================================*
// set show the spring connectivity.

void
PmTorsionPotential::setDisplayConn(const bool flag) {
  show_conn = flag;
  }

//*============================================================*
//*==========              displayConn               ==========*
//*============================================================*
// show the spring connectivity.

void
PmTorsionPotential::displayConnectivity()
  {
  //fprintf (stderr, ">>>>>> PmTorsionPotential::displayConnectivity \n");
  int n;
  PmPotentialPoints *points1, *points2;

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);
  n = points1->number + points2->number;

  if (n != 5) {
    pm_ErrorWarnReport (PM, "potential \"%s\" geometry does not have %d points.", "*", 
                        5, name.c_str());
    return;
    }

  if (!graphics_conn_geometry) { 
    connectivity_verts = new PmVector3[n];
    }

  //  define points as:
  // 
  //           2'-----1'    geometry2
  //           .
  //           .
  //           3
  //           |
  //           |
  //           2------1    geometry1

  connectivity_verts[0] = points1->curr_coordinates[0];
  connectivity_verts[1] = points1->curr_coordinates[1];
  connectivity_verts[2] = points1->curr_coordinates[2];
  connectivity_verts[3] = points2->curr_coordinates[1];
  connectivity_verts[4] = points2->curr_coordinates[0];

  if (!graphics_conn_geometry) { 
    string geom_name = "potential";
    geom_name = geom_name + '[' + name + ']' + "angle_conn";

    PmGraphicsAttributes atts;
    atts.setColor(this->geometry1->color); 
    atts.setLineWidth(this->geometry1->line_width); 

    PmGraphicsLine *geom = new PmGraphicsLine(geom_name, n, connectivity_verts);
    geom->setAttributes(atts);
    geom->display();
    graphics_conn_geometry = geom;
    }

  else {
    PmGraphicsLine *geom = dynamic_cast<PmGraphicsLine*>(graphics_conn_geometry);
    geom->update(n, connectivity_verts);
    }
  }

//*============================================================*
//*==========              compConnectivity          ==========*
//*============================================================*
// compute the spring connectivity.

void
PmTorsionPotential::compConnectivity()
  {
  //fprintf (stderr, "\n>>>>>> PmSpringPotential::compConnectivity \n");
  //PmPotentialPoints *points1, *points2;
  //float dist, r1, r2;
  }

//*============================================================*
//*==========              initConfig                ==========*
//*============================================================*
// initialize the spring configuration.                       

void
PmTorsionPotential::initConfig()
  {
  PmPotentialPoints *points1, *points2;
  PmVector3 v1, v2, v3, v12, v23; 
  float ang, a1, a2; 

  points1 = dynamic_cast<PmPotentialPoints*>(geometry1);
  points2 = dynamic_cast<PmPotentialPoints*>(geometry2);

  v1 = points1->curr_coordinates[1] - points1->curr_coordinates[0];
  v2 = points2->curr_coordinates[1] - points1->curr_coordinates[1];
  v3 = points2->curr_coordinates[0] - points2->curr_coordinates[1];

  v12 = v1.cross(v2);
  v23 = v2.cross(v3);

  a1 = (sqrt(v2*v2)*v1) * v23;
  a2 = v12 * v23;
  init_angle = atan2(a1,a2);
  init = true;
  }

}


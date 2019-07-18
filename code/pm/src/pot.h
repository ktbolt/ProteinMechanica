
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
//* pot:                 p o t e n t i al                      *
//*============================================================*

#ifndef _POTENTIAL_PM_H_ 
#define _POTENTIAL_PM_H_

#include "pm/pm.h"
#include "pm/mth.h"
#include "graphics.h"
#include "sd.h"

namespace ProteinMechanica {

class PmSimulationObj;
class PmGraphicsGeometry;

typedef enum PmDynamicsInterForceGeomType {
  PM_POTENTIAL_GEOM_UNKNOWN,
  PM_POTENTIAL_GEOM_CYLINDER,
  PM_POTENTIAL_GEOM_ELLIPSOID,
  PM_POTENTIAL_GEOM_LINES,
  PM_POTENTIAL_GEOM_PLANE,
  PM_POTENTIAL_GEOM_POINTS,
  PM_POTENTIAL_GEOM_SPHERE,
  PM_POTENTIAL_GEOM_SIZE
  } PmPotentialGeomType;

typedef enum PmPotentialType {
  PM_POTENTIAL_UNKNOWN,
  PM_POTENTIAL_AXIS,
  PM_POTENTIAL_BEND,
  PM_POTENTIAL_CONTACT,
  PM_POTENTIAL_ENTROPIC_SPRING,
  PM_POTENTIAL_MOLECULAR_MECHANICS,
  PM_POTENTIAL_SPRING,
  PM_POTENTIAL_TORSION
  } PmPotentialType;

typedef enum PmEntropicSpringFuncType {
  PM_ENTROPIC_SPRING_FUNC_UNKNOWN,
  PM_ENTROPIC_SPRING_FUNC_GRAPH,
  PM_ENTROPIC_SPRING_FUNC_WLC  
  } PmEntropicSpringFuncType;

class PmPotentialParameters {
  public:
     PmPotentialParameters() { use_sidechains = false; 
                               plane_set = false;
                               scales.set(1,1,1); 
                               hydrogen_bonds = false; }
     void setHydrogenBonds(const bool flag) { hydrogen_bonds = flag; }
     bool areHydrogenBonds() { return hydrogen_bonds; }
     void setUseSideChains(const bool flag) { use_sidechains = flag; }
     bool useSideChains() { return use_sidechains; }
     bool getPlaneData(PmVector3& point, PmVector3& normal) { 
                       point = this->point; normal = this->normal;
                       return plane_set; }
     void setPlaneData(PmVector3& point, PmVector3& normal) { 
                       this->point = point; this->normal = normal;
                       plane_set = true; }
     void setPlaneScales(PmVector3& scales) { 
                         this->scales = scales; } 
     void getPlaneScales(PmVector3& scales) { 
                         scales = this->scales; } 
  private:
     bool use_sidechains;
     bool hydrogen_bonds;
     bool plane_set;
     PmVector3 point, normal, scales;
  };

class PmPotentialMolMechTerms {
  public:
    PmPotentialMolMechTerms() { 
       lj_attractive = false;  lj_repulsive = false;  coulomb = false; 
       }
    bool setTerm(const string term); 
    bool lj_repulsive; 
    bool lj_attractive; 
    bool coulomb;
  };

class PmPotentialMolMechParameters {
  public:
     PmPotentialMolMechParameters(); 
     bool setParameter(const string param);
     float lj_repulsive_power;
     float lj_attractive_power;
     float charge, charge_scale;
     bool surface_charge;
  };


class PmPotential;

//*============================================================*
//*==========         PmPotentialGeomProps           ==========*
//*============================================================*

class PM_EXPORT PmPotentialGeomProps {
   public:
      PmPotentialGeomProps(){ has_charge=false; charge=0.0; };
      bool has_charge; float charge;
   };

//*============================================================*
//*==========         PmPotentialGeom                ==========*
//*============================================================*

class PM_EXPORT PmPotentialGeom {
   public:
      PmPotentialGeom();
      static PmPotentialGeom* create(const string name, const PmPotentialGeomType type, 
                                     PmSimulationObj *sobj, const string rgn,
                                     PmPotentialParameters& params);
      virtual void updateDisplay()=0;
      virtual void display(const bool show)=0;
      virtual void setXform(PmXform& xform)=0;
      virtual void print()=0;
      void setCharge(const float charge);
      void getCharge(float& charge);
      void getName(string& name);
      bool hasName(const string& name);
      void getSimulationObject(PmSimulationObj **sobj);
      void getType(PmPotentialGeomType *type);
      static void convType(const string str, PmPotentialGeomType& type);
      void addPotential(PmPotential *pot);
      void getPotentials(vector<PmPotential*>& list);

      void getXform(PmXform& xform);
      void setColor(PmVector3& color);
      void setDisplayType(PmGeometryDisplayType type);
      void setDisplaySpheres(const bool flag);
      void setLineWidth(float width);
      void setShadingType(PmGeometryShadingType type);

   //protected:
      string name;
      string rgn_name;
      PmPotentialGeomType type;
      PmPotentialGeomProps properties; 
      PmSimulationObj *simulation_obj;
      vector<PmPotential*> potential_list;
      PmXform potential_xform;
      PmPotentialParameters parameters;
      PmGraphicsGeometry *graphics_geometry, *contact_geometry, *graphics_geometry1; 
      PmVector3 color;
      PmGeometryDisplayType display_type;
      PmGeometryShadingType shading_type;
      float line_width;
      bool display_spheres, visible;
   };


// PmPotentialSphere
// -----------------

class PM_EXPORT PmPotentialSphere : public PmPotentialGeom {
   public:
      PmPotentialSphere();
      void setSphereRadius(const float radius);
      float radius;
      PmPcaResults pca;
      PmVector3 center, curr_center;
      bool contact;
      PmVector3 force;
      PmVector3 vel;
      float mass;
      void updateDisplay();
      void display(const bool show);
      void setXform(PmXform& xform);
      void print();
   };

// PmPotentialCylinder
// -------------------

class PM_EXPORT PmPotentialCylinder : public PmPotentialGeom {
   public:
      PmPotentialCylinder();
      float rads[3], max_rad, length;
      PmPcaResults pca;
      PmVector3 u, v, w, curr_u, curr_v, curr_w;
      PmVector3 center, curr_center;
      PmVector3 force;
      bool contact;
      PmVector3 contact_pt;
      void updateDisplay();
      void display(const bool show);
      void getGeometry(PmVector3& pos, PmVector3& u, PmVector3& v, PmVector3& w);
      void setXform(PmXform& xform);
      void print();
   };


// PmPotentialEllipsoid
// --------------------

class PM_EXPORT PmPotentialEllipsoid : public PmPotentialGeom {
   public:
      PmPotentialEllipsoid();
      void setEllipsoidScale (const float scale);
      float rads[3], max_rad, scale;
      PmPcaResults pca;
      PmVector3 u, v, w, curr_u, curr_v, curr_w;
      PmVector3 center, curr_center;
      PmVector3 force;
      bool contact;
      PmVector3 contact_pt;
      void updateDisplay();
      void display(const bool show);
      void getGeometry(PmVector3& pos, PmVector3& u, PmVector3& v, PmVector3& w);
      void setXform(PmXform& xform);
      void print();
   };

// PmPotentialLines 
// ----------------

class PM_EXPORT PmPotentialLines: public PmPotentialGeom {
   public:
      PmPotentialLines();
      void setAxes(const vector<int>& ids);
      void setAxes(const int id, const int aid);
      int number;
      PmVector3 *coordinates, *curr_coordinates;
      PmPcaResults pca;
      PmVector3 origin, curr_origin;
      float scale, scales[3];
      int num_axes;
      PmVector3 axes[3], pca_axes[3];
      int axes_ids[3];
      void updateDisplay();
      void display(const bool show);
      void setXform(PmXform& xform);
      void print();
   };


// PmPotentialPoints
// -----------------

class PM_EXPORT PmPotentialPoints : public PmPotentialGeom {
   public:
      PmPotentialPoints() {}
      void setSphereRadius(const float radius);
      int number;
      PmVector3 *coordinates, *curr_coordinates;
      float *radii, radius, rads[3];
      PmPcaResults pca;
      PmVector3 center, curr_center;
      bool contact;
      bool use_sidechains;
      void updateDisplay();
      void display(const bool show);
      void setXform(PmXform& xform);
      void getCurrentAxes(PmVector3& pos, PmVector3& u, PmVector3& v, PmVector3& w);
      void getPcaScales(float r1[3]);
      void print();
   };

// PmPotentialPlane 
// ----------------

class PM_EXPORT PmPotentialPlane: public PmPotentialGeom {
   public:
      PmPotentialPlane();
      PmVector3 *coordinates, *curr_coordinates;
      PmVector3 center, curr_center, normal;
      float distance;
      bool contact;
      float scale_x, scale_y;
      void setNormal(const PmVector3& vec);
      void setPoint(const PmVector3& pt);
      void setScale(float sx, float sy);
      void updateDisplay();
      void display(const bool show);
      void setXform(PmXform& xform);
      //void getCurrentAxes(PmVector3& pos, PmVector3& u, PmVector3& v, PmVector3& w);
      //void getPcaScales(float r1[3]);
      void print();
   };

//*============================================================*
//*==========             PmPotential                ==========*
//*============================================================*

class PM_EXPORT PmPotential {
   public:
      PmPotential();
      PmPotential(const string name);
      void init();
      static PmPotential* create(const string name, const PmPotentialType type, 
                                 PmPotentialGeom *geom1,  PmPotentialGeom *geom2);
      //void setRegion(string str) { region = str; } 
      virtual void updateDisplay()=0;
      void getEnergy(const float time, const bool active, float& val);
      void getEnergy(float& val);
      void setReportEnergy(const bool flag); 
      void getReportEnergy(bool& flag); 
      virtual void compForces(const bool upd, const float time, const bool active, 
                              const bool verbose)=0;
      void displayGeometry(const bool show);
      void getGeometry(PmPotentialGeom **geom) { *geom = geometry1; }
      void printGeometry();
      void getName(string& name) { name = this->name;  }
      bool hasName(const string& name) { return name == this->name;  }
      virtual int getNumInteractions()=0;
      static void procQuery(PmQuery& query);
      void getStrain(const float time, const bool active, float& val);
      void setMaxStrain(const float val);
      void getMaxStrain(float& val);
      void addTime(const float val);
      bool getTime(const int n, float& t);
      PmPotentialType getType() { return type; }
      static void convType(const string str, PmPotentialType& type);
      static void convType(const PmPotentialType type, string& str);
      void setXform(PmXform& xform);
      void getXform(PmXform& xform);
      void setActivated(const bool flag);

      void setColor(PmVector3& color);
      void setLineWidth(float width);
      virtual void queryState(PmQuery& query)=0;

   protected:
      string name;
      PmPotentialType type;
      PmPotentialGeomType geometry_type1, geometry_type2;
      PmPotentialGeom *geometry1, *geometry2;
      PmSimulationObj *sobj1, *sobj2;
      string region;
      bool verbose;
      bool report_energy;
      float current_energy;
      float current_strain, max_strain;
      bool initialized;
      bool activated; 
      void initSpatialDecomp(PmPotentialPoints *points, 
                             PmSpatialDecomp **p_spatial_decomp, float width);
      PmVector3 color;
      float line_width;
   };


// PmContactPotential
// ------------------

class PM_EXPORT PmContactPotential : public PmPotential {
  public:
     PmContactPotential(){};
     PmContactPotential(const string name);
     void checkContact(bool& flag);
     void setContactStrength(float strength);
     void updateDisplay();
     void compForces(const bool upd, const float time, const bool active, 
                     const bool verbose);
     int getNumInteractions(){ return 1;};

     void setImpulse(const bool val);
     void setRestitutionCoef(const float val);

     void compImpulse (PmContactPotential *potB, PmVector3& cpt, 
                       PmVector3& normal, float *impulse);
     void queryState(PmQuery& query){};

     // force computation functions //
     void no_forces(PmContactPotential *pot, const bool flag){};
     void cylinder_cylinder_forces(PmContactPotential *pot, const bool flag);
     void cylinder_ellipsoid_forces(PmContactPotential *pot, const bool flag);
     void cylinder_plane_forces(PmContactPotential *pot, const bool flag);
     void cylinder_points_forces(PmContactPotential *pot, const bool flag);
     void cylinder_sphere_forces(PmContactPotential *pot, const bool flag);

     void ellipsoid_cylinder_forces(PmContactPotential *pot, const bool flag);
     void ellipsoid_ellipsoid_forces(PmContactPotential *pot, const bool flag);
     void ellipsoid_plane_forces(PmContactPotential *pot, const bool flag);
     void ellipsoid_points_forces(PmContactPotential *pot, const bool flag);
     void ellipsoid_sphere_forces(PmContactPotential *pot, const bool flag);

     void plane_ellipsoid_forces(PmContactPotential *pot, const bool flag);
     void plane_points_forces(PmContactPotential *pot, const bool flag);
     void plane_sphere_forces(PmContactPotential *pot, const bool flag);

     void points_cylinder_forces(PmContactPotential *pot, const bool flag);
     void points_ellipsoid_forces(PmContactPotential *pot, const bool flag);
     void points_plane_forces(PmContactPotential *pot, const bool flag);
     void points_points_forces(PmContactPotential *pot, const bool flag);
     void points_sphere_forces(PmContactPotential *pot, const bool flag);

     void sphere_cylinder_forces(PmContactPotential *pot, const bool flag);
     void sphere_ellipsoid_forces(PmContactPotential *pot, const bool flag);
     void sphere_points_forces(PmContactPotential *pot, const bool flag);
     void sphere_plane_forces(PmContactPotential *pot, const bool flag);
     void sphere_sphere_forces(PmContactPotential *pot, const bool flag);

     float strength;
     float current_energy;
     bool contact;
     bool use_impulse;
     float restitution;
     float maximum_force;
     PmSpatialDecomp *spatial_decomp;
     void initPotential();
  };

// PmEntropicSpringPotential
// -------------------------

class PM_EXPORT PmEntropicSpringPotential : public PmPotential {
  public:
     PmEntropicSpringPotential(){};
     PmEntropicSpringPotential(const string name);
     void updateDisplay();
     void setDisplayConn(const bool flag);
     void compForces(const bool upd, const float time, const bool active, 
                     const bool verbose);
     void setForceConst(const float val);
     void setForceGraphData(const vector<float>& data);
     int getNumInteractions() { return 1; };

     static void convFunctionType(const string str, PmEntropicSpringFuncType& type);
     void setForceFunctionType(PmEntropicSpringFuncType type);
     void setWlcParams(const float lp, const float lc);
     void queryState(PmQuery& query){};

   private:
     float force_const, cutoff;
     bool show_conn;
     PmGraphicsGeometry *graphics_conn_geometry; 
     PmVector3 *connectivity_verts;
     PmVector3 conn_color;
     float conn_line_width;
     PmEntropicSpringFuncType force_function_type;
     void *force_function;
     void displayConnectivity();
   };


// PmSpringPotential
// -----------------

class PM_EXPORT PmSpringPotential : public PmPotential {
  public:
     PmSpringPotential(){};
     PmSpringPotential(const string name);
     void setCutoff(const float val);
     void updateDisplay();
     void setDisplayConn(const bool flag);
     void getDistance(const int i, float& dist);
     void setMaxDistance(const float val);
     void setMapEnergy(const bool flag);
     void setMapScale(const float min, const float max);
     void compForces(const bool upd, const float time, const bool active, 
                     const bool verbose);
     void setForceConst(const float val);
     void getForceConst(float& val);
     void setPairedConnectivity(const bool flag);
     int getNumInteractions();
     void getNumSprings(int& num);
     void setMapStrain(const bool flag);
     void getStrain(const int i, float& strain);
     void queryState(PmQuery& query);
     void setLJSpring(const bool val);
     void getLJSpring(bool& val);

   private:
     float force_const, cutoff;
     vector<float> force_consts; 
     int number_springs;
     bool no_connections_warning;
     bool *spring_active; 
     bool ljspring; 
     bool connectivity_computed;
     int *connectivity;
     float *distances, *strains, *energy, *current_distances;
     float max_distance;
     bool show_conn;
     bool paired_conn;
     PmGraphicsGeometry *graphics_conn_geometry; 
     PmVector3 *connectivity_verts;
     PmVector3 *connectivity_colors;
     void compConnectivity();
     void displayConnectivity();
     bool map_strain, map_energy, map_scale_set;
     float map_scale_min, map_scale_max;
     void getLJSpringConst (const string atom_name1, const string atom_name2, 
                            const float dist, float& k, float& r0);
   };

// PmMolecularMechanicsPotential
// -----------------------------

class PM_EXPORT PmMolecularMechanicsPotential: public PmPotential {
  public:
     PmMolecularMechanicsPotential(){};
     PmMolecularMechanicsPotential(const string name) ;
     void setCutoff(const float val);
     void setTerms(const PmPotentialMolMechTerms& terms);
     void updateDisplay();
     void compForces(const bool upd, const float time, const bool active,
                     const bool verbose);

     void compCoulombForces(const bool upd, const float time, const bool active,
                            const bool verbose);

     void setMaxForce(const float val);
     void getParameters(PmPotentialMolMechParameters& params);
     void setParameters(const PmPotentialMolMechParameters& params);
     void setStrength(const float val);
     int getNumInteractions(){ return 0;};
     void queryState(PmQuery& query);

   private:
     float cutoff, maximum_force;
     float strength;
     int num_contact;
     PmPotentialMolMechTerms terms;
     PmPotentialMolMechParameters parameters;
     PmGraphicsGeometry *graphics_geometry;
     PmSpatialDecomp *spatial_decomp;
     vector<float> charge_data1, charge_data2;
     void initPotential();
     void initChargeData();
   };


// PmBendPotential
// ----------------

class PM_EXPORT PmBendPotential : public PmPotential {
  public:
     PmBendPotential(){};
     PmBendPotential(const string name);
     void updateDisplay();
     void setDisplayConn(const bool flag);
     void compForces(const bool upd, const float time, const bool active, 
                     const bool verbose);
     void setForceConst(const float val);
     int getNumInteractions();
     void setPoint(const PmVector3& point);
     void setOrigin(const PmVector3& origin);
     void queryState(PmQuery& query){};

   private:
     PmVector3 origin, point1, point2;
     bool has_origin; 
     bool init; 
     float init_angle;
     
     float force_const;
     bool show_conn;
     PmGraphicsGeometry *graphics_conn_geometry;
     PmVector3 *connectivity_verts;
     void displayConnectivity();
     void compConnectivity(PmBendPotential *pot);
     PmVector3 conn_color;
     float conn_line_width;
     void initConfig();
   };


// PmAxisPotential
// ---------------

class PM_EXPORT PmAxisPotential : public PmPotential {
  public:
     PmAxisPotential(){};
     PmAxisPotential(const string name);
     void updateDisplay();
     void compForces(const bool upd, const float time, const bool active,
                     const bool verbose);
     void setAxesStrength(const vector<float>& vals);
     void setAxesStrength(const int id, const float val);
     int getNumInteractions(){ return 1;};
     void setPointStrength(const float val);
     void setRestStrength(const float val);
     void setRestCutoff(const float val);
     void setDisplayConn(const bool flag);
     void queryState(PmQuery& query){};

   private:
     float axes_strength[3];
     float restraint_strength, restraint_cutoff;
     vector<PmVector3> points; 
     float point_strength; 
     bool init;

     bool show;
     PmGraphicsGeometry *graphics_geometry;
     PmVector3 *conn_verts;
     PmVector3 color;
     float line_width;
     void displayConnectivity();
   };


// PmTorsionPotential
// ----------------

class PM_EXPORT PmTorsionPotential : public PmPotential {
  public:
     PmTorsionPotential(){};
     PmTorsionPotential(const string name);
     void updateDisplay();
     void setDisplayConn(const bool flag);
     void compForces(const bool upd, const float time, const bool active,
                     const bool verbose);
     void setForceConst(const float val);
     int getNumInteractions();
     void setPoints(const PmVector3& point1, const PmVector3& point2);
     void setOrigin(const PmVector3& origin);
     void queryState(PmQuery& query){};

   private:
     PmVector3 origin, point1, point2;
     bool has_origin;
     bool init;
     float init_angle;

     float force_const;
     bool show_conn;
     PmGraphicsGeometry *graphics_conn_geometry;
     PmVector3 *connectivity_verts;
     void displayConnectivity();
     void compConnectivity();
     PmVector3 conn_color;
     float conn_line_width;
     void initConfig();
   };
}

#endif



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
//* solid:                 s o l i d                           *
//*============================================================*

#ifndef _SOLID_PM_H_ 
#define _SOLID_PM_H_

#include "pm/pm.h"
#include "pobj.h"
#include "graphics.h"

namespace ProteinMechanica {

typedef enum {
  PM_SOLID_UNKNOWN,
  PM_SOLID_CYLINDER,
  PM_SOLID_ELLIPSOID,
  PM_SOLID_SPHERE
  } PmSolidType;

class PmCylinder;
class PmSphere;
class PmEllipsoid;
class PmGraphicsGeometry;

// PmSolidRegion
// -----------

class PM_EXPORT PmSolidRegionParameters {
  public:
     PmSolidRegionParameters() {
       use_distance = false;
       distance = 0.0;
       tolerance = 0.1;
       }

     bool use_distance;
     PmVector3 point;
     float distance;
     float tolerance;
  };

class PM_EXPORT PmSolidRegion : public PmRegion {
  public:
     PmSolidRegion() { };
     void getCoords() { };
     void getData() { };
     void getIndices() { };
     void getRadii() { };
     string descriptor;
  };


// PmSolidAttributes
// --------------------

class PmSolidParameters {
  public:
    bool length_set; float length;
    bool radius_set; float radius;

    PmSolidParameters() { init(); }

    PmSolidParameters& operator=(const PmSolidParameters& src) {
      if (src.length_set) length = src.length;
      if (src.radius_set) radius = src.radius;
      return *this;
      }

    void setLength(const float length) {
      length_set = true;
      this->length = length;
      }

    void init() {
      length_set = false;
      length = 1.0;
      radius_set = false;
      radius = 1.0;
      }
  };


// PmSolid 
// -------

 class PM_EXPORT PmSolid : public PmPhysicalObj {
   public:
     PmSolid(){};
     PmSolid(const string name);
     void getCoordinates(vector<PmVector3>& coords); 
     void setDensity(const float densit) { this->density = density; }
     void getIndices(const string desc, vector<int>& indices);
     void getName(string& name) { name = this->name; }
     void getRadii(vector<float>& rads);
     void defineRegion(const string name, const vector<int>& ids);
     void defineRegion(const string name, const string desc,
                       const PmSolidRegionParameters& params);
     void displayRegion(const string name, PmVector3 color, bool use_spheres);
     void getRegion(const string name, PmRegion **rgn);
     void getXform(PmXform& xform);
     void setXform(PmXform& xform);
     void update(int num_verts, PmVector3 *verts);

     void setColor (PmVector3& color) { this->color = color; }
     void setDisplayType (PmGeometryDisplayType type) { this->display_type = type; }
     void setShadingType (PmGeometryShadingType type) { this->shading_type = type; }
     void addGeometry (PmGraphicsGeometry *geom);
     void buildGeomName (string type, string& grname);
     void getGeometry (string name, PmGraphicsGeometry **geom);

     static PmSolid* create(const string name, const PmSolidType type);
     static void convSolidType(const string str, PmSolidType& type);
     static void getCylinder(const string str, PmCylinder **cyl);
     static void getSphere(const string str, PmSphere **sphere);
     static void getEllipsoid(const string str, PmEllipsoid **ellipsoid);

     virtual void display(bool show) = 0;
     virtual void getDimensions(vector<float>& dims) = 0;
     virtual void getExtent(PmExtent& extent) = 0;
     virtual void getMassProps(PmMassProperties& props) = 0;
     virtual void defineRegion(const string name, const bool use_center) = 0;

   protected:
     string name;
     PmSolidType type;
     PmSolidParameters parameters;
     PmMassProperties mass_properties;
     float density;
     PmXform xform;
     vector<float> dimensions;

     int num_vertices; 
     PmVector3 *vertices, *vertex_normals, *face_normals;
     int num_triangles;
     PmConn *connectivity; 

     PmGraphicsGeometry *graphics_geometry;
     vector<PmGraphicsGeometry*> graphics_geometries;
     PmVector3 color;
     PmGeometryDisplayType display_type;
     PmGeometryShadingType shading_type;
   };
   

// PmCylinder 
// ----------
// cylinder 

 class PM_EXPORT PmCylinder : public PmSolid {
   public:
      PmCylinder(const string name);
      //void create();
      void display(bool show);
      void getDimensions(vector<float>& dims);
      void getMassProps(PmMassProperties& props);
      void setRadius(const float r) { this->radius = r; }
      void setLength(const float len) { this->length = len; }
      void setAxis(const PmVector3& axis);
      void setCenter(const PmVector3& pos) { this->center = pos; }
      void getExtent(PmExtent& extent);
      void getParameters(PmVector3& center, float& radius, float&length, PmVector3& axis) 
        { center = this->center; radius = this->radius; length = this->length; 
          axis = this->axis; }
     void defineRegion(const string name, const bool use_center);

    private:
      float radius, length;
      PmVector3 center;
      PmVector3 axis;
      void genGeometry();
   };


// PmSphere
// --------
// sphere  

 class PM_EXPORT PmSphere: public PmSolid {
   public:
      PmSphere(const string name);
      //void create();
      void display(bool show);
      void getMassProps(PmMassProperties& props);
      void setRadius(const float r) { this->radius = r; }
      void setCenter(const PmVector3& pos) { this->center = pos; }
      void getDimensions(vector<float>& dims);
      void getExtent(PmExtent& extent);
      void getParameters(PmVector3& center, float& radius)
        { center = this->center; radius = this->radius; }
     void defineRegion(const string name, const bool use_center);

    private:
      float radius;
      PmVector3 center;
      void genGeometry();
   };

// PmEllipsoid
// -----------
// ellipsoid

 class PM_EXPORT PmEllipsoid: public PmSolid {
   public:
      PmEllipsoid (string name);
      //void create();
      void display(bool show);
      void getMassProps(PmMassProperties& props);
      void setAxes(const vector<PmVector3>& axes);
      void setRadii(const vector<float>& rads);
      void setCenter(const PmVector3& pos) { this->center = pos; }
      void getDimensions(vector<float>& dims);
      void getExtent(PmExtent& extent);
      void getParameters(PmVector3& center, vector<float>& radii, 
                         vector<PmVector3>& axes);
     void defineRegion(const string name, const bool use_center);

    private:
      vector<float> radii;
      PmVector3 center;
      vector<PmVector3> axes;
      void genGeometry();
   };

}

#endif


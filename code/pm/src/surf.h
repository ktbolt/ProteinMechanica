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
 * surf:             s u r f a c e                            *
 *============================================================*/

#ifndef _SURFACE_PM_H_
#define _SURFACE_PM_H_

#include "pm/pm.h"
#include "pobj.h"
#include "graphics.h"

namespace ProteinMechanica {

class PM_EXPORT PmSurfRegionParameters {
  public:
     PmSurfRegionParameters() {
       use_distance = false;
       use_face = false;
       distance = 0.0;
       tolerance = 0.1;
       radius_set = false;
       radius = 0.0;
       charge_min = -10.0;
       charge_max =  10.0;
       }

     bool use_distance;
     bool use_face, radius_set;
     PmVector3 point;
     float distance;
     float tolerance;
     float radius;
     string data_name;
     float charge_min, charge_max;
  };

// PmSurfRegion
// ------------

class PM_EXPORT PmSurfRegion : public PmRegion {
  public:
     PmSurfRegion() { };
     void getCoords() { };
     void getRadii() { };
     void getData() { };
     void getIndices() { };
     string descriptor;
  };

class PM_EXPORT PmSurfData {
  public:
    PmSurfData() { };
    string name;
    float *values;
    float vmin, vmax; 
  };


// PmSurface 
// ---------

 class PM_EXPORT PmSurface : public PmPhysicalObj {

    public:
       PmSurface(const string name, int num_verts, PmVector3 *coords,
                 int num_cells, int *conn, int nodes_per_cell);
       PmSurface(const string name, const string type); 
       PmSurface(const string name, vector<PmVector3>& coords, float thickness);
       ~PmSurface();

       void getCellCenters(vector<PmVector3>& centers);
       void setMapCharge(bool flag);
       void mapCharge(vector<PmVector3>& colors);
       void setChargeMap(float cmin, float cmax);
       void getConnectivity(int& num, int& nodes_per_cell, vector<int>& conn);
       void getCoordinates(vector<PmVector3>& coords);
       void getDimensions(vector<float>& dims);
       void getIndices(const string desc, vector<int>& indices);
       void getMassProps (PmMassProperties& props);
       void getRadii(vector<float>& rads);

       void addVertexData(string name, float *vals); 
       void getVertexData(string name, vector<float>& vals);

       void display(const bool show);
       void setColor (PmVector3& color) { this->color = color; }
       void setDisplayType (PmGeometryDisplayType type) { this->display_type = type; }
       void setMarkerSize(const float size) { this->marker_size = size; }
       void setShadingType (PmGeometryShadingType type) { this->shading_type = type; }
       void getName (string& name) { name = this->name;  }
       void getExtent (PmExtent& extent ) { extent = this->extent; }
       void setExtent (PmExtent& extent ) { this->extent = extent; }

       void setLighting(bool flag) { this->lighting = flag; }

       void setRegionRadius(float radius);
       void setRegionUseFace(bool flag);
       void defineRegion(const string name, const vector<int>& ids);

       void defineRegion(const string name, const string desc, 
                         const PmSurfRegionParameters& params);

       void displayRegion(const string name, PmVector3 color, bool use_spheres);
       void getRegion(const string name, PmRegion **rgn);

       void getXform (PmXform& xform);
       void setXform (PmXform& xform);
       void xformCoordinates(PmXform& xform);

       void static genPlates(string name, string type, vector<PmSurface*>& surfs);

    private:

      string name;
      int num_coordinates;
      PmVector3 *coordinates;
      int num_cells;
      int *connectivity;
      int nodes_per_cell;
      PmVector3 *cell_normals;
      PmMassProperties mass_properties;
      PmExtent extent;
      PmXform xform;
      vector<float> dimensions;
      bool plate;
      float plate_thickness;

      float rgn_radius;
      bool rgn_use_face;

      vector<PmSurfData*> vertex_data; 

      void genIcosahedron(string name, bool truncated);

      void getCellNormals(PmVector3 **normals);

      PmVector3 color;
      PmGeometryDisplayType display_type;
      PmGeometryShadingType shading_type;
      float marker_size;
      bool map_charge;
      bool lighting;
      float charge_min, charge_max;

      vector<PmGraphicsGeometry*> graphics_geometry;
      void buildGeomName (string type, string& gr_name);
      void addGeometry (PmGraphicsGeometry *geom);
      void getGeometry (string name, PmGraphicsGeometry **geom);

   };

}

#endif


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
//* surf:             s u r f a c e                            *
//*============================================================*

#include "surf.h"
#include "graphics.h"
#include "pm/mth.h"
#include "gc.h"

namespace ProteinMechanica {

#define ndbg_PmSurface

#define X 0.525731112119133606 
#define Z 0.850650808352039932

static float vdata[12][3] = { {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
                              {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
                              {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} };

static int tindices[20][3] = { 
               {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
               {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
               {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
               {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };

static int tedges[12][6] = { 
                             {0, 11, 6, 1, 4, 9}, 
                             {1, 10, 8, 4, 0, 6}, 
                             {2, 11, 9, 5, 3, 7}, 
                             {3, 7, 2, 5,  8, 10}, 
                             {4, 8, 5, 9,  0, 1}, 
                             {5, 9, 4, 8,  3, 2}, 
                             {6, 1, 0,11,  7, 10}, 
                             {7, 2, 3, 10, 6, 11},  
                             {8, 5, 4, 1, 10, 3}, 
                             {9, 4, 5, 2, 11, 0}, 
                             {10, 6, 7, 3, 8, 1}, 
                             {11, 7, 6, 0, 9, 2}};


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSurface::PmSurface(const string name, int num_coords, PmVector3 *coords,
                     int num_cells, int *conn, int nodes_per_cell) 
  {
  #ifdef dbg_PmSurface
  fprintf (stderr, "\n>>>>>> PmSurface:: ctor\n");
  fprintf (stderr, "   >>> num_cells [%d] \n", num_cells);
  fprintf (stderr, "   >>> nodes_per_cell [%d] \n", nodes_per_cell);
  #endif
  this->name = name;
  this->obj_type = PM_OBJECT_SURFACE;
  this->num_coordinates = num_coords;
  this->coordinates = coords;
  this->num_cells = num_cells;
  this->connectivity = conn;
  this->nodes_per_cell = nodes_per_cell;
  this->cell_normals = NULL;
  this->plate = false;
  this->plate_thickness = 0.0;
  this->rgn_radius = 0.02;
  this->rgn_use_face = false;

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  marker_size = 1.0;
  map_charge = false;
  charge_min = -1.8; 
  charge_max =  1.8;
  lighting = false;
  }

// generate a surface from a given type //

PmSurface::PmSurface(const string name, const string type)
  {
  if (type == "icosahedron") {
    genIcosahedron(name, false);
    }
  else if (type == "truncated-icosahedron") {
    genIcosahedron(name, true);
    }
  }

// generate a plate //

PmSurface::PmSurface(const string name, vector<PmVector3>& coords, float thickness)
  {
  PmVector3 *verts;
  int num_coords, *conn, n;
  num_coords = coords.size();

  verts = new PmVector3[num_coords];
  conn = new int[num_coords+1];

  /*
  fprintf (stderr, ">>>>>> PmSurface::ctor \n");
  fprintf (stderr, ">>> name = %s \n", name.c_str());
  fprintf (stderr, ">>> num coords = %d \n", num_coords);
  */

  for (int i = 0; i < num_coords; i++) {
    verts[i][0] = coords[i][0];
    verts[i][1] = coords[i][1];
    verts[i][2] = coords[i][2];
    //fprintf (stderr, "%f %f %f \n", verts[i][0], verts[i][1], verts[i][2]); 
    }

  conn[0] = num_coords; 

  for (int i = 0; i < num_coords; i++) {
    conn[i+1] = i; 
    }

  this->name = name;
  this->obj_type = PM_OBJECT_SURFACE;
  this->num_coordinates = num_coords;
  this->coordinates = verts;
  this->num_cells = 1;
  this->connectivity = conn;
  this->nodes_per_cell = num_coords;
  this->cell_normals = NULL;
  this->plate = true;
  this->plate_thickness = thickness;

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  }

void 
PmSurface::genPlates(string name, string type, vector<PmSurface*>& surfs)
  {
  PmSurface *surf;
  string sname;

  fprintf (stderr, ">>>>>> PmSurface::genPlates \n");
  fprintf (stderr, ">>> type = %s \n", type.c_str());
  surfs.clear();

  if (type == "truncated-icosahedron") {
    PmVector3 v, v1, v2, v3, dv; 
    vector<PmVector3> verts; 
    int i1, i2, i3;
    stringstream dss;

    for (int i = 0; i < 20; i++) {
      verts.clear();
      i1 = tindices[i][2];
      i2 = tindices[i][1];
      i3 = tindices[i][0];

      for (int j = 0; j < 3; j++) {
        v1[j] = vdata[i1][j];
        v2[j] = vdata[i2][j];
        v3[j] = vdata[i3][j];
        }

      dv = (v2 - v1) / 3.0;
      v = v1 + dv; verts.push_back(v);
      v = v + dv;  verts.push_back(v);

      dv = (v3 - v2) / 3.0;
      v = v2 + dv; verts.push_back(v);
      v = v + dv;  verts.push_back(v);

      dv = (v1 - v3) / 3.0;
      v = v3 + dv; verts.push_back(v);
      v = v + dv;  verts.push_back(v);

      dss << name << i+1;
      sname = dss.str();
      dss.str(std::string());
      surf = new PmSurface(sname, verts, 1.0);
      surfs.push_back(surf);
      }

    for (int i = 0; i < 12; i++) {
      verts.clear();
      i1 = tedges[i][0];

      for (int j = 1; j < 6; j++) {
        i2 = tedges[i][j];

        for (int j = 0; j < 3; j++) {
          v1[j] = vdata[i1][j];
          v2[j] = vdata[i2][j];
          }

        v = v1 + (v2 - v1) / 3.0;
        verts.push_back(v);
        }

      dss << "p" << name << i+1;
      sname = dss.str();
      dss.str(std::string());
      surf = new PmSurface(sname, verts, 1.0);
      surfs.push_back(surf);
      }

    }
  }

PmSurface::~PmSurface() {
  }

//*============================================================*
//*==========              getConnectivity           ==========*
//*============================================================*
// get the connectivity for the surface.

void
PmSurface::getConnectivity(int& num, int& np_cell, vector<int>& conn)
  {
  num = num_cells;
  np_cell = nodes_per_cell;
  int n = num_cells*nodes_per_cell;

  for (int i = 0; i < n; i++) {
    conn.push_back(connectivity[i]);
    }
  }

//*============================================================*
//*==========              setMapCharge              ==========*
//*============================================================*
// set the map charge flag to color the surface using charge.

void
PmSurface::setMapCharge(bool flag)
  {
  map_charge = flag;
  }

void 
PmSurface::setChargeMap(float cmin, float cmax)
  {
  charge_min = cmin;
  charge_max = cmax;
  }

//*============================================================*
//*==========              mapCharge                 ==========*
//*============================================================*
// map the charge to colors for the surface.

void
PmSurface::mapCharge(vector<PmVector3>& colors)
  {
  vector<float> vals;
  float vmin, vmax;
  PmVector3 color;
  float qmin, qmax;
  qmin = this->charge_min;
  qmax = this->charge_max;

  this->getVertexData("charge", vals);

  if (vals.empty()) {
    return;
    }

  vmin = vmax = vals[0];

  for (unsigned int i = 1; i < vals.size(); i++) {
    if (vals[i] > vmax) vmax = vals[i];
    if (vals[i] < vmin) vmin = vals[i];
    }

  fprintf (stderr, ">>> PmSurface::mapCharge  \n");      
  fprintf (stderr, ">>> num vals = %d \n", vals.size()); 
  fprintf (stderr, ">>> vmin = %f  vmax = %f \n", vmin, vmax);

  for (unsigned int i = 0; i < vals.size(); i++) {
    if (vals[i] <= qmin) {
      color.set(1,0,0);
      }
    else if (vals[i] >= qmax) {
      color.set(0,0,1);
      }
    else { 
      color.set(1,1,1);
      }

    colors.push_back(color);
    }
  }

//*============================================================*
//*==========              getCoordinates            ==========*
//*============================================================*
// get the coordinates for the surface.

void 
PmSurface::getCoordinates(vector<PmVector3>& coords)
  {
  coords.clear();

  for (int i = 0; i < num_coordinates; i++) {
    coords.push_back(coordinates[i]);
    }
  }

//*============================================================*
//*==========              addVertexData             ==========*
//*============================================================*
// add vertex data.                            

void 
PmSurface::addVertexData(string name, float *vals)
  {
  PmSurfData *sdata;
  sdata = new PmSurfData;
  sdata->name = name;
  sdata->values = vals;
  vertex_data.push_back(sdata);
  }

//*============================================================*
//*==========              getVertexData             ==========*
//*============================================================*
// get vertex data.

void
PmSurface::getVertexData(string name, vector<float>& vals)
  {
  for (unsigned int i = 0; i < vertex_data.size(); i++) {
    if (name == vertex_data[i]->name) { 
      for (int j = 0; j < num_coordinates; j++) {
        vals.push_back(vertex_data[i]->values[j]);
        }

      return;
      }
    }
  }

//*============================================================*
//*==========              getDimensions             ==========*
//*============================================================*
// get surface dimensions.

void
PmSurface::getDimensions(vector<float>& dims)
  {
  //dims.push_back(length);
  //dims.push_back(radius);
  }

//*============================================================*
//*==========              getRadii                  ==========*
//*============================================================*
// get the radii for the surface.

void
PmSurface::getRadii(vector<float>& rads)
  {
  rads.clear();

  for (int i = 0; i < num_coordinates; i++) {
    rads.push_back(0.2);
    }
  }

//*============================================================*
//*==========              xformCoordinates          ==========*
//*============================================================*
// transform the coordinates for the surface.

void
PmSurface::xformCoordinates(PmXform& xform)
  {
  PmMatrix3x3 mat = xform.matrix;
  PmVector3 pos; 

  for (int i = 0; i < num_coordinates; i++) {
    pos = coordinates[i];
    coordinates[i] = mat*(pos-xform.center) + xform.translation + xform.center;
    }

  mass_properties.set = false;
  }

//*============================================================*
//*==========              getCellCenters            ==========*
//*============================================================*
// get the cell centers for the surface.

void
PmSurface::getCellCenters(vector<PmVector3>& centers)
  {
  centers.clear();
  int k;
  PmVector3 center;

  for (int i = 0; i < num_cells; i++) {
    center.set(0,0,0);

    for (int j = 0; j < nodes_per_cell; j++) { 
      k = connectivity[nodes_per_cell*i+j];
      center = center + coordinates[k];
      }

    center = center / nodes_per_cell; 
    centers.push_back(center);
    }
  }

//*============================================================*
//*==========              getCellNormals            ==========*
//*============================================================*
// get the cell normals for the surface.

void
PmSurface::getCellNormals(PmVector3 **normals)
  {

  PmVector3 pverts[40], v1, v2, norm;
  float nx, ny, nz, mag;
  int *conn = connectivity;

  if (cell_normals) {
    *normals = cell_normals;
    return;
    }

  cell_normals = new PmVector3[num_cells];

  for (int i = 0; i < num_cells; i++) {
    for (int j = 0; j < nodes_per_cell; j++) {
      int k = *conn;
      pverts[j] = coordinates[k];
      conn++;
      }

    v1 = pverts[1]   - pverts[0];
    v2 = pverts[nodes_per_cell-1] - pverts[0];
    nx = v1[1]*v2[2] - v1[2]*v2[1];
    ny = v1[2]*v2[0] - v1[0]*v2[2];
    nz = v1[0]*v2[1] - v1[1]*v2[0];
    mag = sqrt(nx*nx + ny*ny + nz*nz);

    if (mag != 0.0) {
      nx = nx / mag;
      ny = ny / mag;
      nz = nz / mag;
      }

    cell_normals[i][0] = nx;
    cell_normals[i][1] = ny;
    cell_normals[i][2] = nz;
    }

  *normals = cell_normals;
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get the mass properties for the surface.

void
PmSurface::getMassProps (PmMassProperties& props)
  {
  #ifdef dbg_PmSurface
  fprintf (stderr, ">>>>>> PmSurface::getMassProps\n");
  #endif

  if (mass_properties.set) {
    props = mass_properties;
    return;
    }

  bool closed;

  pm_MathPolyClosed(num_cells, connectivity, nodes_per_cell, num_coordinates, closed);

  if (!closed) {
    pm_ErrorWarnReport (PM, "surface mesh \"%s\" is not closed.", "*", name.c_str());
    }

  float density, mass, J[3][3];
  PmVector3 com, *normals;
  density = 1.0;

  // get the normals to each cell //

  getCellNormals(&normals);

  // compute mass properties from the polygon mesh //

  if (plate) {
    com.set(0,0,0); 

    for (int i = 0; i < num_coordinates; i++) {
      com = com + coordinates[i]; 
      }

    com = com / num_coordinates; 
    mass = 1.0;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        J[i][j] = 0.0;
        if (i == j) J[i][j] = 1.0;
        }
      }
    }
  else {
    pm_MathPolyMassProps (num_cells, connectivity, nodes_per_cell, num_coordinates,
                          coordinates, normals, density, mass, com, J);
    }

  #ifdef dbg_PmSurface
  fprintf (stderr, "   >>> mass = %g \n", mass);
  fprintf (stderr, "   >>> com = %g %g %g \n", com[0], com[1], com[2]);
  fprintf (stderr, "   >>> J   = %g %g %g \n", J[0][0], J[0][1], J[0][2]);
  fprintf (stderr, "             %g %g %g \n", J[1][0], J[1][1], J[1][2]);
  fprintf (stderr, "             %g %g %g \n", J[2][0], J[2][1], J[2][2]);
  #endif

  mass_properties.set = true;
  mass_properties.mass = mass;
  mass_properties.com = com;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      mass_properties.inertia(i,j) = J[i][j];
      }
    }

  props = mass_properties;
  }

//*============================================================*
//*==========              setRegionRadius           ==========*
//*============================================================*
// set the radius of the spheres used in surface interactions.

void 
PmSurface::setRegionRadius(float radius)
  {
  rgn_radius = radius;
  }

//*============================================================*
//*==========              setRegionUseFace          ==========*
//*============================================================*
// set the flag to use spheres on surface faces.

void 
PmSurface::setRegionUseFace(bool flag)
  {
  rgn_use_face = flag;
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*

void
PmSurface::defineRegion(const string name, const vector<int>& ids)
  {
  PmVector3 center;
  int id, n;
  PmSurfRegion *region = new PmSurfRegion;
  region->name = name;

  n = 0;
  center.set(0,0,0);

  // define the region on the surface cell faces //

  if (rgn_use_face) {
    PmVector3 pts[3], c;
    float r, area;

    for (int i = 0; i < ids.size(); i++) {
      id = ids[i]-1;

      if ((id < 0) || (id >= num_cells)) {
        continue;
        }

      pts[0] = coordinates[connectivity[3*id]];
      pts[1] = coordinates[connectivity[3*id+1]];
      pts[2] = coordinates[connectivity[3*id+2]];

      pm_GcTriInscribedCircle(pts, c, r, area);
      region->coords.push_back(c);
      region->coord_indexes.push_back(id);
      region->radii.push_back(r);
      center = center + c;
      n += 1;
      }
    }

  // define the region on the surface vertices //

  else {
    for (int i = 0; i < ids.size(); i++) {
      id = ids[i]-1;

      if ((id < 0) || (id >= num_coordinates)) {
        continue;
        }

      region->coords.push_back(coordinates[id]);
      region->coord_indexes.push_back(id);
      region->radii.push_back(rgn_radius);
      center = center + coordinates[id];
      n += 1;
      }
    }

  if (!n) {
    return;
    }

  n = region->coords.size();
  region->center = center / n;
  regions.push_back(region);
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*

void 
PmSurface::defineRegion(const string name, const string desc,
                        const PmSurfRegionParameters& params)
  {
  #define dbgPmSurface_defineRegion
  #ifdef dbgPmSurface_defineRegion
  fprintf (stderr, "\n>>>>> PmSurface::defineRegion \n"); 
  fprintf (stderr, ">>> desc = %s \n", desc.c_str()); 
  fprintf (stderr, ">>> num cells = %d \n", num_cells); 
  fprintf (stderr, ">>> num coordinates = %d \n", num_coordinates); 
  fprintf (stderr, ">>> data = %s \n", params.data_name.c_str()); 
  #endif

  int n, id;
  PmVector3 center;
  vector<int> indices; 
  float radius;

  if (params.radius_set) {
    radius = params.radius;
    }
  else {
    radius = rgn_radius;
    }

  // descriptor using indices //
 
  if (desc.size()) {
    getIndices(desc, indices);
    }

  // distance from a point //

  else if (params.use_distance) {
    PmVector3 v, pt, point;
    float dist, d;

    point = params.point;
    dist = params.distance;

    if (params.use_face) {

      for (int i = 0; i < num_cells; i++) {
        pt = (coordinates[connectivity[3*i]] + coordinates[connectivity[3*i+1]] + 
              coordinates[connectivity[3*i+2]]) / 3.0;
        v = point - pt;
        d = v.length();

        if (d <= dist) {
          indices.push_back(i+1);
          }
        }
      }
    else {
      for (int i = 0; i < num_coordinates; i++) {
        pt = coordinates[i]; 
        v = point - pt;
        d = v.length();

        if (d <= dist) {
          indices.push_back(i+1);
          }
        }
      }
    }

  // all vertices or faces //

  else {
    if (params.use_face) {
      for (int i = 0; i < num_cells; i++) {
        indices.push_back(i+1);
        }
      }
    else {
      for (int i = 0; i < num_coordinates; i++) {
        indices.push_back(i+1);
        }
      }
    }

  if (!indices.size()) {
    pm_ErrorWarnReport (PM, "no surface coordinates found for region = \"%s\".", "*",
                        name.c_str());
    return;
    }

  PmSurfRegion *region = new PmSurfRegion;
  region->name = name;

  n = 0;
  center.set(0,0,0);

  // center spheres on surface triangles //

  if (params.use_face) {
    PmVector3 pts[3], c;
    float r, area;

    for (int i = 0; i < indices.size(); i++) {
      id = indices[i]-1;

      if ((id < 0) || (id >= num_cells)) {
        continue;
        }

      pts[0] = coordinates[connectivity[3*id]];
      pts[1] = coordinates[connectivity[3*id+1]];
      pts[2] = coordinates[connectivity[3*id+2]];

      pm_GcTriInscribedCircle(pts, c, r, area);

      if (r >= radius) {
        region->coords.push_back(c);
        region->coord_indexes.push_back(id);
        region->radii.push_back(r);
        center = center + c;
        n += 1;
        }
      }
    }

  // center spheres on surface vertices //

  else {

    if (params.data_name == "charge") { 
      vector<float> vals;
      this->getVertexData(params.data_name, vals);

      if (vals.empty()) {
        pm_ErrorWarnReport (PM, "surface \"%s\" does not have charge data. ", "*", 
                            name.c_str());
        return;
        }

      for (int i = 0; i < indices.size(); i++) {
        id = indices[i]-1;

        if ((id < 0) || (id >= num_coordinates)) {
          continue;
          }

        if ((vals[id] >= params.charge_min) && (vals[id] <= params.charge_max)) { 
          region->coords.push_back(coordinates[id]);
          region->coord_indexes.push_back(id);
          region->radii.push_back(rgn_radius);
          region->data.push_back(vals[id]);
          center = center + coordinates[id]; 
          n += 1;
          }
        }

      fprintf (stderr, ">>> charge region has n = %d \n", n);
      }
    else {
      for (int i = 0; i < indices.size(); i++) {
        id = indices[i]-1;

        if ((id < 0) || (id >= num_coordinates)) {
          continue;
          }

        region->coords.push_back(coordinates[id]);
        region->coord_indexes.push_back(id);
        region->radii.push_back(rgn_radius);
        center = center + coordinates[id]; 
        n += 1;
        }
      }
    }

  if (!n) {
    return;
    }

  n = region->coords.size();
  region->center = center / n;
  regions.push_back(region);
  }

//*============================================================*
//*==========              getRegion                 ==========*
//*============================================================*
// get a region for a molecule.

void
PmSurface::getRegion(const string name, PmRegion **rgn)
  {
  *rgn = NULL;

  if (!regions.size()) {
    return;
    }

  for (unsigned int i = 0; i < regions.size(); i++) {
    if (regions[i]->name == name) {
      *rgn = regions[i];
      return;
      }
    }
  }

//*============================================================*
//*==========              getIndices                ==========*
//*============================================================*

void
PmSurface::getIndices(const string desc, vector<int>& indices) 
  {
  //fprintf (stderr, "\n>>>>> PmSurface::getIndices \n");
  int i, j, n, m;
  int num, begin, end, rindex, pnum, size;
  char *str, val[2][80], name[20];

  str = (char*)desc.c_str();
  n = 0;
  m = 0;
  pnum = 0;
  //fprintf (stderr, "   >>> parse %s \n", str);
  str++;
  str++;

  while (*str != '\0') {
    char c = *str;
    char nc = *(str+1);

    if (!isdigit(c) || (nc == '\0')) {
      if (nc == '\0') {
        val[m][n++] = c;
        }

      val[m][n] = '\0';

      if (c == '-') {
        n = 0;
        m++;
        }

      else if ((c == ',') || (c == ']') || (nc == '\0')) {
        begin = atoi (val[0]); 

        if (!m) {
          end = begin; 
          }
        else {
          end = atoi (val[1]); 
          }

        for (i = begin; i <= end; i++) {
          indices.push_back(i);
          }

        m = 0;
        n = 0;

        if (c == ']') {
          break;
          }
        }
      }
    else {
      val[m][n++] = c;
      }
   
    str++;
    }

  /*
  fprintf (stderr, ">>> indices =  ");
  for (int i = 0; i < indices.size(); i++) {
    fprintf (stderr, " %d ", indices[i]);
    }
  fprintf (stderr, "\n");
  */
  }

//*============================================================*
//*==========              getXform                  ==========*
//*============================================================*

void
PmSurface::getXform (PmXform& xform) {
  xform = this->xform;
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the transformation for a surface.    

void 
PmSurface::setXform (PmXform& xform)
  {
  //fprintf (stderr, ">>>>>> PmSurface::setXform \n");
  this->xform = xform;

  for (int i = 0; i < graphics_geometry.size(); i++) {
    PmGraphicsGeometry *geom = graphics_geometry[i];
    geom->setXform(xform);
    geom->display();
    }
  }

////////////////////////////////////////////////////////////////
//                    g r a p h i c s                        //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             addGeometry                ==========*
//*============================================================*
// add graphics geometry object.

void
PmSurface::addGeometry (PmGraphicsGeometry *geom) {
  graphics_geometry.push_back (geom);
  }

//*============================================================*
//*==========               getGeometry              ==========*
//*============================================================*
// get a graphics geometry object.

void
PmSurface::getGeometry (string name, PmGraphicsGeometry **geom)
  {
  *geom = NULL;

  for (int i = 0; i < graphics_geometry.size(); i++) {
    if (name == graphics_geometry[i]->name) {
      *geom = graphics_geometry[i];
      }
    }
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display the surface.

void
PmSurface::display(const bool show)
  {
  #define ndbg_PmSurface
  #ifdef dbg_PmSurface
  fprintf (stderr, ">>>>>> PmSurface::display \n");
  fprintf (stderr, "   >>> name [%s] \n", name.c_str());
  fprintf (stderr, "   >>> num_cells [%d] \n", num_cells);
  fprintf (stderr, "   >>> nodes_per_cell [%d] \n", nodes_per_cell);
  fprintf (stderr, "   >>> map charge     [%d] \n", map_charge);
  #endif

  if (!pmSystem.useGraphics()) return;
  string geom_name;
  buildGeomName ("", geom_name);
  #ifdef dbg_PmSurface
  fprintf (stderr, "   >>> geom name [%s] \n", geom_name.c_str());
  #endif

  PmGraphicsGeometry *geom;
  PmGraphicsSurface *sgeom;
  getGeometry (geom_name, &geom);

  if (!geom) {
    int size = nodes_per_cell * num_cells;
    int hgeom;

    if (this->plate) {
      size += 1;
      hgeom = 0;
      }
    else {
      hgeom = nodes_per_cell;
      }

    PmConn *pconn = new PmConn(size);

    for (int i = 0; i < size; i++) {
      (*pconn)[i] = connectivity[i];
      }

    sgeom = new PmGraphicsSurface (geom_name, num_coordinates, coordinates, num_cells,
                                   pconn, hgeom);

    if (map_charge) {
      vector<PmVector3> colors;
      this->mapCharge(colors);
      sgeom->setColors(colors);
      }

    addGeometry (sgeom);
    geom = sgeom;
    }

  PmGraphicsAttributes atts;
  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setColor(color);
  atts.setLighting(lighting);
  atts.setScale(marker_size);
  geom->setAttributes(atts);
  geom->display();
  }

//*============================================================*
//*==========              displayRegion             ==========*
//*============================================================*
// display a region for a molecule.

void 
PmSurface::displayRegion(const string name, PmVector3 color, bool use_spheres)
  {
  PmRegion *rgn;

  getRegion(name, &rgn);

  if (!rgn) {
    return;
    }

  string geom_name;
  PmGraphicsPoint *points;
  PmGraphicsSphere *spheres;
  PmGraphicsGeometry *geom;
  int num;
  PmVector3 *verts;
  float *rads;
  PmGraphicsAttributes atts;

  num = rgn->coords.size();

  if (!num) {
    return;
    }

  buildGeomName (rgn->name, geom_name);
  getGeometry(geom_name, &geom);

  if (!geom) {
    verts = new PmVector3[num];

    for (int i = 0; i < num; i++) {
      verts[i] = rgn->coords[i];
      }

    if (use_spheres) {
      rads = new float[num];

      for (int i = 0; i < num; i++) {
        rads[i] = rgn->radii[i];
        }

      spheres = new PmGraphicsSphere(geom_name, num, verts, rads);
      addGeometry (spheres);
      geom = spheres;
      }
    else {
      points = new PmGraphicsPoint(geom_name, num, verts);
      addGeometry (points);
      geom = points;
      }
    }

  atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
  atts.setMarker(false);
  atts.setScale(2.0);
  atts.setColor(color);

  geom->setAttributes(atts);
  geom->display();
  }

//*============================================================*
//*==========              genIcosahedron            ==========*
//*============================================================*
// generate a surface for an icosahedron.

void
PmSurface::genIcosahedron (const string name, bool truncated)
  {

  PmVector3 *verts;
  int *conn, n; 
  int num_verts, num_polys;

  if (truncated) {

    }
  else {
    num_verts = 12;
    num_polys = 20;
    verts = new PmVector3[num_verts];
    conn = new int[3*num_polys];

    for (int i = 0; i < num_verts; i++) {
      verts[i][0] = vdata[i][0];
      verts[i][1] = vdata[i][1];
      verts[i][2] = vdata[i][2];
      }

    n = 0;

    for (int i = 0; i < num_polys; i++) {
      conn[n++] = tindices[i][2];
      conn[n++] = tindices[i][1];
      conn[n++] = tindices[i][0];
      }

    int i1, i2, i3;

    for (int i = 0; i < num_polys; i++) {
      i1 = tindices[i][2];
      i2 = tindices[i][1];
      i3 = tindices[i][0];

      fprintf (stderr, " ---- %d ---- \n", i+1); 
      fprintf (stderr, " %f %f %f \n", vdata[i1][0], vdata[i1][1], vdata[i1][2]);
      fprintf (stderr, " %f %f %f \n", vdata[i2][0], vdata[i2][1], vdata[i2][2]);
      fprintf (stderr, " %f %f %f\n \n", vdata[i3][0], vdata[i3][1], vdata[i3][2]);
      }
    }

  this->name = name;
  this->obj_type = PM_OBJECT_SURFACE;
  this->num_coordinates = num_verts;
  this->coordinates = verts;
  this->num_cells = num_polys;
  this->connectivity = conn;
  this->nodes_per_cell = 3;
  this->cell_normals = NULL;

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  }

////////////////////////////////////////////////////////////////
//                    p r i v a t e                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              buildGeomName             ==========*
//*============================================================*

void
PmSurface::buildGeomName (string type, string& grname) {
  grname = "surface";
  grname += '[' + name + ']';

  if (type.size()) {
    grname += type;
    }
  }

}



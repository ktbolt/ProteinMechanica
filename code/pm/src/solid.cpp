
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
//* solid:           s o l i d                                 *
//*============================================================*

#include "solid.h"
#include "graphics.h"
#include "pm/mth.h"
#include "gc.h"

namespace ProteinMechanica {

///////////////////////////////////////////////////////////////
//                       s o l i d                          //
/////////////////////////////////////////////////////////////

PmSolid::PmSolid(const string name) 
  {
  this->name = name;
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a solid object of a particular type.

PmSolid*
PmSolid::create(const string name, const PmSolidType type)
  {
  switch (type) {
    case PM_SOLID_CYLINDER:
      {
      PmCylinder *cyl = new PmCylinder(name);
      cyl->type = type;
      cyl->obj_type = PM_OBJECT_SOLID;
      return cyl;
      }
    break;

    case PM_SOLID_SPHERE:
      {
      PmSphere *sphere = new PmSphere(name);
      sphere->type = type;
      sphere->obj_type = PM_OBJECT_SOLID;
      return sphere;
      }
    break;

    case PM_SOLID_ELLIPSOID:
      {
      PmEllipsoid *ellipsoid = new PmEllipsoid(name);
      ellipsoid->type = type;
      ellipsoid->obj_type = PM_OBJECT_SOLID;
      return ellipsoid ;
      }
    break;

    default:
      return NULL;
    }
  }

//*============================================================*
//*==========               convSolidType            ==========*
//*============================================================*
// convert a srting solid type to a symbol. 

void 
PmSolid::convSolidType(const string str, PmSolidType& type)
  {
  type = PM_SOLID_UNKNOWN;

  if (str == "cylinder") {
    type = PM_SOLID_CYLINDER;
    }
  else if (str == "ellipsoid") {
    type = PM_SOLID_ELLIPSOID;
    }
  else if (str == "sphere") {
    type = PM_SOLID_SPHERE;
    }
  }

//*============================================================*
//*==========               getCylinder              ==========*
//*============================================================*
// get a cylinder object.                          

void      
PmSolid::getCylinder(const string str, PmCylinder **cyl)
  {
  *cyl = NULL;
  PmSolid *solid;
  pmSystem.getSolid(str, &solid);

  if (solid) {
    *cyl = dynamic_cast<PmCylinder*>(solid);
    }
  }

//*============================================================*
//*==========               getEllipsoid             ==========*
//*============================================================*
// get a sphere object.

void
PmSolid::getEllipsoid(const string str, PmEllipsoid **ellipsoid)
  {
  *ellipsoid = NULL;
  PmSolid *solid;
  pmSystem.getSolid(str, &solid);

  if (solid) {
    *ellipsoid = dynamic_cast<PmEllipsoid*>(solid);
    }
  }

//*============================================================*
//*==========               getSphere                ==========*
//*============================================================*
// get a sphere object.

void
PmSolid::getSphere(const string str, PmSphere **sphere)
  {
  *sphere = NULL;
  PmSolid *solid;
  pmSystem.getSolid(str, &solid);

  if (solid) {
    *sphere = dynamic_cast<PmSphere*>(solid);
    }
  }

//*============================================================*
//*==========               getCoordinates           ==========*
//*============================================================*
// get the coordinates for a solid.

void
PmSolid::getCoordinates(vector<PmVector3>& coords)
  {
  coords.clear();

  for (int i = 0; i < num_vertices; i++) {
    coords.push_back(vertices[i]);
    }
  }

//*============================================================*
//*==========               getRadii                 ==========*
//*============================================================*
// get the radii for a solid.

void
PmSolid::getRadii(vector<float>& rads)
  {
  rads.clear();

  for (int i = 0; i < num_vertices; i++) {
    rads.push_back(0.2);
    }
  }

//*============================================================*
//*==========               update                   ==========*
//*============================================================*
// update geometry with new vertices.

void
PmSolid::update(int num_verts, PmVector3 *verts)
  {
  }

//*============================================================*
//*==========              getXform                  ==========*
//*============================================================*

void
PmSolid::getXform (PmXform& xform) {
  xform = this->xform;
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the transformation for a solid.  

void
PmSolid::setXform (PmXform& xform)
  {
  this->xform = xform;

  if (!graphics_geometry) {
    return;
    }

  graphics_geometry->setXform(xform);
  graphics_geometry->display();
  }


//*============================================================*
//*==========              getIndices                ==========*
//*============================================================*

void
PmSolid::getIndices(const string desc, vector<int>& indices) 
  {
  /*
  fprintf (stderr, "\n>>>>>> PmSolid::getIndices \n");
  fprintf (stderr, ">>> desc=%s \n", desc.c_str()); 
  */
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
  for (int i = 0; i < indices.size(); i++) {
    fprintf (stderr, " %d ", indices[i]);
    }
  */
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a solid.

void
PmSolid::defineRegion(const string name, const string desc,
                       const PmSolidRegionParameters& params)
  {
  fprintf (stderr, "\n>>>>>> PmSolid::defineRegion  \n");
  fprintf (stderr, ">>> region name=%s  \n", name.c_str());
  int n, id;
  PmVector3 center;
  vector<int> indices; 

  if (desc.size()) {
    getIndices(desc, indices);
    }

  else if (params.use_distance) {
    PmVector3 v, pt, point;
    float dist, d;

    point = params.point;
    dist = params.distance;

    for (int i = 0; i < num_vertices; i++) {
      pt = vertices[i]; 
      v = point - pt;
      d = v.length();

      if (d <= dist) {
        indices.push_back(i+1);
        }
      }

    if (!indices.size()) {
      pm_ErrorWarnReport (PM, "no solid coordinates found for region = \"%s\".", "*",
                          name.c_str());
      return;
      }
    }

  fprintf (stderr, ">>> num indices=%d  \n", indices.size());
 
  PmSolidRegion *region = new PmSolidRegion;
  region->name = name;

  n = 0;
  center.set(0,0,0);

  for (int i = 0; i < indices.size(); i++) {
    id = indices[i];

    if ((id < 0) || (id >= num_vertices)) {
      continue;
      }

    region->coords.push_back(vertices[id-1]);
    region->coord_indexes.push_back(id-1);
    region->radii.push_back(0.02);
    center = center + vertices[id-1]; 
    n += 1;
    }

  fprintf (stderr, ">>> region has n = %d \n", n);

  if (!n) {
    return;
    }

  n = region->coords.size();
  region->center = center / n;
  regions.push_back(region);
  }

void 
PmSolid::defineRegion(const string name, const vector<int>& ids)
  {
  PmVector3 center;
  int id, n;
  PmSolidRegion *region = new PmSolidRegion;
  region->name = name;

  n = 0;
  center.set(0,0,0);

  for (int i = 0; i < ids.size(); i++) {
    id = ids[i]-1;

    if ((id < 0) || (id >= num_vertices)) {
      continue;
      }

    region->coords.push_back(vertices[id]);
    region->coord_indexes.push_back(id);
    region->radii.push_back(0.02);
    center = center + vertices[id];
    n += 1;
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
// get a region for a solid.     

void
PmSolid::getRegion(const string name, PmRegion **rgn)
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
//*==========              displayRegion             ==========*
//*============================================================*
// display a region for a solid.

void
PmSolid::displayRegion(const string name, PmVector3 color, bool use_spheres)
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

/////////////////////////////////////////////////////////////////
//            s o l i d   c y l i n d e r                     //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmCylinder::PmCylinder(const string name)
  {
  this->name = name;
  length = 1.0;
  radius = 0.5;
  axis.set(0,0,1);
  center.set(0,0,0);
  density = 1.0;

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  graphics_geometry = NULL;

  num_vertices = 0;
  connectivity = NULL;
  num_triangles = 0;
  vertices = NULL;
  vertex_normals = NULL; 
  face_normals = NULL;
  }

//*============================================================*
//*==========              setAxis                   ==========*
//*============================================================*
// set cylinder axis.

void 
PmCylinder::setAxis(const PmVector3& vec) 
  { 
  PmVector3 v = vec;

  if (v.length() == 0.0) { 
    return;
    }

  this->axis = vec;
  this->axis.normalize(); 
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display cylinder.

void
PmCylinder::display(bool show)
  {
  PmGraphicsAttributes atts;

  if (!graphics_geometry) {
    genGeometry();
    graphics_geometry = new PmGraphicsPolygon(name, num_triangles, connectivity, 3,
                                              num_vertices, vertices); 
    }

  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setColor(color);
  graphics_geometry->setAttributes(atts);
  graphics_geometry->display(); 
  }

//*============================================================*
//*==========              getDimensions             ==========*
//*============================================================*
// get cylinder dimensions.

void
PmCylinder::getDimensions(vector<float>& dims)
  {
  dims.push_back(length); 
  dims.push_back(radius); 
  }

//*============================================================*
//*==========              getExtent                 ==========*
//*============================================================*
// get cylinder extent.

void 
PmCylinder::getExtent(PmExtent& extent)
  {
  PmVector3 v;
  v = center + 0.5*length*axis;
  extent.update(v);

  v = center - 0.5*length*axis;
  extent.update(v);
  }

//*============================================================*
//*==========              genGeometry               ==========*
//*============================================================*
// generate cylinder geometry.

void
PmCylinder::genGeometry()
  {

  if (num_vertices) {
    return;
    }

  int i, j, k;
  PmVector3 u, v, w, p, pos, norm, v1, v2;
  float h, tl, t, dt, ct, st, dtl;
  float tlen;
  int n, i1, i2, i3, j1, j2, j3, ndiv, offset;

  n = 16;
  h = length / 2.0;
  dt = 2.0*M_PI / n;
  tlen = 2.0*M_PI*radius / n;
  ndiv = (int)(length / (4.0*tlen) + 0.5);
  if (ndiv == 0) ndiv = 1;
  dtl = length / ndiv;

  //fprintf (stderr, ">>>>>> PmCylinder::genCylinderGeom \n");
  //fprintf (stderr, "   >>> ndiv[%d] \n", ndiv); 

  // get cylinder basis //

  u = axis;
  pm_MathBasisCompute (u, v, w);

  // first create vertices //

  num_vertices = n*(ndiv+1) + 2;
  vertices = new PmVector3[num_vertices];
  vertex_normals = new PmVector3[num_vertices];
  num_vertices = 0; 

  tl = 0.0;
  pos = center - h*u;
  /*
  fprintf (stderr, "   >>> radius[%f] \n", radius);
  fprintf (stderr, "   >>> u (%f %f %f)  \n", u[0], u[1], u[2]);
  fprintf (stderr, "   >>> v (%f %f %f)  \n", v[0], v[1], v[2]);
  fprintf (stderr, "   >>> w (%f %f %f)  \n", w[0], w[1], w[2]);
  */

  for (k = 0; k <= ndiv; k++) {
    t = 0.0;

    for (i = 0; i < n; i++) {
      ct = radius*cos(t);
      st = radius*sin(t);
      //fprintf (stderr, "   >>> vertices[%d] ", num_vertices);

      for (j = 0; j < 3; j++) {
        vertices[num_vertices][j]         = pos[j] + ct*v[j] + st*w[j];
        vertex_normals[num_vertices][j]   = ct*v[j] + st*w[j];
        //fprintf (stderr, " %f ", vertices[num_vertices][j]); 
        }
 
      //fprintf (stderr, "\n");
      num_vertices += 1;
      t += dt;
      }

    tl += dtl;
    pos = pos + dtl*u;
    }

  //fprintf (stderr, "   >>> num_vertices[%d] \n", num_vertices); 

  // end vertices //

  for (j = 0; j < 3; j++) {
    vertices[num_vertices][j]   = center[j] - h*u[j];
    vertices[num_vertices+1][j] = center[j] + h*u[j];
    vertex_normals[num_vertices][j]  = -axis[j];
    vertex_normals[num_vertices+1][j] = axis[j];
    }

  num_vertices += 2;

  // now generate connectivity //

  num_triangles = 2*n*ndiv + 2*n;

  connectivity = new PmConn(3*num_triangles);
  int *pconn = connectivity->vals;
  face_normals = new PmVector3[num_triangles];
  num_triangles = 0;
  offset = 0;

  for (k = 0; k < ndiv; k++) {
    for (i = 0; i < n; i++) {
      i1 = i; i2 = i+1; i3 = i+n;

      if (i2 == n) {
        i2 = 0;
        }

      i1 += offset;
      i2 += offset;
      i3 += offset;
      //fprintf (stderr, "   >>> i1 i2 i3 %d %d %d \n", i1, i2, i3); 

      pconn[3*num_triangles+0] = i1;
      pconn[3*num_triangles+1] = i2;
      pconn[3*num_triangles+2] = i3;
      face_normals[num_triangles] = (1/3.0)*(vertex_normals[i1] + vertex_normals[i2] + 
                                       vertex_normals[i3]);
      num_triangles += 1;

      j1 = i2;
      j2 = i2+n;
      j3 = i3;
      //fprintf (stderr, "       j1 j2 j3 %d %d %d \n", j1, j2, j3); 

      pconn[3*num_triangles+0] = j1;
      pconn[3*num_triangles+1] = j2;
      pconn[3*num_triangles+2] = j3;
      face_normals[num_triangles] = (1.0/3.)*(vertex_normals[j1] + vertex_normals[j2] + 
                                            vertex_normals[j3]);
      num_triangles += 1;
      }

    offset += n;
    }

  // end caps //

  for (i = 0; i < n; i++) {
    i1 = i+1; i2 = i; i3 = num_vertices-2;
    if (i1 == n) i1 = 0;
    pconn[3*num_triangles+0] = i1;
    pconn[3*num_triangles+1] = i2;
    pconn[3*num_triangles+2] = i3;
    v1 = vertices[i2] - vertices[i1];
    v2 = vertices[i3] - vertices[i1];
    norm = v1.cross(v2);
    norm.normalize();
    face_normals[num_triangles] = norm;
    num_triangles += 1;
    //fprintf (stderr, "   >>> i1 i2 i3 %d %d %d \n", i1, i2, i3); 
    //fprintf (stderr, "   >>> fnorm  %f %f %f \n",  norm[0], norm[1], norm[2]);

    j1 = i2 + offset;
    j2 = i1 + offset;
    j3 = num_vertices-1;
    //fprintf (stderr, "       j1 j2 j3 %d %d %d \n", j1, j2, j3); 

    pconn[3*num_triangles+0] = j1;
    pconn[3*num_triangles+1] = j2;
    pconn[3*num_triangles+2] = j3;
    v1 = vertices[j2] - vertices[j1];
    v2 = vertices[j3] - vertices[j1];
    norm = v1.cross(v2);
    norm.normalize();
    face_normals[num_triangles] = norm;
    num_triangles += 1;
    //fprintf (stderr, "       fnorm  %f %f %f \n",  norm[0], norm[1], norm[2]);
    }

  //fprintf (stderr, "   >>> num_triangles[%d] \n", num_triangles); 
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get the mass properties of a cylinder.

void
PmCylinder::getMassProps(PmMassProperties& props)
  {
  //fprintf (stderr, "\n>>>>>> PmCylinder::getMassProps  %s \n", name.c_str());

  if (mass_properties.set) { 
    props = mass_properties;
    return;
    }

  PmMatrix3x3 I, cmat, It;
  float mass, r2, h2, J[3][3], sum;
  PmVector3 com, a1, a2, a3, b1, b2, b3;

  /*
  mass = density * M_PI*radius*radius * length;
  fprintf (stderr, "   >>> analytic mass = %g \n", mass);
  */

#define use_analytic_PmCylinder_getMassProps
#ifdef use_analytic_PmCylinder_getMassProps

  mass = density * M_PI*radius*radius*length;
  r2 = radius*radius;
  h2 = length*length;

  I(0,0) = 0.5*mass*r2;
  I(1,1) = mass*(3*r2 + h2) / 12.0;
  I(2,2) = mass*(3*r2 + h2) / 12.0;

  /*
  fprintf (stderr, "   >>> analytic I = %g \n", I(0,0));
  fprintf (stderr, "                    %g \n", I(1,1));
  fprintf (stderr, "                    %g \n", I(2,2));
  */

  if (mass_properties.set) { 
    props = mass_properties;
    return;
    }

  a1.set(1,0,0); 
  a2.set(0,1,0); 
  a3.set(0,0,1); 

  b1 = axis;
  pm_MathBasisCompute(b1, b2, b3);

  cmat(0,0) = b1*a1;
  cmat(0,1) = b1*a2;
  cmat(0,2) = b1*a3;

  cmat(1,0) = b2*a1;
  cmat(1,1) = b2*a2;
  cmat(1,2) = b2*a3;

  cmat(2,0) = b3*a1;
  cmat(2,1) = b3*a2;
  cmat(2,2) = b3*a3;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum = 0.0;

      for (int k = 0; k < 3; k++) {
        sum += I(i,k)*cmat(k,j);
        }

      It(i,j) = sum;
      }
    }

  cmat.transpose();

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum = 0.0;

      for (int k = 0; k < 3; k++) {
        sum += cmat(i,k)*It(k,j);
        }

      J[i][j] = sum;
      }
    }

  com = center;

  /*
  fprintf (stderr, "   >>> J   = %g %g %g \n", J[0][0], J[0][1], J[0][2]);
  fprintf (stderr, "             %g %g %g \n", J[1][0], J[1][1], J[1][2]);
  fprintf (stderr, "             %g %g %g \n", J[2][0], J[2][1], J[2][2]);
  fprintf (stderr, "----------------------\n");
  */

#else

  genGeometry();

  pm_MathPolyMassProps (num_triangles, connectivity->vals, 3, num_vertices, vertices, 
                        face_normals, density, mass, com, J);

  #define dbg_PmCylinder_getMassProps
  #ifdef dbg_PmCylinder_getMassProps
  fprintf (stderr, "   >>> mass = %g \n", mass);
  fprintf (stderr, "   >>> center = %g %g %g \n", center[0], center[1], center[2]);
  fprintf (stderr, "   >>> com = %g %g %g \n", com[0], com[1], com[2]);
  fprintf (stderr, "   >>> J   = %f %f %f \n", J[0][0], J[0][1], J[0][2]);
  fprintf (stderr, "             %f %f %f \n", J[1][0], J[1][1], J[1][2]);
  fprintf (stderr, "             %f %f %f \n", J[2][0], J[2][1], J[2][2]);
  #endif

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
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a cylinder.

void
PmCylinder::defineRegion(const string name, const bool use_center)
  {
  fprintf (stderr, "\n>>>>>> PmCylinder::defineRegion  name = \"%s\" \n", name.c_str());
  PmMassProperties props;
  PmSolidRegion *region = new PmSolidRegion;
  region->name = name;
  genGeometry();

  for (int i = 0; i < num_vertices; i++) {
    region->coords.push_back(vertices[i]);
    region->coord_indexes.push_back(i);
    region->radii.push_back(0.1);
    }

  getMassProps(props);
  region->center = props.com; 
  regions.push_back(region);
  }

/////////////////////////////////////////////////////////////////
//            s o l i d   s p h e r e                         //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSphere::PmSphere(const string name)
  {
  this->name = name;
  radius = 0.5;
  center.set(0,0,0);
  density = 0.1;

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  graphics_geometry = NULL;

  num_vertices = 0;
  connectivity = NULL;
  num_triangles = 0;
  vertices = NULL;
  vertex_normals = NULL;
  face_normals = NULL;
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a sphere.  

void
PmSphere::defineRegion(const string name, const bool use_center)
  {
  //fprintf (stderr, "\n>>>>>> PmSphere::defineRegion  name = \"%s\" \n", name.c_str());
  PmMassProperties props;
  PmSolidRegion *region = new PmSolidRegion;
  region->name = name;
  genGeometry();
  getMassProps(props);

  if (use_center) {
    region->coords.push_back(props.com);
    region->coord_indexes.push_back(0);
    region->radii.push_back(0.1);
    }
  else {
    for (int i = 0; i < num_vertices; i++) {
      region->coords.push_back(vertices[i]);
      region->coord_indexes.push_back(i);
      region->radii.push_back(0.1);
      }
    }

  getMassProps(props);
  region->center = props.com;
  regions.push_back(region);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display sphere.  

void
PmSphere::display(bool show)
  {
  PmGraphicsAttributes atts;

  genGeometry();

  if (!graphics_geometry) {
    genGeometry();
    graphics_geometry = new PmGraphicsPolygon(name, num_triangles, connectivity, 3,
                                              num_vertices, vertices); 
    }

  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setColor(color);
  graphics_geometry->setAttributes(atts);
  graphics_geometry->display(); 
  }

//*============================================================*
//*==========              getDimensions             ==========*
//*============================================================*
// get sphere dimensions. 

void
PmSphere::getDimensions(vector<float>& dims)
  {
  dims.push_back(radius);
  }

//*============================================================*
//*==========              getExtent                 ==========*
//*============================================================*
// get sphere extent.

void 
PmSphere::getExtent(PmExtent& extent)
  {
  PmVector3 v;

  for (int i = 0; i < 3; i++) {
    v = center;
    v[i] += radius;
    extent.update(v);

    v[i] -= 2.0*radius;
    extent.update(v);
    }
  }

//*============================================================*
//*==========              genGeometry               ==========*
//*============================================================*
// generate sphere geometry.

void
PmSphere::genGeometry()
  {

  if (num_vertices) {
    return;
    }

  PmGcTri tri[1000];
  int num_tri, num_verts, num_iverts, *conn;
  PmVector3 *verts, *iverts, *fnorms, poly[3], *vnorms;

  // generate sphere vertices //

  pm_GcSphereGen (tri, 3, &num_tri);
  //fprintf (stderr, "\n>>>>> PmSphere::num tri = %d \n", num_tri); 

  num_verts = 3 * num_tri;
  verts = new PmVector3[num_verts];
  fnorms = new PmVector3[num_tri];
  int n = 0;

  // compute face normals //

  for (int i = 0; i < num_tri; i++) {
    verts[n++] = radius*tri[i].p1 + center;
    verts[n++] = radius*tri[i].p2 + center;
    verts[n++] = radius*tri[i].p3 + center;
    poly[0] = tri[i].p1;
    poly[1] = tri[i].p2;
    poly[2] = tri[i].p3;
    pm_GcPolyNormComp(3, poly, fnorms[i]);
    }

  // convert into indexed polygons //

  pm_GcTriGetIndexes (num_tri, verts, &conn, &num_iverts, &iverts);

  // compute vertex normals //

  pm_GcPolyVertexNormsComp(num_tri, conn, 0, num_iverts, iverts, fnorms, &vnorms);

  // set object data //

  num_triangles = num_tri;
  num_vertices = num_iverts;
  vertices = iverts; 
  vertex_normals = vnorms; 
  face_normals = fnorms;

  connectivity = new PmConn(3*num_tri);
  int *ivals = connectivity->vals;
  n = 0;

  for (int i = 0; i < num_tri; i++) {
    ivals[n++] = conn[4*i+1]; 
    ivals[n++] = conn[4*i+2]; 
    ivals[n++] = conn[4*i+3]; 
    }
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get the mass properties of a sphere.  

void
PmSphere::getMassProps(PmMassProperties& props)
  {

  if (mass_properties.set) { 
    props = mass_properties;
    return;
    }

  PmMatrix3x3 I, cmat, It;
  float mass, r2;
  PmVector3 com, a1, a2, a3, b1, b2, b3;

  mass = density * (4.0/3.0)*M_PI*radius*radius*radius;
  r2 = radius*radius;

  I(0,0) = 2.0/5.0*mass*r2;
  I(1,1) = I(0,0); 
  I(2,2) = I(0,0); 

  mass_properties.set = true;
  mass_properties.mass = mass;
  mass_properties.com = center;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      mass_properties.inertia(i,j) = I(i,j);
      }
    }

  props = mass_properties;
  }

/////////////////////////////////////////////////////////////////
//            s o l i d   e l l i p s o i d                   //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmEllipsoid::PmEllipsoid(const string name)
  {
  this->name = name;
  center.set(0,0,0);
  density = 0.1;
  PmVector3 v; 

  for (int i = 0; i < 3; i++) {
    radii.push_back(1.0);
    }

  v.set(1,0,0); axes.push_back(v); 
  v.set(0,1,0); axes.push_back(v); 
  v.set(0,0,1); axes.push_back(v); 

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  graphics_geometry = NULL;

  num_vertices = 0;
  connectivity = NULL;
  num_triangles = 0;
  vertices = NULL;
  vertex_normals = NULL;
  face_normals = NULL;
  }

//*============================================================*
//*==========              defineRegion              ==========*
//*============================================================*
// define a region for a sphere.  

void
PmEllipsoid::defineRegion(const string name, const bool use_center)
  {
  //fprintf (stderr, "\n>>>>> PmEllipsoid::defineRegion  name = \"%s\" \n", name.c_str());
  PmMassProperties props;
  PmSolidRegion *region = new PmSolidRegion;
  region->name = name;
  genGeometry();
  getMassProps(props);

  if (use_center) {
    region->coords.push_back(props.com);
    region->coord_indexes.push_back(0);
    region->radii.push_back(0.1);
    }
  else {
    for (int i = 0; i < num_vertices; i++) {
      region->coords.push_back(vertices[i]);
      region->coord_indexes.push_back(i);
      region->radii.push_back(0.1);
      }
    }

  getMassProps(props);
  region->center = props.com;
  regions.push_back(region);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display ellipsoid.

void
PmEllipsoid::display(bool show)
  {
  PmGraphicsAttributes atts;

  genGeometry();

  if (!graphics_geometry) {
    genGeometry();
    graphics_geometry = new PmGraphicsPolygon(name, num_triangles, connectivity, 3,
                                              num_vertices, vertices); 
    }

  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setColor(color);
  graphics_geometry->setAttributes(atts);
  graphics_geometry->display(); 
  }

//*============================================================*
//*==========              setAxes                   ==========*
//*============================================================*
// set ellipsoid axes.

void 
PmEllipsoid::setAxes(const vector<PmVector3>& axes)
  {
  this->axes = axes;
  }

//*============================================================*
//*==========              getDimensions             ==========*
//*============================================================*
// get ellipsoid dimensions. 

void
PmEllipsoid::getDimensions(vector<float>& dims)
  {
  dims = radii; 
  }

//*============================================================*
//*==========              getExtent                 ==========*
//*============================================================*
// get ellipsoid extent.

void 
PmEllipsoid::getExtent(PmExtent& extent)
  {
  PmVector3 v;

  for (int i = 0; i < 3; i++) {
    v = center + radii[i]*axes[i];
    extent.update(v);

    v = center - radii[i]*axes[i];
    extent.update(v);
    }
  }

//*============================================================*
//*==========              genGeometry               ==========*
//*============================================================*
// generate ellipsoid geometry.

void
PmEllipsoid::genGeometry()
  {

  if (num_vertices) {
    return;
    }

  PmGcTri tri[1000];
  int num_tri, num_verts, num_iverts, *conn;
  PmVector3 *verts, *iverts, *fnorms, poly[3], *vnorms;

  // generate vertices //

  pm_GcSphereGen (tri, 3, &num_tri);

  num_verts = 3 * num_tri;
  verts = new PmVector3[num_verts];
  fnorms = new PmVector3[num_tri];
  int n = 0;

  // compute face normals //

  for (int i = 0; i < num_tri; i++) {
    verts[n++] = tri[i].p1;
    verts[n++] = tri[i].p2;
    verts[n++] = tri[i].p3;
    poly[0] = tri[i].p1;
    poly[1] = tri[i].p2;
    poly[2] = tri[i].p3;
    pm_GcPolyNormComp(3, poly, fnorms[i]);
    }

  // convert into indexed polygons //

  pm_GcTriGetIndexes (num_tri, verts, &conn, &num_iverts, &iverts);

  // compute vertex normals //

  pm_GcPolyVertexNormsComp(num_tri, conn, 0, num_iverts, iverts, fnorms, &vnorms);

  // set object data //

  num_triangles = num_tri;
  num_vertices = num_iverts;
  vertices = iverts; 
  vertex_normals = vnorms; 
  face_normals = fnorms;

  connectivity = new PmConn(3*num_tri);
  int *ivals = connectivity->vals;
  n = 0;

  for (int i = 0; i < num_tri; i++) {
    ivals[n++] = conn[4*i+1]; 
    ivals[n++] = conn[4*i+2]; 
    ivals[n++] = conn[4*i+3]; 
    }

  // transform sphere geometry into ellipsoid //

  float a  = radii[0];
  float b  = radii[1];
  float c  = radii[2];

  PmVector3 u = axes[0];
  PmVector3 v = axes[1];
  PmVector3 w = axes[2];
  PmVector3 p;
  float cx, cy, cz, dpu, dpv, dpw;
  cx = center[0];
  cy = center[1];
  cz = center[2];

  for (int i = 0; i < num_iverts; i++) {
    p = iverts[i];
    dpu = p * u;
    dpv = p * v;
    dpw = p * w;
    iverts[i][0] = a*dpu*u[0] + b*dpv*v[0] + c*dpw*w[0] + cx;
    iverts[i][1] = a*dpu*u[1] + b*dpv*v[1] + c*dpw*w[1] + cy;
    iverts[i][2] = a*dpu*u[2] + b*dpv*v[2] + c*dpw*w[2] + cz;
    }
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get the mass properties of an ellipsoid.  

void
PmEllipsoid::getMassProps(PmMassProperties& props)
  {

  if (mass_properties.set) { 
    props = mass_properties;
    return;
    }

  PmMatrix3x3 I, cmat, It;
  float mass, r2, a, b, c, sum, J[3][3];
  PmVector3 com, a1, a2, a3, b1, b2, b3;

#define use_analytic_PmEllipsoid_getMassProps
#ifdef use_analytic_PmEllipsoid_getMassProps

  a = radii[0];
  b = radii[1];
  c = radii[2];

  mass = (4.0/3.0)*density*M_PI*a*b*c;

  I(0,0) = mass*(b*b + c*c) / 5.0;;
  I(1,1) = mass*(a*a + c*c) / 5.0;;
  I(2,2) = mass*(b*b + a*a) / 5.0;;

  a1.set(1,0,0); 
  a2.set(0,1,0); 
  a3.set(0,0,1); 

  b1 = axes[0];
  b2 = axes[1];
  b3 = axes[2];

  cmat(0,0) = b1*a1;
  cmat(0,1) = b1*a2;
  cmat(0,2) = b1*a3;

  cmat(1,0) = b2*a1;
  cmat(1,1) = b2*a2;
  cmat(1,2) = b2*a3;

  cmat(2,0) = b3*a1;
  cmat(2,1) = b3*a2;
  cmat(2,2) = b3*a3;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum = 0.0;

      for (int k = 0; k < 3; k++) {
        sum += I(i,k)*cmat(k,j);
        }

      It(i,j) = sum;
      }
    }

  cmat.transpose();

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum = 0.0;

      for (int k = 0; k < 3; k++) {
        sum += cmat(i,k)*It(k,j);
        }

      J[i][j] = sum;
      }
    }

  com = center;

  /*
  fprintf (stderr, "   >>> J   = %g %g %g \n", J[0][0], J[0][1], J[0][2]);
  fprintf (stderr, "             %g %g %g \n", J[1][0], J[1][1], J[1][2]);
  fprintf (stderr, "             %g %g %g \n", J[2][0], J[2][1], J[2][2]);
  fprintf (stderr, "----------------------\n");
  */

#else

  genGeometry();

  pm_MathPolyMassProps (num_triangles, connectivity->vals, 3, num_vertices, vertices, 
                        face_normals, density, mass, com, J);

  #define dbg_PmCylinder_getMassProps
  #ifdef dbg_PmCylinder_getMassProps
  fprintf (stderr, "   >>> mass = %g \n", mass);
  fprintf (stderr, "   >>> center = %g %g %g \n", center[0], center[1], center[2]);
  fprintf (stderr, "   >>> com = %g %g %g \n", com[0], com[1], com[2]);
  fprintf (stderr, "   >>> J   = %f %f %f \n", J[0][0], J[0][1], J[0][2]);
  fprintf (stderr, "             %f %f %f \n", J[1][0], J[1][1], J[1][2]);
  fprintf (stderr, "             %f %f %f \n", J[2][0], J[2][1], J[2][2]);
  #endif

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
//*==========              getParameters             ==========*
//*============================================================*
// get the parameters of an ellipsoid.

void
PmEllipsoid::getParameters(PmVector3& center, vector<float>& radii, 
                           vector<PmVector3>& axes)
  {
  center = this->center;
  radii = this->radii;
  axes = this->axes;
  }

//*============================================================*
//*==========              getParameters             ==========*
//*============================================================*

void 
PmEllipsoid::setRadii(const vector<float>& radii)
  {
  this->radii = radii;
  }


//*============================================================*
//*==========              buildGeomName             ==========*
//*============================================================*

void
PmSolid::buildGeomName (string type, string& grname) {
  grname = "solid";
  grname += '[' + name + ']';

  if (type.size()) {
    grname += type;
    }
  }

//*============================================================*
//*==========             addGeometry                ==========*
//*============================================================*
// add graphics geometry object.

void
PmSolid::addGeometry (PmGraphicsGeometry *geom) {
  graphics_geometries.push_back (geom);
  }

//*============================================================*
//*==========               getGeometry              ==========*
//*============================================================*
// get a graphics geometry object.

void
PmSolid::getGeometry (string name, PmGraphicsGeometry **geom)
  {
  *geom = NULL;

  for (int i = 0; i < graphics_geometries.size(); i++) {
    if (name == graphics_geometries[i]->name) {
      *geom = graphics_geometries[i];
      }
    }
  }

}

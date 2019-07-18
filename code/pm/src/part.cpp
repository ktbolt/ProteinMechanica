
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
 * part:             p a r t i c l e   s y s t e m            *
 *============================================================*/

#include "part.h"

namespace ProteinMechanica {


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmParticle::PmParticle(const string name, int num_coords, PmVector3 *coords) 
  {
  //fprintf (stderr, "\n>>>>>> PmParticle  ctor  name  %s \n", name.c_str());
  this->name = name;
  this->obj_type = PM_OBJECT_PARTICLE;
  this->num_coordinates = num_coords;
  this->coordinates = coords;
  this->masses = new float[num_coords];
  this->radii = new float[num_coords];

  for (int i = 0; i < num_coords; i++) {
    this->masses[i] = 1.0;
    this->radii[i] = 0.2;
    }

  color.set(1,1,1);
  display_spheres = false;
  map_mass = false;
  map_colors = NULL;
  }

PmParticle::~PmParticle() {
  }

//*============================================================*
//*==========          setColor                      ==========*
//*============================================================*

void
PmParticle::setColor(const PmVector3& color) {
  this->color = color;
  }

//*============================================================*
//*==========          get/set name                  ==========*
//*============================================================*

void 
PmParticle::getName (string& name) { 
  name = this->name;  
  }

void 
PmParticle::setName (const string name) { 
  this->name = name;  
  }

//*============================================================*
//*==========          setMapMass                    ==========*
//*============================================================*

void
PmParticle::setMapMass (const bool val)
  {
  map_mass = val;
  }

//*============================================================*
//*==========          getMassProps                  ==========*
//*============================================================*

void 
PmParticle::getMassProps (PmMassProperties& props)
  {
  float mass;
  PmMatrix3x3 inertia;
  float cmx, cmy, cmz, cx, cy, cz;
  float x, y, z, xx, xy, xz, yy, yz, zz;

  float scale = 1.0 / pmSystem.getUnitsMassScale();
  float total_mass = 0.0;
  cmx = cmy = cmz = 0.0;
  cx = cy = cz = 0.0;

  //  first compute center of mass //

  for (int i = 0; i < num_coordinates; i++) {
    x = coordinates[i][0];
    y = coordinates[i][1];
    z = coordinates[i][2];
    mass = scale * masses[i];
    total_mass += mass;
    cmx += x*mass;
    cmy += y*mass;
    cmz += z*mass;
    cx += x;
    cy += y;
    cz += z;
    }

  cmx /= total_mass; cmy /= total_mass; cmz /= total_mass;
  cx  /= num_coordinates;  cy  /= num_coordinates;  cz  /= num_coordinates;

  //  compute inertia tensor //

  xx = xy = xz = yy = yz = zz = 0.0;

  for (int i = 0; i < num_coordinates; i++) {
    x = coordinates[i][0] - cmx;
    y = coordinates[i][1] - cmy;
    z = coordinates[i][2] - cmz;
    mass = scale * masses[i];
    xx += mass*x*x;
    xy += mass*x*y;
    xz += mass*x*z;
    yy += mass*y*y;
    yz += mass*y*z;
    zz += mass*z*z;
    }

  inertia(0,0) = yy + zz;
  inertia(0,1) = -xy;
  inertia(0,2) = -xz;
  inertia(1,0) = -xy;
  inertia(1,1) = xx + zz;
  inertia(1,2) = -yz;
  inertia(2,0) = -xz;
  inertia(2,1) = -yz;
  inertia(2,2) = xx + yy;

  props.mass = total_mass;
  props.com.set(cmx, cmy, cmz);
  props.inertia = inertia;
  }

//*============================================================*
//*==========          setMasses                     ==========*
//*============================================================*

void 
PmParticle::setMasses(float *masses)
  {
  /*
  float mass, vmin, vmax;
  float maxr, r, dr; 
  int n = 10;

  vmax = vmin = masses[0];

  for (int i = 0; i < num_coordinates; i++) {
    this->masses[i] = masses[i];
    mass = masses[i];
    if (mass > vmax) vmax = mass;
    if (mass < vmin) vmin = mass;
    }

  maxr = 0.2;
  dr = maxr / (vmax - vmin);

  for (int i = 0; i < num_coordinates; i++) {
    mass = this->masses[i];
    this->radii[i] = (mass - vmin)*dr; 
    }
  */

  for (int i = 0; i < num_coordinates; i++) {
    this->masses[i] = masses[i];
    }
  }

//*============================================================*
//*==========          display                       ==========*
//*============================================================*

void 
PmParticle::display(const bool show) 
  {
  fprintf (stderr, "\n>>>>>> PmParticle  display name  %s \n", name.c_str());

  if (!num_coordinates) {
    return;
    }

  string geom_name, sgeom_name;
  PmGraphicsAttributes atts;

  atts.setColor(color);
  atts.setMarker(false);
  atts.setScale(2.0);
  atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);

  if (!point_geom) {
    geom_name = "grid[" + name + ']' + "points";
    point_geom = new PmGraphicsPoint(geom_name, num_coordinates, coordinates);

    if (map_mass) {
      float mass, vmin, vmax;
      vmin = vmax = masses[0];
      map_colors = new PmVector3[num_coordinates];

      for (int i = 0; i < num_coordinates; i++) {
        mass = masses[i];
        if (mass > vmax) vmax = mass;
        if (mass < vmin) vmin = mass;
        }

      PmGraphicsInterface::mapDataToColors(vmin, vmax, num_coordinates, masses,
                                            map_colors); 
      point_geom->setColors(num_coordinates, map_colors);
      }

    }

  point_geom->setAttributes(atts); 
  point_geom->display(); 

  if (display_spheres) { 
    if (!sphere_geom) { 
      sgeom_name = "grid[" + name + ']' + "spheres";
      sphere_geom = new PmGraphicsSphere(sgeom_name, num_coordinates, coordinates, radii);

      if (map_colors) {
        sphere_geom->setColors(num_coordinates, map_colors);
        }
      }
    }

  if (sphere_geom) { 
    sphere_geom->setAttributes(atts);
    sphere_geom->display();
    }
  }

//*============================================================*
//*==========          setDisplaySpheres             ==========*
//*============================================================*

void
PmParticle::setDisplaySpheres(const bool spheres) {
  display_spheres = spheres;
  }

//*============================================================*
//*==========              getXform                  ==========*
//*============================================================*

void
PmParticle::getXform (PmXform& xform) {
  xform = this->xform;
  }

//*============================================================*
//*==========              setXform                  ==========*
//*============================================================*
// set the transformation for a surface.

void
PmParticle::setXform (PmXform& xform)
  {
  this->xform = xform;
  PmGraphicsGeometry *geom;

  if (point_geom) {  
    geom = point_geom;
    geom->setXform(xform);
    geom->display();
    }

  if (sphere_geom) {  
    geom = sphere_geom;
    geom->setXform(xform);
    geom->display();
    }
  }

}



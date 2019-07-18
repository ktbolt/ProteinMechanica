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
//* rest:                     r e s t r a i n t                *
//*============================================================*

#include "rest.h"
#include "sobj.h"
#include "pobj.h"
#include "graphics.h"


//===== debug symbols =====//
#define ndbg_PmRestraint_compForces

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmRestraint::PmRestraint(const string name) 
  {
  this->name = name;
  type = PM_RESTRAINT_CENTER;
  initialized = false;
  line_geometry = NULL;
  active = false;
  has_power_params = false;
  power_params[0] = 0.0;
  power_params[1] = 0.0;
  distance = 0.0;
  use_absolute_distance = false;

  has_ramp = false;
  ramp_init = true;
  ramp_start = 0.0; 
  ramp_time = 0.0;
  dramp = 0.0; 
  compute_energy = false;
  current_energy = 0.0;
  current_residual = 0.0;

  color.set(1,0,1);
  width = 1;
  show = false;
  }

//*============================================================*
//*==========              convType                  ==========*
//*============================================================*
// convert a string restraint type to a symbol.

void
PmRestraint::convType(const string str, PmRestraintType& type)
  {
  type = PM_RESTRAINT_UNKNOWN;

  if (str == "all") {
    type = PM_RESTRAINT_ALL;
    }
  else if (str == "center") {
    type = PM_RESTRAINT_CENTER;
    }
  }

//*============================================================*
//*==========              setType                   ==========*
//*============================================================*
// set restraint type.

void
PmRestraint::setType(const PmRestraintType type) {
  this->type = type;
  }

//*============================================================*
//*==========               hasEnergy                ==========*
//*============================================================*

bool 
PmRestraint::hasEnergy() {
  return compute_energy;
  }

//*============================================================*
//*==========               setCompEnergy            ==========*
//*============================================================*

void 
PmRestraint::setCompEnergy(const bool flag) {
  compute_energy = flag;
  }

//*============================================================*
//*==========               getEnergy                ==========*
//*============================================================*

void 
PmRestraint::getEnergy(float& energy) { 
  energy = this->current_energy; 
  }

//*============================================================*
//*==========               getName                  ==========*
//*============================================================*

void 
PmRestraint::getName(string& name) { 
  name = this->name; 
  }

//*============================================================*
//*==========               getResidual              ==========*
//*============================================================*

void
PmRestraint::getResidual(float& residual)
  {
  residual = current_residual;
  }

//*============================================================*
//*==========               isActive                 ==========*
//*============================================================*
// check if a force is active for the given time.

bool 
PmRestraint::isActive(float t)
  {
  //fprintf (stderr, "\n>>>>>>> PmRestraint::isActive  %s t[%f] \n", name.c_str(), t);
  bool has_time = time_interval.hasTime(t);
  //fprintf (stderr, "   >>>> has time[%d] \n", has_time);

  // if the restraint is active then deactivate it //

  if (active && !has_time) { 
    if (line_geometry && show) {
      PmGraphicsAttributes atts;
      atts.setVisible(false);
      PmGraphicsLine *line = dynamic_cast<PmGraphicsLine*>(line_geometry);
      line->setAttributes(atts);
      line->display();
      PmGraphicsPoint *point = dynamic_cast<PmGraphicsPoint*>(point_geometry);
      point->setAttributes(atts);
      point->display();
      }

    active = false;
    }

  else if (!active && has_time) { 
    if (line_geometry && show) {
      PmGraphicsAttributes atts;
      atts.setVisible(true);
      PmGraphicsLine *line = dynamic_cast<PmGraphicsLine*>(line_geometry);
      line->setAttributes(atts);
      line->display();
      PmGraphicsPoint *point = dynamic_cast<PmGraphicsPoint*>(point_geometry);
      point->setAttributes(atts);
      point->display();
      }

    active = true;
    }

  return (active);
  }

//*============================================================*
//*==========               compForces               ==========*
//*============================================================*

void 
PmRestraint::compForces(const bool upd, const float time)
  {
  #ifdef dbg_PmRestraint_compForces
  fprintf (stderr, ">>>>>> PmRestraint::compForces[%s] \n", name.c_str());
  #endif

  bool has_time = time_interval.hasTime(time);

  if (!has_time) {
    this->current_energy = 0.0;
    return;
    }

  if (!initialized) {
    this->init();
    }

  if (points1.size() == 0) {
    return;
    }

  //===== transform restraint points =====//

  this->xformPoints();

  //===== compute force on bodies due to restraint spring =====//

  PmVector3 dl, dir, force;
  float d, ds, d0, k, dist, fmag, total_energy, residual;
  float pmax, pmin;
  float rampf;
  bool abs_dist;

  k = this->force_const;
  dist = this->distance;
  total_energy = 0.0;
  residual = 0.0;

  for (unsigned int i = 0; i < points1.size(); i++) {
    #ifdef dbg_PmRestraint_compForces
    fprintf (stderr, "----- %d -----\n", i);
    #endif
    dl = points1[i] - points2[i];
    d0 = dl.length();
    dir = current_points1[i] - current_points2[i];
    d = dir.length();
    residual += d;
    fmag = 0.0;

    #ifdef dbg_PmRestraint_compForces
    fprintf (stderr, ">>> d0=%g \n", d0);
    fprintf (stderr, ">>> dir=%g %g %g\n", dir[0], dir[1], dir[2]);
    #endif

    if (use_absolute_distance) {
      ds = dist;
      }
    else {
      ds = d0*(1.0 - dist);
      }

    if (this->has_power_params) { 
      pmin = this->power_params[0];
      pmax = this->power_params[1];
      ds = exp(-10.0*d/d0);

      if (pmax*ds < pmin) {
        fmag = -pmin;
        }
      else {
        fmag = -pmax*ds;
        }
      }
    else {
      if (d == 0.0) {
        fmag = -k;
        }
      else {
        fmag = k*(ds / d - 1.0);
        }

      if (this->getRampParams(time, rampf)) {
        fmag = fmag * rampf;
        }
      }

    #ifdef dbg_PmRestraint_compForces
    fprintf (stderr,  ">>> fmag=%g d=%g ds=%g d-ds=%g \n", fmag, d, ds, d-ds);
    #endif
    force = fmag*dir;
    this->sobjs[0]->addForce(current_points1[i], force);
    force = -force;
    this->sobjs[1]->addForce(current_points2[i], force);

    if (compute_energy) { 
      //total_energy += 0.5*k*(d-ds)*(d-ds);
      total_energy += 0.5*k*d*d;
      }

    this->current_residual = residual;
    }

  if (compute_energy) { 
    this->current_energy = total_energy;
    }

  //===== show the restraint =====//

  this->display();
  }

//*============================================================*
//*==========               xformPoints              ==========*
//*============================================================*

void
PmRestraint::xformPoints()
  {
  PmXform xform;
  PmPhysicalObj *pobj;
  PmSimulationObj *sobj;

  sobj = sobjs[0];
  sobj->getPhysicalObj(&pobj);
  pobj->getXform(xform);
  PmMatrix3x3 mat = xform.matrix;

  for (unsigned int i = 0; i < points1.size(); i++) {
    current_points1[i] = mat*(points1[i]-xform.center) + xform.translation + xform.center;
    }

  sobj = sobjs[1];
  sobj->getPhysicalObj(&pobj);
  pobj->getXform(xform);
  mat = xform.matrix;

  for (unsigned int i = 0; i < points2.size(); i++) {
    current_points2[i] = mat*(points2[i]-xform.center) + xform.translation + xform.center;
    }
  }

//*============================================================*
//*==========               get/setSimObjs           ==========*
//*============================================================*

void 
PmRestraint::getSimObjs(PmSimulationObj **sobj1, PmSimulationObj **sobj2) 
  {
  *sobj1 = this->sobjs[0];
  *sobj2 = this->sobjs[1];
  }

void 
PmRestraint::setSimObjs (PmSimulationObj *sobj1, PmSimulationObj *sobj2)
   {
   this->sobjs[0] = sobj1;
   this->sobjs[1] = sobj2;
   }

//*============================================================*
//*==========               setTime                  ==========*
//*============================================================*

void 
PmRestraint::setTime (PmTimeInterval& time) { 
  this->time_interval = time; 
  }

//*============================================================*
//*==========               getPointDistance         ==========*
//*============================================================*
// get current distance between restraint point.

void 
PmRestraint::getPointDistance(float& dist, float& pdist)
  {
  dist = current_residual;
  }

//*============================================================*
//*==========               getRampParams            ==========*
//*============================================================*
// get ramp params.

bool 
PmRestraint::getRampParams(const float time, float& dramp)
  {
  dramp = 0.0;

  if (!has_ramp) {
    return false;
    }

  if (ramp_init) {
    ramp_start = time;
    ramp_init = false;
    }

  else if (ramp_start + ramp_time < time) {
    return false;
    }

  dramp = (time - ramp_start) / ramp_time;
  return true;
  }

//*============================================================*
//*==========               setRampParams            ==========*
//*============================================================*
// set ramp params.

void
PmRestraint::setRampParams(const float param)
  {

  if (param == 0.0) {
    return;
    }

  ramp_time = param;
  ramp_init = true;
  has_ramp = true;
  }

//*============================================================*
//*==========               setRegions               ==========*
//*============================================================*
// set regions.             

void 
PmRestraint::setRegions(const string rgn1, const string rgn2) {
  this->regions[0] = rgn1;
  this->regions[1] = rgn2;
  }

//*============================================================*
//*==========               init                     ==========*
//*============================================================*

void
PmRestraint::init()
  {
  #ifdef dbg_PmRestraint_init
  fprintf (stderr, ">>>>>> PmRestraint::init[%s] \n", name.c_str());
  fprintf (stderr, ">>> type=%d \n", type); 
  #endif
  PmPhysicalObj *pobj;
  PmSimulationObj *sobj;
  PmVector3 center;

  if (!sobjs[0] || !sobjs[1]) {
    pm_ErrorWarnReport (PM, "restraint \"%s\" doesn't have simulation objects.", "*",
                        name.c_str());
    return;
    }

  //===== compute the center of each region =====//

  if (type == PM_RESTRAINT_CENTER) {
    #ifdef dbg_PmRestraint_init
    fprintf (stderr, ">>> type=center \n");
    #endif
    sobjs[0]->getPhysicalObj(&pobj);
    pobj->getRegionCenter(regions[0], center);
    points1.push_back(center);
    current_points1.push_back(center);

    sobjs[1]->getPhysicalObj(&pobj);
    pobj->getRegionCenter(regions[1], center);
    points2.push_back(center);
    current_points2.push_back(center);
    }

  else if (type == PM_RESTRAINT_ALL) {
    #ifdef dbg_PmRestraint_init
    fprintf (stderr, ">>> type=all\n");
    #endif
    vector<PmVector3> coords1, coords2;
    sobjs[0]->getRegionCoords(regions[0], coords1);
    sobjs[1]->getRegionCoords(regions[1], coords2);

    if (coords1.size() != coords2.size()) {
      pm_ErrorWarnReport (PM, "restraint \"%s\" regions don't have the same number of coordinates.", "*", name.c_str());
      pm_ErrorWarnReport (PM, "rgn 1 coords size %d ", "*", coords1.size()); 
      pm_ErrorWarnReport (PM, "rgn 2 coords size %d ", "*", coords2.size()); 
      return;
      }

    #ifdef dbg_PmRestraint_init
    fprintf (stderr, ">>> coords1.size=%d \n", coords1.size());
    #endif

    for (unsigned int i = 0; i < coords1.size(); i++) {
      points1.push_back(coords1[i]);
      current_points1.push_back(coords1[i]);
      points2.push_back(coords2[i]);
      current_points2.push_back(coords2[i]);
      }
    }

  initialized = true;
  }

//*============================================================*
//*==========               setPowerParams           ==========*
//*============================================================*
// set power params.       

void 
PmRestraint::setPowerParams(const float params[2])
  {
  power_params[0] = params[0];
  power_params[1] = params[1];
  has_power_params = true; 
  }

//*============================================================*
//*==========               getPowerParams           ==========*
//*============================================================*
// get power params.

bool 
PmRestraint::getPowerParams(float params[2])
  {
  if (!has_power_params) { 
    params[0] = 0; 
    params[1] = 0; 
    }
  else {
    params[0] = power_params[0];
    params[1] = power_params[1];
    }

  return has_power_params;
  } 

//*============================================================*
//*==========               get/setSpringData        ==========*
//*============================================================*

void 
PmRestraint::getSpringData(float *k, float *dist) 
  {
  *k = this->force_const;
  *dist = this->distance;
  }

void 
PmRestraint::setSpringData(const float k, const float dist) 
  {
  this->force_const = k;
  this->distance = dist;
  }

//*============================================================*
//*==========           getAbsoluteDistance          ==========*
//*============================================================*

void 
PmRestraint::getAbsoluteDistance(bool& flag) {
  flag = this->use_absolute_distance;
  }

void 
PmRestraint::setAbsoluteDistance(const bool flag) {
  this->use_absolute_distance = flag;
  }

//*============================================================*
//*==========                 setColor               ==========*
//*============================================================*

void 
PmRestraint::setColor(const PmVector3 color) {
  this->color = color;
  }

//*============================================================*
//*==========                setWidth                ==========*
//*============================================================*

void 
PmRestraint::setWidth(const float val) {
  this->width = val;
  }

//*============================================================*
//*==========                setShow                 ==========*
//*============================================================*

void 
PmRestraint::setShow(const bool val) {
  this->show = val;
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a restraint.

void 
PmRestraint::display()
  {
  PmGraphicsLine *line;
  PmGraphicsPoint *point;
  PmGraphicsAttributes patts;
  int n = current_points1.size();

  if (!line_geometry) {
    string geom_name;
    geom_name = "force[" + name + "]";
    PmVector3 *verts = new PmVector3[2*n];

    for (unsigned int i = 0; i < current_points1.size(); i++) {
      verts[2*i] = current_points1[i];
      verts[2*i+1] = current_points2[i];
      }

    line = new PmGraphicsLine(geom_name, 2*n, verts);
    patts.setVisible(this->show);
    patts.setColor(this->color);
    patts.setLineWidth(this->width);
    patts.setDisjoint(true);
    line->setAttributes(patts);
    line_geometry = dynamic_cast<PmGraphicsGeometry*>(line);
    line->display();

    patts.setScale(0.1);
    patts.setMarker(false);
    point = new PmGraphicsPoint(geom_name, 2*n, verts);
    point->setAttributes(patts);
    point_geometry = dynamic_cast<PmGraphicsGeometry*>(point);
    point->display();
    geometry_vertices = verts;
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(line_geometry);
    point = dynamic_cast<PmGraphicsPoint*>(point_geometry);

    for (unsigned int i = 0; i < current_points1.size(); i++) {
      geometry_vertices[2*i] = current_points1[i];
      geometry_vertices[2*i+1] = current_points2[i];
      }

    line->update(2*n, geometry_vertices);
    point->update(2*n, geometry_vertices);
    }
  }

//*============================================================*
//*==========          getRegionCenter               ==========*
//*============================================================*
// get the center of a region for a restraint.

void 
PmRestraint::getRegionCenter(const int rgn_num, PmVector3& center)
  {
  string rgn;
  PmPhysicalObj *pobj;
  PmSimulationObj *sobj;

  if (rgn_num == 1) {
    rgn = regions[0];
    sobj = sobjs[0];
    }
  else {
    rgn = regions[1];
    sobj = sobjs[1];
    }
 
  if (!sobj) {
    return;
    }

  sobj->getPhysicalObj(&pobj);
  pobj->getRegionCenter(rgn, center);
  }

}





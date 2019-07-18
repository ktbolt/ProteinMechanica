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
//* force:                    f o r c e                        *
//*============================================================*

#include "force.h"
#include "graphics.h"
#include "body.h"
#include "pm/mth.h"

namespace ProteinMechanica {

// debug symbols //

#define ndbg_PmExplicitForce_display

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmForce::PmForce()
  {
  has_ramp = false;
  ramp_init = true;
  ramp_start = 0.0;
  ramp_time = 0.0;
  dramp = 0.0;
  show = true;
  color.set(1,0,1);
  }

//*============================================================*
//*==========            convForceType               ==========*
//*============================================================*

void
PmForce::convForceType(const string str, PmForceType& type)
  {
  if (str == "explicit") {
    type = PM_FORCE_EXPLICIT;
    }
  else if (str == "random") {
    type = PM_FORCE_RANDOM;
    }
  else {
    type = PM_FORCE_UNKNOWN;
    }
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a simulation object of a particular type.

PmForce*
PmForce::create(const string name, const PmForceType type)
  {
  switch (type) {
    case PM_FORCE_EXPLICIT:
      {
      PmExplicitForce *eforce = new PmExplicitForce(name);
      eforce->type = type;
      return eforce;
      }
    break;

    case PM_FORCE_RANDOM:
      {
      PmRandomForce *rforce = new PmRandomForce(name);
      rforce->type = type;
      return rforce;
      }
    break;

    default:
      return NULL;
    break;
    }
  }

//*============================================================*
//*==========              procQuery                 ==========*
//*============================================================*

void
PmForce::procQuery(PmQuery& query)
  {
  fprintf (stderr, ">>>>>> PmForce::procQuery \n"); 
  string force_name, stype;
  unsigned int n, n1, n2, n3;
  //vector<PmPotential*> pot_list;
  PmForce *force;
  //PmPotentialGeom *pot_geom;
  float energy;
  PmForceType type;

  n = query.name.length();
  n1 = query.name.find("[");
  n2 = query.name.find("]");

  if ((n1 == string::npos) || (n2 == string::npos)) {
    return;
    }

  force_name.assign(query.name, n1+1, n2-n1-1);
  fprintf (stderr, ">>> force_name=\"%s\" \n", force_name.c_str());

  pmSystem.getForce(force_name, &force);

  if (!force) {
    return;
    }

  type = force->getType();

  if (PM_FORCE_EXPLICIT) {
    fprintf (stderr, ">>> force type=explicit \n");
    PmExplicitForce *eforce = dynamic_cast<PmExplicitForce*>(force);
    eforce->procAxesSelect(query);
    }
  }

//*============================================================*
//*==========               getRampParams            ==========*
//*============================================================*
// get ramp params for the given time.

bool
PmForce::getRampParams(const float time, float& dr)
  {
  dr = 0.0;

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

  dr = (time - ramp_start) / ramp_time;
  dramp = dr; 
  return true;
  }

//*============================================================*
//*==========               setRampParams            ==========*
//*============================================================*
// set ramp parameters.

void 
PmForce::setRampParams(const float param)
  {
  if (param == 0.0) {
    return;
    }

  ramp_time = param;
  ramp_init = true;
  has_ramp = true;
  }

//*============================================================*
//*==========               setColor                 ==========*
//*============================================================*
// set color.             

void
PmForce::setColor(PmVector3& color)
  {
  this->color = color;
  }

//*============================================================*
//*==========               setShow                  ==========*
//*============================================================*
// set show.

void
PmForce::setShow(bool show)
  {
  this->show = show;
  }

//*============================================================*
//*==========               scale                    ==========*
//*============================================================*
// set display scale.

void
PmForce::setScale(const float val) {
  this->scale = val;
  }


///////////////////////////////////////////////////////////////
//       e x p l i c i t   f o r c e  s                     //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmExplicitForce::PmExplicitForce(const string name) 
  {
  this->name = name;
  scale = 1.0;
  sobj = NULL;
  line_geometry = NULL;
  axes_line_geometry = NULL;
  axes_sphere_geometry = NULL;
  picking_sphere_geometry = NULL;
  point_geometry = NULL;
  active = false;
  global_frame = false;
  torque = false;
  interactive = false;
  strength = 1.0;
  use_picking_sphere = false;
  }

//*============================================================*
//*==========               isActive                 ==========*
//*============================================================*
// check if a force is active for the given time.

bool 
PmExplicitForce::isActive(float t)
  {
  bool has_time = time_interval.hasTime(t);
  //fprintf (stderr, "\n>>>>>>> PmExplicitForce::isActive  t[%f] \n", t);
  //fprintf (stderr, "   >>>> has time[%d] \n", has_time);

  // if the force is active then deactivate it //

  if (active && !has_time) { 
    if (line_geometry) {
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
    if (line_geometry) {
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

  //fprintf (stderr, "   >>>> active  [%d] \n", active);
  return active;
  }

//*============================================================*
//*==========               setUsePickingSphere      ==========*
//*============================================================*

void 
PmExplicitForce::setUsePickingSphere(const bool flag)
  {
  use_picking_sphere = flag;
  }

void 
PmExplicitForce::setPickingSphereRadius(const float val)
  {
  picking_sphere_radius = val;
  }

//*============================================================*
//*==========               getInteractive           ==========*
//*============================================================*
// get if the force is interactive. 

void 
PmExplicitForce::getInteractive(bool& flag) {
  flag = interactive;
  }

void 
PmExplicitForce::setInteractive(const bool flag) {
  interactive = flag;
  }

void 
PmExplicitForce::setStrength(const float val)
  {
  strength = val;
  }

//*============================================================*
//*==========               setAxes                  ==========*
//*============================================================*

void 
PmExplicitForce::setAxes(PmVector3& axis1, PmVector3& axis2, PmVector3& axis3)
  {
  interactive_axes.push_back(axis1);
  interactive_axes.push_back(axis2);
  interactive_axes.push_back(axis3);
  }

void 
PmExplicitForce::setXform(PmXform& xform) {
  this->xform = xform;
  }

//*============================================================*
//*==========               procAxesSelect           ==========*
//*============================================================*

void
PmExplicitForce::procAxesSelect(PmQuery& query)
  {
  //fprintf (stderr, "\n>>>>>>> PmExplicitForce::procAxesSelect \n");

  if (interactive) {
    PmVector3 axis, taxis;

    if (query.name.find("pick_sphere") != string::npos) {
      axis = query.point - current_point;
      axis.normalize(); 
      taxis = strength*axis;
      }
    else {
      int axis_num = query.entity-1;

      if (axis_num == 0) {
        axis = interactive_axes[0];
        }
      else if (axis_num == 1) {
        axis = interactive_axes[1];
        }
      else if (axis_num == 2) {
        axis = interactive_axes[2];
        }
      else if (axis_num == 3) {
        axis = -interactive_axes[0];
        }
      else if (axis_num == 4) {
        axis = -interactive_axes[1];
        }
      else if (axis_num == 5) {
        axis = -interactive_axes[2];
        }

      if (query.name.find("ends") != string::npos) {
        fprintf (stderr, ">>> select torque \n");
        this->setTorque(true);
        }
      else { 
        this->setTorque(false);
        }

      PmMatrix3x3 mat;
      mat = xform.matrix;
      axis.normalize();
      taxis = strength*(mat*axis);
      }

    fprintf (stderr, ">>> taxs=%f %f %f \n", taxis[0], taxis[1], taxis[2]);
    this->setDirection(taxis);
    this->display(current_point);
    }
  }

//*============================================================*
//*==========               getDirection             ==========*
//*============================================================*
// get the direction of a force. 

void
PmExplicitForce::getDirection(const float time, PmVector3& dir)
  {
  float dr;

  if (this->getRampParams(time, dr)) {
    dir = dr*this->direction;
    }
  else {
    dir = this->direction;
    }
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a force.

void 
PmExplicitForce::display(PmVector3& loc)
  {
  #define ndbg_PmExplicitForce_display
  #ifdef dbg_PmExplicitForce_display
  fprintf (stderr, "\n>>>>>>> PmExplicitForce::display \n");
  fprintf (stderr, "   >>>> loc %g %g %g \n", loc[0], loc[1], loc[2]);
  fprintf (stderr, "   >>>> dir %g %g %g \n", direction[0], direction[1], direction[2]);
  fprintf (stderr, "   >>>> scale %g \n", scale); 
  #endif

  PmGraphicsLine *line, *axes_lines=NULL;
  PmGraphicsPoint *point;
  PmGraphicsSphere *axes_spheres=NULL;
  PmGraphicsSphere *pick_sphere=NULL;
  PmGraphicsAttributes latts, patts, alatts;
  current_point = loc;

  if (!line_geometry) {
    string geom_name;
    geom_name = "force[" + name + "]";
    PmVector3 *verts = new PmVector3[2];
    verts[0] = loc;
    verts[1] = loc + scale*direction;
    line = new PmGraphicsLine(geom_name, 2, verts);

    latts.setColor(color);
    latts.setLineWidth(2.0);
    line->setAttributes(latts);
    line->display();
    line_geometry = dynamic_cast<PmGraphicsGeometry*>(line);

    //patts.setScale(scale);
    patts.setColor(color);
    patts.setScale(0.1);
    patts.setMarker(true);
    point = new PmGraphicsPoint(geom_name, 1, verts);
    point->setAttributes(patts);
    point->display();
    point_geometry = dynamic_cast<PmGraphicsGeometry*>(point);

    if (interactive) {
      string geom_name;
      PmVector3 acolor;
      geom_name = "force[" + name + "]axes";
      PmVector3 *verts = new PmVector3[12];

      for (int i = 0; i < 3; i++) {
        verts[2*i] = loc;
        verts[2*i+1] = loc + interactive_axes[i];
        verts[2*i+6] = loc;
        verts[2*i+7] = loc - interactive_axes[i];
        }

      axes_lines = new PmGraphicsLine(geom_name, 12, verts);
      acolor.set(1,1,1);
      alatts.setDisjoint(true);
      alatts.setLineWidth(1.0);
      alatts.setColor(acolor);
      axes_lines->setAttributes(alatts);
      axes_lines->display();
      axes_line_geometry = dynamic_cast<PmGraphicsGeometry*>(axes_lines);

      geom_name = "force[" + name + "]ends";
      verts = new PmVector3[6];
      float *rads = new float[6];

      for (int i = 0; i < 3; i++) {
        verts[i] = loc + interactive_axes[i];
        verts[i+3] = loc - interactive_axes[i];
        rads[i] = 0.1;
        rads[i+3] = 0.1;
        }

      axes_spheres = new PmGraphicsSphere(geom_name, 6, verts, rads);
      alatts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
      axes_spheres->setAttributes(alatts);
      axes_spheres->display();
      axes_sphere_geometry = dynamic_cast<PmGraphicsGeometry*>(axes_spheres);

      if (use_picking_sphere) { 
        geom_name = "force[" + name + "]pick_sphere";
        verts = new PmVector3[1];
        verts[0] = loc;
        rads = new float[1];
        rads[0] = picking_sphere_radius;
        pick_sphere = new PmGraphicsSphere(geom_name, 1, verts, rads);
        //alatts.setDisplayType(PM_GEOMETRY_DISPLAY_POINT);
        alatts.setDisplayType(PM_GEOMETRY_DISPLAY_LINE);
        acolor.set(0.4,0.4,0.4);
        alatts.setColor(acolor);
        pick_sphere->setAttributes(alatts);
        pick_sphere->display();
        picking_sphere_geometry = dynamic_cast<PmGraphicsGeometry*>(pick_sphere);
        }
      }
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(line_geometry);
    point = dynamic_cast<PmGraphicsPoint*>(point_geometry);
    axes_lines = dynamic_cast<PmGraphicsLine*>(axes_line_geometry);
    axes_spheres = dynamic_cast<PmGraphicsSphere*>(axes_sphere_geometry);
    pick_sphere = dynamic_cast<PmGraphicsSphere*>(picking_sphere_geometry);
    }

  PmVector3 verts[12]; 

  if (axes_lines) {
    PmBody *body;
    //PmXform xform;
    //this->getObj((void**)&body);
    //fprintf (stderr, ">>>> body=%x \n", body); 
    //body->getXform(xform);

    // transform axes //

    PmMatrix3x3 mat;
    PmVector3 axis;
    mat = xform.matrix;

    for (int i = 0; i < 3; i++) {
      axis = mat*interactive_axes[i];
      //axis = interactive_axes[i];
      verts[2*i] = loc;
      verts[2*i+1] = loc + axis;
      verts[2*i+6] = loc;
      verts[2*i+7] = loc - axis;
      }

    axes_lines->update(12, verts);

    for (int i = 0; i < 3; i++) {
      axis = mat*interactive_axes[i];
      //axis = interactive_axes[i];
      verts[i] = loc + axis;
      verts[i+3] = loc - axis;
      }

    axes_spheres->update(6, verts);

    if (use_picking_sphere) { 
      verts[0] = loc;
      pick_sphere->update(1, verts);
      }

    // update force //

    verts[0] = loc;
    verts[1] = loc + scale*direction;
    }
  else {
    verts[0] = loc;
    verts[1] = loc + scale*direction;

    if (has_ramp) {
      verts[1] = loc + dramp*scale*direction;
      }
    else {
      verts[1] = loc + scale*direction;
      }
    }

  line->update(2, verts);
  point->update(1, verts);
  }

//*============================================================*
//*==========               torque                   ==========*
//*============================================================*
// set/get torque flag.

void
PmExplicitForce::setTorque(const bool flag) {
  torque = flag;
  }

bool 
PmExplicitForce::isTorque() { return torque; }

///////////////////////////////////////////////////////////////
//           r a n d o m     f o r c e s                    //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmRandomForce::PmRandomForce(const string name)
  {
  this->name = name;
  scale = 1.0;
  sobj = NULL;
  line_geometry = NULL;
  point_geometry = NULL;
  active = false;
  //global_frame = false;
  standard_deviation = 1.0;
  mean = 0.0;
  seed = 5319321;
  current_seed = seed;
  }

//*============================================================*
//*==========               isActive                 ==========*
//*============================================================*
// check if a force is active for the given time.

bool
PmRandomForce::isActive(float t)
  {
  bool has_time = time_interval.hasTime(t);
  //fprintf (stderr, "\n>>>>>>> PmExplicitForce::isActive  t[%f] \n", t);
  //fprintf (stderr, "   >>>> has time[%d] \n", has_time);

  // if the force is active then deactivate it //

  if (active && !has_time) {
    if (line_geometry) {
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
    if (line_geometry) {
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

  //fprintf (stderr, "   >>>> active  [%d] \n", active);
  return active;
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a force.

void
PmRandomForce::display(PmVector3& loc) { 
  }

//*============================================================*
//*==========               getData                  ==========*
//*============================================================*
// get the data generating a random force.

void
PmRandomForce::getData(float& sd, float& mean, long int& seed)
  {
  sd = this->standard_deviation;
  mean = this->mean;
  seed = this->current_seed;
  }

//*============================================================*
//*==========               setMean                  ==========*
//*============================================================*
// set the mean for generating a random force.

void
PmRandomForce::setMean(const float mean)
  {
  this->mean = mean;
  }

//*============================================================*
//*==========               setStandardDev           ==========*
//*============================================================*
// set the standard deviation for generating a random force.

void
PmRandomForce::setStandardDev(const float sd)
  {
  this->standard_deviation = sd;
  }

//*============================================================*
//*==========               setSeed                  ==========*
//*============================================================*
// set the seed for generating a random force.

void
PmRandomForce::setSeed(const long int seed)
  {
  this->current_seed = seed;
  }

//*============================================================*
//*==========               genForce                 ==========*
//*============================================================*
// generate a random force.                       

void
PmRandomForce::genForce(PmVector3& rvec)
  {
  pm_MathRandGaussVector (mean, standard_deviation, rvec, &current_seed);
  /*
  fprintf (stderr, "   >>>> random vec  %f %f %f   %g  \n", rvec[0], rvec[1], rvec[2],
           rvec.length());
  */
  }

}





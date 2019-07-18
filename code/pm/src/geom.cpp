
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
//* geom:                 g e o m e t r y                      *
//*============================================================*

#include "geom.h"
#include "graphics.h"

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGeometry::PmGeometry() {
  init();
  }

PmGeometry::PmGeometry(const string name) {
  this->name = name;
  init();
  }

void
PmGeometry::init()
  {
  body = NULL;
  sobj = NULL;
  global_frame = false;
  display_type = PM_GEOMETRY_DISPLAY_LINE;
  line_geometry = NULL;
  point_geometry = NULL;
  color.set(1,1,1);
  width = 1.0;
  }

//*============================================================*
//*==========               set/getBody              ==========*
//*============================================================*

void
PmGeometry::getBody(PmBody **body) {
  *body = this->body;
  }

void
PmGeometry::setBody(PmBody *body) {
  this->body = body;
  }

//*============================================================*
//*==========               setDisplayType           ==========*
//*============================================================*
// set display type.

void 
PmGeometry::setDisplayType(const PmGeometryDisplayType dtype) {
  display_type = dtype;
  }

///////////////////////////////////////////////////////////////
//                    v e c t o r                           //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGeometryVector::PmGeometryVector(const string name, PmVector3 origin, PmVector3 dir)
  {
  this->name = name;
  this->origin = origin;
  this->direction = dir;
  this->scale = 1.0;

  points = new PmVector3[2];
  points[0] = origin; 
  points[1] = origin + direction; 

  curr_points = new PmVector3[2];
  curr_points[0] = points[0];
  curr_points[1] = points[1];
  }

//*============================================================*
//*==========               setScale                 ==========*
//*============================================================*
// set vector scale.

void 
PmGeometryVector::setScale(const float scale) 
  {
  this->scale = scale;
  points[0] = this->origin; 
  points[1] = this->origin + scale*this->direction; 
  }

//*============================================================*
//*==========               update                   ==========*
//*============================================================*
// update  vector.

void 
PmGeometryVector::update(PmXform& xform)
  {
  PmMatrix3x3 mat = xform.matrix;

  for (int i = 0; i < 2; i++) {
    curr_points[i] = mat*(points[i]-xform.center) + xform.translation + xform.center;
    }
  }

//*============================================================*
//*==========               getCurrPoints            ==========*
//*============================================================*
// get the current (transformed) points for a vector.

void 
PmGeometryVector::getCurrPoints(vector<PmVector3>& pts)
  {
  pts.push_back(curr_points[0]);
  pts.push_back(curr_points[1]);
  }

//*============================================================*
//*==========               getCurrVector            ==========*
//*============================================================*
// get the current (transformed) vector.

void 
PmGeometryVector::getCurrVector(PmVector3& v) 
  {
  v = curr_points[1] - curr_points[0];
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a vector.

void 
PmGeometryVector::display(const bool show)
  {
  PmGraphicsLine *line;
  PmGraphicsPoint *mpoint;
  PmGraphicsAttributes atts, patts;
  //fprintf (stderr, "\n>>>>>> PmGeometryVector::display  name = %s \n", name.c_str());
  //fprintf (stderr, ">>> line_geometry = %x \n", line_geometry); 

  if (!line_geometry) {
    string geom_name;
    geom_name = "vector[" + name + "]line";
    line = new PmGraphicsLine(geom_name, 2, points, false);
    atts.setColor(color);
    atts.setLineWidth(width);
    atts.setDisplayType(display_type);
    line->setAttributes(atts);
    line->display();
    line_geometry = dynamic_cast<PmGraphicsGeometry*>(line);

    geom_name = "vector[" + name + "]marker";
    patts.setScale(0.1*width);
    patts.setMarker(true);
    patts.setColor(color);
    mpoint = new PmGraphicsPoint(geom_name, 1, points);
    mpoint->setAttributes(patts);
    mpoint->display();
    point_geometry = dynamic_cast<PmGraphicsGeometry*>(mpoint);
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(line_geometry);
    mpoint = dynamic_cast<PmGraphicsPoint*>(point_geometry);
    }

  /*
  fprintf (stderr, ">>> curr points = %f %f %f \n", curr_points[0][0], 
           curr_points[0][1], curr_points[0][2]);
  fprintf (stderr, ">>> curr points = %f %f %f \n", curr_points[1][0], 
           curr_points[1][1], curr_points[1][2]);
  */

  line->update(2, curr_points);
  mpoint->update(1, curr_points);
  }

///////////////////////////////////////////////////////////////
//                    c u r v e                             //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGeometryCurve::PmGeometryCurve(const string name, const vector<PmVector3>& pts) 
  {
  fprintf (stderr, "\n>>>>>> PmGeometryCurve::ctor name=%s\n", name.c_str());
  this->name = name;
  number_of_points = pts.size();
  points = new PmVector3[number_of_points];
  curr_points = new PmVector3[number_of_points];

  for (int i = 0; i < number_of_points; i++) {
    points[i] = pts[i];
    curr_points[i] = pts[i];
    }
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a curve. 

void
PmGeometryCurve::display(const bool show)
  {
  fprintf (stderr, "\n>>>>>> PmGeometryCurve::display name=%s\n", name.c_str());
  PmGraphicsLine *line;
  PmGraphicsSphere *sphere;
  PmGraphicsPoint *mpoint;
  PmGraphicsAttributes atts, patts;

  if (!line_geometry) {
    string geom_name;
    geom_name = "curve[" + name + "]line";
    line = new PmGraphicsLine(geom_name, number_of_points, points, false);
    atts.setColor(color);
    atts.setLineWidth(width);
    atts.setDisplayType(display_type);
    line->setAttributes(atts);
    line->display();
    line_geometry = dynamic_cast<PmGraphicsGeometry*>(line);
   
    // points //
    geom_name = "curve[" + name + "]points";
    float *rads = new float[number_of_points];

    for (int i = 0; i < number_of_points; i++) {
      rads[i] = 0.05;
      }

    sphere = new PmGraphicsSphere(geom_name, number_of_points, points, rads);
    sphere->display();
    sphere_geometry = dynamic_cast<PmGraphicsGeometry*>(sphere);
    }

  else {
    line = dynamic_cast<PmGraphicsLine*>(line_geometry);
    sphere = dynamic_cast<PmGraphicsSphere*>(sphere_geometry);
    }

  line->update(number_of_points, curr_points);
  sphere->update(number_of_points, curr_points);
  }

//*============================================================*
//*==========               update                   ==========*
//*============================================================*

void
PmGeometryCurve::update(PmXform& xform)
  {
  fprintf (stderr, "\n>>>>>> PmGeometryCurve::update name=%s \n", name.c_str());
  PmMatrix3x3 mat = xform.matrix;

  for (int i = 0; i < number_of_points; i++) {
    curr_points[i] = mat*(points[i]-xform.center) + xform.translation + xform.center;
    }
  }

//*============================================================*
//*==========               getCurrPoints            ==========*
//*============================================================*
// get the current (transformed) points for a curve.

void
PmGeometryCurve::getCurrPoints(vector<PmVector3>& pts)
  {
  for (int i = 0; i < number_of_points; i++) {
    pts.push_back(curr_points[i]); 
    }
  }

}



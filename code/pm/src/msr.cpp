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
//* msr:                   m e a s u r e m e n t               *
//*============================================================*

#include "pm/mth.h"
#include "msr.h"
#include "geom.h"
#include "graphics.h"

namespace ProteinMechanica {


//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmMeasurement::PmMeasurement()
  {
  color.set(1,1,1);
  file_ptr = NULL;
  binary_format = true;
  write_init = true;
  last_written = 0;
  print = false;
  }

//*============================================================*
//*==========         convMeasurementType            ==========*
//*============================================================*

void
PmMeasurement::convMeasurementType(const string str, PmMeasurementType& type)
  {
  if (str == "angle") {
    type = PM_MEASUREMENT_ANGLE;
    }
  else if (str == "distance") {
    type = PM_MEASUREMENT_DISTANCE;
    }
  else if (str == "position") {
    type = PM_MEASUREMENT_POSITION;
    }
  else {
    type = PM_MEASUREMENT_UNKNOWN;
    }
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create a measurement object of a particular type.

PmMeasurement*
PmMeasurement::create(const string name, const PmMeasurementType type)
  {
  switch (type) {
    case PM_MEASUREMENT_ANGLE:
      {
      PmAngleMeasurement *amsr = new PmAngleMeasurement(name);
      amsr->type = type;
      return amsr;
      }
    break;

    case PM_MEASUREMENT_DISTANCE:
      {
      PmDistMeasurement *dmsr = new PmDistMeasurement(name);
      dmsr->type = type;
      return dmsr;
      }
    break;

    case PM_MEASUREMENT_POSITION:
      {
      PmPositionMeasurement *pmsr = new PmPositionMeasurement(name);
      pmsr->type = type;
      return pmsr;
      }
    break;

    default:
      return NULL;
    break;
    }
  }

//*============================================================*
//*==========               setColor                 ==========*
//*============================================================*
// set color.                   

void 
PmMeasurement::setColor(const PmVector3& color) {
  this->color = color;
  }

//*============================================================*
//*==========               getGlobalFrames          ==========*
//*============================================================*
// get global frame flags.

void 
PmMeasurement::getGlobalFrames(vector<bool>& frames) {
  frames = global_frames;
  }

void 
PmMeasurement::setGlobalFrames(const vector<bool> frames) {
  global_frames = frames;
  }

//*============================================================*
//*==========               getObjs                  ==========*
//*============================================================*
// get physical objects.                      

void 
PmMeasurement::getObjs(vector<void*>& objs) {
  objs = sobjs;
  }

void 
PmMeasurement::setObjs(const vector<void*>objs) {
  sobjs = objs;
  }

//*============================================================*
//*==========               getPoints                ==========*
//*============================================================*
// get points.              

void 
PmMeasurement::getPoints(vector<PmVector3>& pts) {
  pts = points;
  }

void 
PmMeasurement::setPoints(const vector<PmVector3> pts) {
  points = pts;
  current_points = pts;
  }

//*============================================================*
//*==========             setPrintFlag               ==========*
//*============================================================*
// get/set print flag.                            

void
PmMeasurement::getPrintFlag(bool& flag) {
  flag = print;
  }

void 
PmMeasurement::setPrintFlag(const bool flag) {
  print = flag;
  }

//*============================================================*
//*==========       setWriteResultsParams            ==========*
//*============================================================*
// set parameterts for writing results.

void 
PmMeasurement::setWriteResultsParams(const string file_name, const bool binary)
  {
  this->write_file_name = file_name;
  }

///////////////////////////////////////////////////////////////
//                d i s t a n c e                           //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmDistMeasurement::PmDistMeasurement(const string name) 
  {
  this->name = name;
  geometry = NULL;
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a distance measurement. 

void 
PmDistMeasurement::display()
  {
  PmGraphicsLine *line;
  PmGraphicsAttributes atts;

  if (!geometry) {
    string geom_name;
    geom_name = "distance[" + name + "]";
    PmVector3 *verts = new PmVector3[2];
    verts[0] = points[0];
    verts[1] = points[1];
    line = new PmGraphicsLine(geom_name, 2, verts);
    atts.setColor(color);
    line->setAttributes(atts);
    line->display();
    geometry = dynamic_cast<PmGraphicsGeometry*>(line);
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(geometry);
    }

  PmVector3 verts[2]; 
  verts[0] = current_points[0];
  verts[1] = current_points[1];
  line->update(2, verts);
  }

//*============================================================*
//*==========               updatePoints             ==========*
//*============================================================*
// update points.

void
PmDistMeasurement::updatePoints(const vector<PmVector3> pts) 
  {
  //fprintf (stderr, "\n>>>>>> PmDistMeasurement::updatePoints \n");
  PmVector3 v;

  if (distances.empty()) {
    v = points[1] - points[0];
    distances.push_back(v.length());
    }

  current_points = pts;
  v = current_points[1] - current_points[0];
  distances.push_back(v.length());
  }

//*============================================================*
//*==========               writeResults             ==========*
//*============================================================*
// write distance measurement results.

void 
PmDistMeasurement::writeResults(const vector<float>& time, const int lastw, 
                                const int num)
  {
  char fname[100];

  if (this->write_file_name.empty()) {
    return;
    }

  // write file header //

  if (this->write_init) {
    sprintf (fname, "%s.pm", this->write_file_name.c_str());
    this->file_ptr = fopen (fname, "w");

    if (!this->file_ptr) {
      return;
      }

    fprintf (this->file_ptr, "# Protein Mechanica distance measurement file \n");
    fprintf (this->file_ptr, "# time distance           \n");
    fprintf (this->file_ptr, "distance measurement \"%s\" \n", name.c_str());
    this->write_init = false;
    }

  //===== write data =====//

  if (distances.empty()) {
    return;
    }

  for (int n = lastw; n < num; n++) {
    fprintf (this->file_ptr, "%f %f \n", time[n], distances[n]);

    if (print) {
      fprintf (stderr, ">>> time=%f distance=%f \n", time[n], distances[n]);
      }
    }
  }

///////////////////////////////////////////////////////////////
//                      a n g l e                           //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmAngleMeasurement::PmAngleMeasurement(const string name)
  {
  this->name = name;
  geometry = NULL;
  vector_geom1 = NULL; 
  vector_geom2 = NULL;
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display an angle measurement. 

void
PmAngleMeasurement::display() 
  { 

  if (vector_geom1) {
    return;
    }

  PmGraphicsLine *line;
  PmGraphicsAttributes atts;

  if (!geometry) {
    string geom_name;
    geom_name = "angle[" + name + "]";
    PmVector3 *verts = new PmVector3[3];
    verts[0] = points[0];
    verts[1] = points[1];
    verts[2] = points[2];
    atts.setColor(color);
    line = new PmGraphicsLine(geom_name, 3, verts);
    line->setAttributes(atts);
    line->display();
    geometry = dynamic_cast<PmGraphicsGeometry*>(line);
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(geometry);
    }

  PmVector3 verts[3];
  verts[0] = current_points[0];
  verts[1] = current_points[1];
  verts[2] = current_points[2];
  line->update(3, verts);
  }

//*============================================================*
//*==========               updatePoints             ==========*
//*============================================================*
// update points.

void
PmAngleMeasurement::updatePoints(const vector<PmVector3> pts) 
  {
  //fprintf (stderr, "\n>>>>>> updatePoints \n");
  PmVector3 v1, v2, cp;
  float dp, angle;

  if (vector_geom1) {
    vector_geom1->getCurrVector(v1);
    vector_geom2->getCurrVector(v2);
    }

  else { 
    if (angles.empty()) {
      v1 = points[1] - points[0];
      v2 = points[2] - points[0];
      }
    else {
      v1 = pts[1] - pts[0];
      v2 = pts[2] - pts[0];
      }
    }

  v1.normalize();
  v2.normalize();
  dp = v1*v2;
  angle = 180.0*acos(dp) / M_PI;
  angles.push_back(angle);
  /*
  fprintf (stderr, ">>> v1 = %f %f %f \n", v1[0], v1[1], v1[2]);
  fprintf (stderr, ">>> v2 = %f %f %f \n", v2[0], v2[1], v2[2]);
  fprintf (stderr, ">>> dp = %f \n", dp);
  */

  //cp = v1.cross(v2);
  //angles.push_back(asin(cp.length()));

  current_points = pts;
  }

//*============================================================*
//*==========               setVectors               ==========*
//*============================================================*
// set the vectors for an angle measurement.

void 
PmAngleMeasurement::setVectors(PmGeometryVector *vgeom1, PmBody *body1, 
                               PmGeometryVector *vgeom2, PmBody *body2)
  {
  vector_geom1 = vgeom1; 
  vector_body1 = body1;
  vector_geom2 = vgeom2;
  vector_body2 = body2;
  }

//*============================================================*
//*==========               writeResults             ==========*
//*============================================================*
// write angle measurment results.

void
PmAngleMeasurement::writeResults(const vector<float>& time, const int lastw,
                                 const int num)
  {
  char fname[100];

  if (this->write_file_name.empty()) {
    return;
    }

  // write file header //

  if (this->write_init) {
    sprintf (fname, "%s.pm", this->write_file_name.c_str());
    this->file_ptr = fopen (fname, "w");

    if (!this->file_ptr) {
      return;
      }

    fprintf (this->file_ptr, "# Protein Mechanica angle measurement file \n");
    fprintf (this->file_ptr, "# time angle              \n");
    fprintf (this->file_ptr, "angle measurement \"%s\" \n", name.c_str());
    this->write_init = false;
    }

  // write data //

  //fprintf (stderr, "##### angles.size = %d \n", angles.size());
  //fprintf (stderr, "##### time  .size = %d \n", time.size());

  if (print) {
    fprintf (stderr, "angle measurement \"%s\" \n", name.c_str());
    }

  if (angles.size() != 0) {
    for (int n = lastw; n < num; n++) {
      fprintf (this->file_ptr, "%f %f \n", time[n], angles[n]);

      if (print) {
        fprintf (stderr, "%f %f \n", time[n], angles[n]);
        }
      }
    }
  }

///////////////////////////////////////////////////////////////
//                p o s i t i o n                           //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmPositionMeasurement::PmPositionMeasurement(const string name)
  {
  this->name = name;
  geometry = NULL;
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a position measurement. 

void
PmPositionMeasurement::display() 
  { 
  PmGraphicsPoint *pt;
  PmGraphicsAttributes atts;

  if (!geometry) {
    string geom_name;
    geom_name = "position[" + name + "]";
    PmVector3 *verts = new PmVector3[1];
    verts[0] = points[0];
    atts.setColor(color);
    atts.setScale(0.1);
    pt = new PmGraphicsPoint(geom_name, 1, verts);
    pt->setAttributes(atts);
    pt->display();
    geometry = dynamic_cast<PmGraphicsGeometry*>(pt);
    }
  else {
    pt = dynamic_cast<PmGraphicsPoint*>(geometry);
    }

  PmVector3 verts[1];
  verts[0] = current_points[0];
  pt->update(1, verts);
  }

//*============================================================*
//*==========               updatePoints             ==========*
//*============================================================*
// update points.

void
PmPositionMeasurement::updatePoints(const vector<PmVector3> pts) 
  {
  if (positions.empty()) {
    positions.push_back(points[0]);
    }
  else {
    positions.push_back(pts[0]);
    }

  current_points = pts;
  }

//*============================================================*
//*==========               writeResults             ==========*
//*============================================================*
// write position measurment results.

void
PmPositionMeasurement::writeResults(const vector<float>& time, const int lastw,
                                    const int num)
  {
  char fname[100];

  if (this->write_file_name.empty()) {
    return;
    }

  // write file header //

  if (this->write_init) {
    sprintf (fname, "%s.pm", this->write_file_name.c_str());
    this->file_ptr = fopen (fname, "w");

    if (!this->file_ptr) {
      return;
      }

    fprintf (this->file_ptr, "# Protein Mechanica position measurement file \n");
    fprintf (this->file_ptr, "# time position           \n");
    fprintf (this->file_ptr, "position measurement \"%s\" \n", name.c_str());
    this->write_init = false;
    }

  // write data //

  for (int n = lastw; n < num; n++) {
    fprintf (this->file_ptr, "%f %f %f %f \n", time[n], positions[n][0], 
             positions[n][1], positions[n][2]);
    }
  }

}






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
//* trace:                    t r a c e                        *
//*============================================================*

#include "trace.h"
#include "graphics.h"

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmTrace::PmTrace(const string name) {
  this->name = name;
  init();
  }

PmTrace::PmTrace(const string name, PmBody *body, PmVector3& point) {
  init();
  this->name = name;
  this->sobj = body;
  this->point = point;
  }

void
PmTrace::init()
  {
  sobj = NULL;
  global_frame = false;

  num_points = 0;
  size_points = 100;
  points = (PmVector3*)malloc(sizeof(PmVector3)*size_points);
  time = (float*)malloc(sizeof(float)*size_points);

  file_ptr = NULL;
  binary_format = false;
  write_init = true;
  last_written = 0;

  display_type = PM_GEOMETRY_DISPLAY_LINE;
  line_geometry = NULL;
  point_geometry = NULL;
  color.set(1,1,1);
  line_width = 1.0;
  }

//*============================================================*
//*==========               findInterval             ==========*
//*============================================================*
// find an interval within the trace of a given tolerance.

void 
PmTrace::findInterval(const float tol, const bool show, int& ipt, float& t, 
                      float& min_dl, float& max_dl)
  {
  float dl, dt;
  PmVector3 dv;
  ipt = -1;
  t = -1.0;

  max_dl = 0.0;
  min_dl = 1e6;

  for (int i = 0; i < num_points-1; i++) {
    dv = points[i+1] - points[i];
    dl = dv.length();

    if (dl < min_dl) {
      min_dl = dl;
      }
    else if (dl > max_dl) {
      max_dl = dl;
      }

    if (dl <= tol) {
      ipt = i+1;
      break;
      }
    }

  if (!show || (ipt == -1)) {
    return;
    }

  PmGraphicsPoint *point;
  PmGraphicsAttributes atts;

  if (!point_geometry) {
    string geom_name;
    geom_name = "trace[" + name + "]point";
    point = new PmGraphicsPoint(geom_name, 1, &points[ipt-1]);
    atts.setColor(color);
    atts.setScale(0.1);
    point->setAttributes(atts);
    point->display();
    }
  else {
    point = dynamic_cast<PmGraphicsPoint*>(point_geometry);
    }

  point->update(1, &points[ipt-1]);
  }

//*============================================================*
//*==========               compLength               ==========*
//*============================================================*
// compute the length of a trace.

float
PmTrace::compLength(PmVector3& dir, int ipt)
  {
  float length = 0.0;
  PmVector3 dl;
  float dp;

  if ((ipt > 0) && (ipt <=  num_points)) {
    ipt = ipt - 1;
    }
  else {
    ipt = num_points;
    }

  // project end points of trace onto direction //

  if (dir.length() != 0.0) {
    dl = points[ipt-1] - points[0];
    dp = dl*dir;
    length = dp;
    }

  // compute total length // 

  else {
    for (int i = 0; i < ipt-1; i++) {
      dl = points[i+1] - points[i];
      length += dl.length();
      }
    }

  return length;
  }

//*============================================================*
//*==========               addPoint                 ==========*
//*============================================================*
// add a point along a trace.

void
PmTrace::addPoint(PmVector3& pt)
  {
  //fprintf (stderr, "\n>>>>>>> PmTrace::add pt %f %f %f \n", pt[0], pt[1], pt[2]);
  points[num_points] = pt;
  num_points += 1;

  if (num_points == size_points) {
    size_points += 100;
    points = (PmVector3*)realloc(points, sizeof(PmVector3)*size_points);
    }
  }

void
PmTrace::addPoint(float t, PmVector3& pt)
  {
  points[num_points] = pt;
  time[num_points] = t;
  num_points += 1;

  if (num_points == size_points) {
    size_points += 100;
    points = (PmVector3*)realloc(points, sizeof(PmVector3)*size_points);
    time = (float*)realloc(time, sizeof(float)*size_points);
    }
  }

//*============================================================*
//*==========          setWriteResultsParams         ==========*
//*============================================================*
// set writting results.        

void
PmTrace::setWriteResultsParams(const string file_name, const bool binary)
  {
  //fprintf (stderr, ">>>>>> PmTrace::setWriteResultsParams \n"); 
  this->write_file_name = file_name;
  this->write_init = true;
  this->last_written = 0;
  }

//*============================================================*
//*==========               writeResults             ==========*
//*============================================================*
// write results.

void
PmTrace::writeResults(const vector<float>& time, const int lastw, const int num)
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

    fprintf (this->file_ptr, "# Protein Mechanica trace file \n");
    fprintf (this->file_ptr, "# time position           \n");
    fprintf (this->file_ptr, "trace \"%s\" \n", name.c_str()); 
    this->write_init = false;
    }

  // write data //

  for (int n = lastw; n < num; n++) {
    fprintf (this->file_ptr, "%f  %f %f %f \n", time[n], points[n][0], points[n][1], 
             points[n][2]);
    }

  fflush (this->file_ptr);
  }

//*============================================================*
//*==========               setDisplayType           ==========*
//*============================================================*
// set display type.

void 
PmTrace::setDisplayType(const PmGeometryDisplayType dtype) {
  display_type = dtype;
  }

void PmTrace::setLineWidth(float width)
  {
  line_width = width;
  }

//*============================================================*
//*==========               display                  ==========*
//*============================================================*
// display a trace.

void 
PmTrace::display(const bool show)
  {
  #ifdef dbg_PmTrace_display
  fprintf (stderr, "\n>>>>>>> PmTrace::display  num_points = %d \n", num_points);
  #endif

  if (!num_points) {
    return;
    }

  PmGraphicsLine *line;
  PmGraphicsAttributes atts;

  if (!line_geometry) {
    string geom_name;
    geom_name = "trace[" + name + "]line";
    line = new PmGraphicsLine(geom_name, num_points, points, false);
    atts.setColor(color);
    atts.setDisplayType(display_type);
    atts.setLineWidth(line_width);
    line->setAttributes(atts);
    line->display();
    line_geometry = dynamic_cast<PmGraphicsGeometry*>(line);
    }
  else {
    line = dynamic_cast<PmGraphicsLine*>(line_geometry);
    }

  line->update(num_points, points);
  }

}



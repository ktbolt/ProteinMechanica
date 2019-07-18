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

#ifndef _TRACE_PM_H_
#define _TRACE_PM_H_

#include "pm/pm.h"
#include "graphics.h"

namespace ProteinMechanica {

class PmGraphicsGeometry;


// PmTrace
// -------

class PM_EXPORT PmTrace {

   public:
      PmTrace(const string name);
      PmTrace(const string name, PmBody *body, PmVector3& point);
      ~PmTrace(){};
      void setColor (PmVector3& color) { this->color = color; }
      void getExtent (PmExtent& extent ) { extent = this->extent; }
      void setExtent (PmExtent& extent ) { this->extent = extent; }
      void setGlobalFrame(const bool flag) { this->global_frame = flag; }
      void getGlobalFrame(bool& flag) { flag = this->global_frame; }
      void findInterval(const float tol, const bool show, int& ipt, float& t, 
                        float& min_dl, float& max_dl);
      float compLength(PmVector3& dir, int ipt);
      void getName(string& name) { name = this->name; }
      void getObj(void **obj) { *obj = this->sobj; }
      void setObj(void *obj) { this->sobj = obj; }
      void addPoint(PmVector3& point);
      void addPoint(float t, PmVector3& pt);
      void getPoint(PmVector3& point) { point = this->point; }
      void setPoint(const PmVector3& point) { this->point = point; }
      void setDisplayType(const PmGeometryDisplayType dtype);
      void setLineWidth(float width);
      void setWriteResultsParams(const string file_name, const bool binary);
      void writeResults(const vector<float>& time, const int last_written, const int num);
      void display(const bool show);

   private:
      string name;
      PmVector3 point;
      bool global_frame;
      void *sobj;
      PmExtent extent;

      // data //
      int num_points, size_points;
      PmVector3 *points;
      float *time;

      // graphics object //
      PmVector3 color;
      PmGeometryDisplayType display_type;
      PmGraphicsGeometry *line_geometry; 
      PmGraphicsGeometry *point_geometry; 
      float line_width;

      // results i/o //
      string write_file_name;
      FILE *file_ptr; 
      bool binary_format;
      bool write_init;
      int last_written;

      void init();
   };

}

#endif



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
//* geom:                    g e o m                           *
//*============================================================*

#ifndef _GEOMETRY_PM_H_
#define _GEOMETRY_PM_H_

#include "pm/pm.h"
#include "graphics.h"

namespace ProteinMechanica {

class PmGraphicsGeometry;

//*============================================================*
//*==========              PmGeometry                ==========*
//*============================================================*

class PM_EXPORT PmGeometry {
   public:
      PmGeometry();
      PmGeometry(const string name);
      ~PmGeometry(){};
      void getBody(PmBody **body);
      void setBody(PmBody *body);
      void setColor (PmVector3& color) { this->color = color; }
      void getExtent (PmExtent& extent ) { extent = this->extent; }
      void setGlobalFrame(const bool flag) { this->global_frame = flag; }
      void getGlobalFrame(bool& flag) { flag = this->global_frame; }
      void getName(string& name) { name = this->name; }
      void getObj(void **obj) { *obj = this->sobj; }
      void setObj(void *obj) { this->sobj = obj; }
      void setDisplayType(const PmGeometryDisplayType dtype);
      void setWidth(const float width) { this->width = width; }
      virtual void display(const bool show)=0;
      virtual void update(PmXform& xform)=0;
      virtual void getCurrPoints(vector<PmVector3>& pts)=0;

   protected:
      string name;
      bool global_frame;
      void *sobj;
      PmBody *body;
      PmExtent extent;

      // graphics object //
      float width;
      PmVector3 color;
      PmGeometryDisplayType display_type;
      PmGraphicsGeometry *line_geometry; 
      PmGraphicsGeometry *point_geometry; 
      void init();
   };

//*============================================================*
//*==========         PmGeometryVector               ==========*
//*============================================================*

class PM_EXPORT PmGeometryVector: public PmGeometry {
   public:
      PmGeometryVector(const string name, PmVector3 origin, PmVector3 dir) ;
      ~PmGeometryVector(){};
      void display(const bool show);
      void setScale(const float scale);
      void update(PmXform& xform);
      void getCurrPoints(vector<PmVector3>& pts);
      void getCurrVector(PmVector3& v);

   private:
      PmVector3 origin; 
      PmVector3 direction;
      PmVector3 *points, *curr_points; 
      float scale;
   };

//*============================================================*
//*==========         PmGeometryCurve                ==========*
//*============================================================*

class PM_EXPORT PmGeometryCurve: public PmGeometry {
   public:
      PmGeometryCurve(const string name, const vector<PmVector3>& pts) ;
      ~PmGeometryCurve(){};
      void display(const bool show);
      void update(PmXform& xform);
      void getCurrPoints(vector<PmVector3>& pts);

   private:
      int number_of_points;
      PmVector3 *points, *curr_points;
      PmGraphicsGeometry *sphere_geometry; 
   };

}

#endif


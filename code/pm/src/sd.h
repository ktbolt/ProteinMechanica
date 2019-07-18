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
//* sd:         s p a t i a l  d e c o m p o s i t i o n       *
//*============================================================*

#ifndef _SPATIAL_DECOMP_PM_H_ 
#define _SPATIAL_DECOMP_PM_H_

#include "pm/pm.h"
#include "graphics.h"

namespace ProteinMechanica {

typedef struct PmSpatialCell {
  vector<int> ids;
  } PmSpatialCell;


// PmSpatialDecomp
// ---------------

class PM_EXPORT PmSpatialDecomp {
   public:
      PmSpatialDecomp(const string name);
      void setCoordSys(PmVector3& origin, PmVector3& u, PmVector3& v, PmVector3& w);
      void setCellSize(const float du, const float dv, const float dw);
      void setAxesLength(const float lu, const float lv, const float lw);
      void addPoint(const int id, const PmVector3& pt, const float r);
      void getPoints(const PmVector3& pt, const float r, vector<int>& ids);
      void mapPoint(const PmVector3& pt, const float r, int& imin, int& imax,
                    int& jmin, int& jmax, int& kmin, int& kmax);
      void printInfo();
      void setXform(PmXform& xform);
      void display();

   private:
      string name;
      bool coord_sys_set;
      PmVector3 origin, u_axis, v_axis, w_axis;
      PmVector3 curr_origin, curr_u_axis, curr_v_axis, curr_w_axis;
      float u_length, v_length, w_length; 
      float cell_du, cell_dv, cell_dw;
      int num_cells, num_ucell, num_vcell, num_wcell;
      PmSpatialCell *cells;
      PmXform xform;
      PmGraphicsLine *line, *line1;
      void initialize();
      void mapLocalCoord(const float c, const float d, const int nc, int& i);
   };

   
}

#endif




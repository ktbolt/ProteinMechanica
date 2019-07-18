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

#include "sd.h"
#include "graphics.h"
#include "pm/mth.h"

namespace ProteinMechanica {

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmSpatialDecomp::PmSpatialDecomp(const string str) 
  {
  name = str;
  cell_du = 0.0;
  cell_dv = 0.0; 
  cell_dw = 0.0;
  num_cells = 0; 
  num_ucell = 0; 
  num_vcell = 0; 
  num_wcell = 0;
  coord_sys_set = false;
  u_length = 0.0;
  v_length = 0.0;
  w_length = 0.0;
  line = NULL;
  line1 = NULL;
  }

//*============================================================*
//*==========             setCoordSys                ==========*
//*============================================================*
// set the coordinate system.

void 
PmSpatialDecomp::setCoordSys(PmVector3& origin, PmVector3& u, PmVector3& v, PmVector3& w)
  {
  this->origin = origin;
  this->u_axis = u;
  this->v_axis = v;
  this->w_axis = w;

  this->curr_origin = origin;
  this->curr_u_axis = u;
  this->curr_v_axis = v;
  this->curr_w_axis = w;

  this->coord_sys_set = true;
  this->cells = NULL;
  }

//*============================================================*
//*==========             setCellSize                ==========*
//*============================================================*
// set the size of the cells.

void
PmSpatialDecomp::setCellSize(const float du, const float dv, const float dw) 
  {
  cell_du = du;
  cell_dv = dv; 
  cell_dw = dw;
  }

//*============================================================*
//*==========             setAxesLength              ==========*
//*============================================================*
// set axes lengths.           

void 
PmSpatialDecomp::setAxesLength(const float lu, const float lv, const float lw)
  {
  u_length = lu;
  v_length = lv;
  w_length = lw;
  }

//*============================================================*
//*==========                mapLocalCoord           ==========*
//*============================================================*
// map a local coordinate into a cell index.

void
PmSpatialDecomp::mapLocalCoord(const float c, const float d, const int nc, int& i)
  {
  i = static_cast<int>(c / d);
  if (i < 0) i = 0;
  if (i > nc-1) i = nc-1;
  }

//*============================================================*
//*==========                addPoint                ==========*
//*============================================================*
// add a point.     

void
PmSpatialDecomp::addPoint(const int id, const PmVector3& pt, const float r)
  {
  PmSpatialCell *cell;
  int imin, imax, jmin, jmax, kmin, kmax;

  if (!cells) {
    initialize();
    }

  mapPoint(pt, r, imin, imax, jmin, jmax, kmin, kmax);

  for (int k = kmin; k <= kmax; k++) {
    for (int j = jmin; j <= jmax; j++) {
      for (int i = imin; i <= imax; i++) {
        cell = &cells[i + num_ucell*j + num_ucell*num_vcell*k];
        cell->ids.push_back(id);
        }
      }
    }
  }

//*============================================================*
//*==========                getPoints               ==========*
//*============================================================*
// get points at the given point within radius r.

void 
PmSpatialDecomp::getPoints(const PmVector3& pt, const float r, vector<int>& ids)
  {
  //fprintf (stderr, "\n>>>>>> PmSpatialDecomp::getPoints \n");
  PmVector3 p;
  int imin, imax, jmin, jmax, kmin, kmax;
  int id;
  bool found;
  PmSpatialCell *cell;

  ids.clear();

  // get indexes into cells //

  mapPoint(pt, r, imin, imax, jmin, jmax, kmin, kmax);

  /*
  fprintf (stderr, ">>> imin imax  %d %d \n", imin, imax);
  fprintf (stderr, ">>> jmin jmax  %d %d \n", jmin, jmax);
  fprintf (stderr, ">>> kmin kmax  %d %d \n", kmin, kmax);
  fprintf (stderr, ">>> indexes  %d %d   %d %d  %d %d \n", imin, imax, jmin, jmax,
           kmin, kmax);
  */

  for (int k = kmin; k <= kmax; k++) {
    for (int j = jmin; j <= jmax; j++) {
      for (int i = imin; i <= imax; i++) {
        cell = &cells[i + num_ucell*j + num_ucell*num_vcell*k]; 

        if (cell->ids.size() != 0) {
          for (unsigned int n = 0; n < cell->ids.size(); n++) {
            id = cell->ids[n]; 
            found = false;

            for (unsigned int m = 0; m < ids.size(); m++) {
              if (id == ids[m]) {
                found = true;
                break;
                }
              }

            if (!found) {
              ids.push_back(id); 
              }
            }
          }
        }
      }
    }

  // debug //

  #define ndbg_PmSpatialDecomp_getPoints
  #ifdef dbg_PmSpatialDecomp_getPoints
  string geom_name; 
  PmGraphicsAttributes atts, atts1;
  PmVector3 color;
  PmVector3 verts[100], min_pt, max_pt;
  PmGraphicsLine *line;
  int n;

  min_pt = curr_origin + imin*cell_du*curr_u_axis + jmin*cell_dv*curr_v_axis +
                         kmin*cell_dw*curr_w_axis;

  max_pt = curr_origin + (imax+1)*cell_du*curr_u_axis + (jmax+1)*cell_dv*curr_v_axis +
                         (kmax+1)*cell_dw*curr_w_axis;

  verts[0] = min_pt; 
  verts[1] = verts[0] + (imax-imin+1)*cell_du*curr_u_axis;

  verts[2] = verts[1]; 
  verts[3] = verts[2] + (jmax-jmin+1)*cell_dv*curr_v_axis;

  verts[4] = verts[3]; 
  verts[5] = verts[4] - (imax-imin+1)*cell_du*curr_u_axis;

  verts[6] = verts[5]; 
  verts[7] = verts[6] - (jmax-jmin+1)*cell_dv*curr_v_axis;

  // -- //

  verts[8] = min_pt + (kmax-kmin+1)*cell_dw*curr_w_axis; 
  verts[9] = verts[8] + (imax-imin+1)*cell_du*curr_u_axis;

  verts[10] = verts[9]; 
  verts[11] = verts[10] + (jmax-jmin+1)*cell_dv*curr_v_axis;

  verts[12] = verts[11]; 
  verts[13] = verts[12] - (imax-imin+1)*cell_du*curr_u_axis;

  verts[14] = verts[13]; 
  verts[15] = verts[14] - (jmax-jmin+1)*cell_dv*curr_v_axis;

  // -- //

  verts[16] = verts[0]; 
  verts[17] = verts[8];

  verts[18] = verts[2]; 
  verts[19] = verts[10];

  verts[20] = verts[4]; 
  verts[21] = verts[12];

  verts[22] = verts[6]; 
  verts[23] = verts[14];

  n = 24;

  geom_name = "mm_pot_cell_q";
  color.set(1,1,0);
  atts.setColor(color);
  atts.setDisjoint(true);
  atts.setLineWidth(3.0);
  line = new PmGraphicsLine(geom_name, n, verts);
  line->setAttributes(atts);
  line->display();
  #endif

  }

//*============================================================*
//*==========                mapPoint                ==========*
//*============================================================*
// map a point into a cell index.

void
PmSpatialDecomp::mapPoint(const PmVector3& pt, const float r, int& imin, int& imax, 
                         int& jmin, int& jmax, int& kmin, int& kmax) 
  {
  PmVector3 p;
  float uc, vc, wc;
  float umin, vmin, wmin, umax, vmax, wmax;

  // get indexes into cells //

  p = pt - curr_origin;
  uc = p*curr_u_axis;
  umin = uc - r;
  umax = uc + r;
  mapLocalCoord(umin, cell_du, num_ucell, imin);
  mapLocalCoord(umax, cell_du, num_ucell, imax);

  vc = p*curr_v_axis;
  vmin = vc - r;
  vmax = vc + r;
  mapLocalCoord(vmin, cell_dv, num_vcell, jmin);
  mapLocalCoord(vmax, cell_dv, num_vcell, jmax);

  wc = p*curr_w_axis;
  wmin = wc - r;
  wmax = wc + r;
  mapLocalCoord(wmin, cell_dw, num_wcell, kmin);
  mapLocalCoord(wmax, cell_dw, num_wcell, kmax);
  }

//*============================================================*
//*==========                initialize              ==========*
//*============================================================*

void
PmSpatialDecomp::initialize()
  {

  // compute # cells along each axis //

  num_ucell = static_cast<int>(ceil(u_length / cell_du)); 
  num_vcell = static_cast<int>(ceil(v_length / cell_dv)); 
  num_wcell = static_cast<int>(ceil(w_length / cell_dw)); 
  num_cells = num_ucell * num_vcell * num_wcell;
  cells = new PmSpatialCell[num_cells];

  /*
  fprintf (stderr, "\n>>>>>> PmSpatialDecomp::initialize %s \n", name.c_str());
  fprintf (stderr, ">>> axes lengths = %f %f %f\n", u_length, v_length, w_length); 
  fprintf (stderr, ">>> cell widths  = %f %f %f\n", cell_du, cell_dv, cell_dw); 
  fprintf (stderr, ">>> cell nums = %d %d %d \n", num_ucell, num_vcell, num_wcell); 
  fprintf (stderr, ">>> num cells = %d \n", num_cells); 
  */
  }

//*============================================================*
//*==========                printInfo               ==========*
//*============================================================*

void
PmSpatialDecomp::printInfo()
  {
  fprintf (stderr, "\n>>>>>> PmSpatialDecomp::printInfo \n");
  int n = 0;

  for (int k = 0; k < num_wcell; k++) {
    for (int j = 0; j < num_vcell; j++) {
      for (int i = 0; i < num_ucell; i++) {
        fprintf (stderr, "%d %d %d: ", i, j, k);

        if (cells[n].ids.size() != 0) {
          fprintf (stderr, "%d ", cells[n].ids.size());
          }

        n += 1;
        fprintf (stderr, "\n");
        }
      }
    }
  }

//*============================================================*
//*==========                setXform                ==========*
//*============================================================*

void 
PmSpatialDecomp::setXform(PmXform& xform)
  {
  PmMatrix3x3 mat;

  // transform origin and axes //

  mat = xform.matrix;
  curr_origin = mat*(origin-xform.center) + xform.translation + xform.center;

  // transform axes //

  curr_u_axis = mat*u_axis;
  curr_v_axis = mat*v_axis;
  curr_w_axis = mat*w_axis;
  this->xform = xform;
  //display();
  }

//*============================================================*
//*==========                setXform                ==========*
//*============================================================*

void
PmSpatialDecomp::display()
  {
  PmVector3 verts[100], pt, verts1[10000], lvec, dvec;
  PmVector3 lu, lv, lw;
  int n, n1;

  lu = u_length*curr_u_axis;
  lv = v_length*curr_v_axis;
  lw = w_length*curr_w_axis;

  //fprintf (stderr, " curr_u_axis*curr_v_axis = %f \n", curr_u_axis*curr_v_axis); 
  //fprintf (stderr, " curr_u_axis*curr_w_axis = %f \n", curr_u_axis*curr_w_axis); 
  //fprintf (stderr, " curr_v_axis*curr_w_axis = %f \n", curr_v_axis*curr_w_axis); 

  verts[0] = curr_origin; 
  verts[1] = verts[0] + lu; 
  verts[2] = verts[1]; 
  verts[3] = verts[2] + lv; 
  verts[4] = verts[3]; 
  verts[5] = verts[4] - lu; 
  verts[6] = verts[5]; 
  verts[7] = verts[6] - lv; 

  verts[8] = curr_origin + lw; 
  verts[9] = verts[8] + lu;
  verts[10] = verts[9];
  verts[11] = verts[10] + lv;
  verts[12] = verts[11];
  verts[13] = verts[12] - lu;
  verts[14] = verts[13];
  verts[15] = verts[14] - lv;

  verts[16] = verts[0]; 
  verts[17] = verts[8]; 
  verts[18] = verts[1]; 
  verts[19] = verts[9]; 
  verts[20] = verts[3]; 
  verts[21] = verts[11]; 
  verts[22] = verts[5]; 
  verts[23] = verts[13]; 
  verts[24] = verts[7]; 
  verts[25] = verts[15]; 
  n = 24;

  n1 = 0;

  for (int k = 0; k < num_wcell; k++) {
    pt = curr_origin + k*cell_dw*curr_w_axis;

    for (int i = 0; i < num_ucell; i++) {
      verts1[n1++] = pt; 
      verts1[n1++] = pt + lv; 
      pt = pt + cell_du*curr_u_axis;
      }
    }

  for (int k = 0; k < num_wcell; k++) {
    pt = curr_origin + k*cell_dw*curr_w_axis;

    for (int j = 0; j < num_vcell; j++) {
      verts1[n1++] = pt; 
      verts1[n1++] = pt + lu; 
      pt = pt + cell_dv*curr_v_axis;
      }
    }

  for (int i = 0; i < num_ucell; i++) {
    pt = curr_origin + i*cell_du*curr_u_axis;

    for (int j = 0; j < num_vcell; j++) {
      verts1[n1++] = pt;
      verts1[n1++] = pt + lw;
      pt = pt + cell_dv*curr_v_axis;
      }
    }

  if (!line) {
    string geom_name, geom_name1;
    PmGraphicsAttributes atts, atts1;
    PmVector3 color;

    geom_name = "mm_pot_coord_sys";
    color.set(1,0.6,0.6);
    atts.setColor(color);
    atts.setDisjoint(true);
    atts.setLineWidth(3.0);
    line = new PmGraphicsLine(geom_name, n, verts);
    line->setAttributes(atts);
    line->display();

    geom_name1 = "mm_pot_cells";
    color.set(0.6, 0.6, 0.6);
    atts1.setColor(color);
    atts1.setDisjoint(true);
    line1 = new PmGraphicsLine(geom_name1, n1, verts1);
    line1->setAttributes(atts1);
    line1->display();
    }

  else {
    line->update(n, verts);
    line->display();
    line1->update(n1, verts1);
    line1->display();
    }
  }

}


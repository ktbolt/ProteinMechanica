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
//* grid:                 g r i d                              *
//*============================================================*

#include "grid.h"
#include "grid_prv.h"
#include "pm/mth.h"

namespace ProteinMechanica {


////////////////////////////////////////////////////////////////
//                    p u b l i c                            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGrid::PmGrid(const string name) {
  this->name = name;
  this->init();
  }

PmGrid::PmGrid(const string name, const int idim, const int jdim, const int kdim,
               const float dx, const float dy, const float dz, float *data)
  {
  this->name = name;
  this->init();
  this->idim = idim;
  this->jdim = jdim;
  this->kdim = kdim;

  this->cell_dx = dx;
  this->cell_dy = dy;
  this->cell_dz = dz;

  this->grid_data = data;
  data_vmin = data_vmax = data[0];

  for (int i = 0; i < idim*jdim*kdim; i++) {
    if (data[i] > data_vmax) data_vmax = data[i];
    if (data[i] < data_vmin) data_vmin = data[i];
    }

  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  }

void
PmGrid::init() 
  {
  origin.set(0,0,0);
  basis_u.set(1,0,0); 
  basis_v.set(0,1,0); 
  basis_w.set(0,0,1);
  idim = 0;
  jdim = 0;
  kdim = 0;
  cell_dx = cell_dy = cell_dz = 0.0;
  color.set(1,1,1);
  display_type = PM_GEOMETRY_DISPLAY_SOLID;
  shading_type = PM_GEOMETRY_SHADING_FLAT;
  }

PmGrid::~PmGrid() {
  }

//*============================================================*
//*==========          getDataRange                  ==========*
//*============================================================*

void 
PmGrid::getDataRange(float& vmin, float& vmax)
  {
  vmin = data_vmin; 
  vmax = data_vmax;
  }

//*============================================================*
//*==========          get/set extent                ==========*
//*============================================================*

void 
PmGrid::getExtent (PmExtent& extent) { extent = this->extent; }

void 
PmGrid::setExtent (const PmExtent& extent) { this->extent = extent; }

//*============================================================*
//*==========          get/set name                  ==========*
//*============================================================*

void 
PmGrid::getName (string& name) { name = this->name;  }

//*============================================================*
//*==========               interpScalar             ==========*
//*============================================================*
// interpolate the grid to the given point.

void
PmGrid::interpScalar(PmVector3& pt, int& i, int& j, int& k, float& ival)
  {
  float xmin, ymin, zmin;

  float x0, y0, z0;

  float fx, fy, fz;

  float v04, v37, v26, v15;

  float vl, vr;

  float cell_data[8];

  //fprintf (stderr, "\n>>>>>> PmGrid::interpScalar \n");

  // get the cell index containing the point //

  getIndex(pt, i, j, k);

  xmin = origin[0];
  ymin = origin[1];
  zmin = origin[2];

  x0 = xmin + i*cell_dx;
  y0 = ymin + j*cell_dy;
  z0 = zmin + k*cell_dz;

  fx = (pt[0] - x0) / cell_dx;
  fy = (pt[1] - y0) / cell_dy;
  fz = (pt[2] - z0) / cell_dz;
  //fprintf (stderr, ">>> dx %f  dy %f  dz %f \n", dx, dy, dz);
  //fprintf (stderr, ">>> fx %f  fy %f  fz %f \n", fx, fy, fz);

  this->getCellData(i, j, k, cell_data);
  /*
  fprintf (stderr, ">>> cell data: ");

  for (int i = 0; i < 8; i++) {
    fprintf (stderr, "%f ", cell_data[i]);
    }

  fprintf (stderr, "\n");
  */

  v04 = cell_data[0] + fz * (cell_data[4] - cell_data[0]);
  v37 = cell_data[3] + fz * (cell_data[7] - cell_data[3]);
  v15 = cell_data[1] + fz * (cell_data[5] - cell_data[1]);
  v26 = cell_data[2] + fz * (cell_data[6] - cell_data[2]);

  vl = v04 + fy * (v37 - v04);
  vr = v15 + fy * (v26 - v15);
  ival = vl + fx * (vr - vl);
  }

//*============================================================*
//*==========               getCellData              ==========*
//*============================================================*
// get the data for a cell.                                  

void 
PmGrid::getCellData(int i, int j, int k, float cell_data[8])
  {

  if ((i < 0) || (i >= idim-1)) {
    return;
    }

  if ((j < 0) || (j >= jdim-1)) {
    return;
    }

  if ((k < 0) || (k >= kdim-1)) {
    return;
    }

  cell_data[0] = grid_data[i + idim*(j + jdim*(k))];
  cell_data[1] = grid_data[(i+1) + idim*(j + jdim*(k))];
  cell_data[2] = grid_data[(i+1) + idim*((j+1) + jdim*(k))];
  cell_data[3] = grid_data[(i) + idim*((j+1) + jdim*(k))];

  cell_data[4] = grid_data[i + idim*(j + jdim*(k+1))];
  cell_data[5] = grid_data[(i+1) + idim*(j + jdim*(k+1))];
  cell_data[6] = grid_data[(i+1) + idim*((j+1) + jdim*(k+1))];
  cell_data[7] = grid_data[(i) + idim*((j+1) + jdim*(k+1))];
  }

//*============================================================*
//*==========               getIndex                 ==========*
//*============================================================*
// get an index from a point.     

void 
PmGrid::getIndex(PmVector3& pt, int& i, int& j, int& k) 
  {

  float x0, y0, z0;
  float x, y, z;
  PmVector3 origin;
  PmVector3 u, v, w;
  float dx, dy, dz, dp;

  i = -1;
  j = -1;
  k = -1;

  x = pt[0];
  y = pt[1];
  z = pt[2];

  u.set(1,0,0);
  v.set(0,1,0);
  w.set(0,0,1);

  x0 = 0.0; 
  y0 = 0.0;
  z0 = 0.0;
  dx = x - x0;
  dy = y - y0;
  dz = z - z0;

  /*  compute indexes into mesh.  */

  dp = dx*u[0] + dy*u[1] + dz*u[2];
  i = (int)(dp / this->cell_dx);
  if (i < 0) i = 0;
  if (i > this->idim - 1) i = this->idim - 1;

  dp = dx*v[0] + dy*v[1] + dz*v[2];
  j = (int)(dp / this->cell_dy);
  if (j < 0) j = 0;
  if (j > this->jdim - 1) j = this->jdim - 1;

  dp = dx*w[0] + dy*w[1] + dz*w[2];
  k = (int)(dp / this->cell_dz);
  if (k < 0) k = 0;
  if (k > this->kdim - 1) k = this->kdim - 1;
  }

////////////////////////////////////////////////////////////////
//                    g r a p h i c s                        //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              buildGeomName             ==========*
//*============================================================*

void
PmGrid::buildGeomName (string& grname) {
  grname = "grid";
  grname += '[' + name + ']';
  }

//*============================================================*
//*==========             addGeometry                ==========*
//*============================================================*
// add graphics geometry object.

void
PmGrid::addGeometry (PmGraphicsPoint *geom) {
  graphics_point_geometry.push_back (geom);
  }

//*============================================================*
//*==========               getGeometry              ==========*
//*============================================================*
// get a graphics geometry object.

void
PmGrid::getGeometry (string name, PmGraphicsPoint **geom)
  {
  *geom = NULL;

  for (unsigned int i = 0; i < graphics_point_geometry.size(); i++) {
    if (name == graphics_point_geometry[i]->name) {
      *geom = graphics_point_geometry[i];
      }
    }
  }

//*============================================================*
//*==========          displayPoints                 ==========*
//*============================================================*

void 
PmGrid::displayPoints(const float vmin, const float vmax, PmDbPmInterface *db)
  {
  fprintf (stderr, "\n>>>>>> PmGrid::displayPoints \n");
  int num_verts;
  PmVector3 *verts;
  float val, x, y, z;
  string geom_name;
  PmGraphicsPoint *geom; 
  PmGraphicsAttributes atts;
  PmVector3 color, *gcolors;
  vector<PmVector3> colors;
  vector<float> vals;

  num_verts = 0;
  fprintf (stderr, "   >>> idim %d \n", idim);
  fprintf (stderr, "   >>> jdim %d \n", jdim);
  fprintf (stderr, "   >>> kdim %d \n", kdim);
  //fprintf (stderr, "   >>> data %x \n", grid_data);
  fprintf (stderr, "   >>> vmin %f   vmax %f \n", vmin, vmax);
  fprintf (stderr, "   >>> data_vmin %f   data_vmax %f \n", data_vmin, data_vmax);

  for (int i = 0; i < idim*jdim*kdim; i++) {
    val = grid_data[i];
    //fprintf (stderr, "   >>> %d val %f \n", i, val);

    if ((val <= vmax) && (val >= vmin)) {
      num_verts += 1;
      }
    }

  fprintf (stderr, "   >>> num verts = %d \n", num_verts);

  if (!num_verts) {
    return;
    }

  buildGeomName (geom_name);
  getGeometry (geom_name, &geom);

  verts = new PmVector3[num_verts];
  num_verts = 0;
  z = 0.0;

  for (int k = 0; k < kdim; k++) {
    y = 0.0;

    for (int j = 0; j < jdim; j++) {
      x = 0.0;

      for (int i = 0; i < idim; i++) {
        //val = grid_data[i+(jdim-1)*j + (jdim-1)*(kdim-1)*k];
        val = grid_data[i+(idim)*j + (idim)*(jdim)*k];

        if ((val <= vmax) && (val >= vmin)) {
          verts[num_verts][0] = x;
          verts[num_verts][1] = y;
          verts[num_verts][2] = z;
          vals.push_back(val);
          num_verts += 1;
          }

        x += cell_dx;
        }

      y += cell_dy;
      }

    z += cell_dz;
    }

  // map data values into a color map //

  grid_MapData(vmin, vmax, vals, colors);
  gcolors = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gcolors[i] = colors[i];
    }

  color.set(1,1,1);
  atts.setColor(color);
  atts.setMarker(false);
  atts.setScale(1.0);

  geom = new PmGraphicsPoint(geom_name, num_verts, verts);
  geom->setAttributes(atts); 
  geom->setColors(num_verts, gcolors); 
  geom->display(); 

  // write points //

  if (db) {
    vector<PmVector3> pverts; 

    for (int i = 0; i < num_verts; i++) {
      pverts.push_back(verts[i]);
      }

    db->writeParticle(this->name, pverts, vals);
    }
  }

//*============================================================*
//*==========             displayVolume              ==========*
//*============================================================*

void 
PmGrid::displayVolume(const float vmin, const float vmax, const int num_samples)
  {
  fprintf (stderr, "\n>>>>>> PmGrid::displayVolume\n");
  int num_verts;
  vector<PmVector3> verts;
  float val, ox, oy, oz, x0, y0, z0, rx, ry, rz;
  float ux, uy, uz, vx, vy, vz, wx, wy, wz;
  float v04, v37, v26, v15, vl, vr;
  string geom_name;
  PmGraphicsPoint *geom; 
  PmGraphicsAttributes atts;
  PmVector3 color, *gcolors;
  vector<PmVector3> colors;
  vector<float> vals;
  float cell_data[8]; 
  int long seed = 1353971;
  PmVector3 pt, *gverts;

  ux = cell_dx;
  uy = 0.0;
  uz = 0.0;

  vx = 0.0;
  vy = cell_dy;
  vz = 0.0;

  wx = 0.0;
  wy = 0.0;
  wz = cell_dz;

  for (int k = 0; k < kdim-1; k++) {
    for (int j = 0; j < jdim-1; j++) {
      for (int i = 0; i < idim-1; i++) {
        this->getCellData(i, j, k, cell_data);
        x0 = ox + i*ux + vx*j + wx*k;
        y0 = oy + i*uy + vy*j + wy*k;
        z0 = oz + i*uz + vz*j + wz*k;

        for (int n = 0; n < num_samples; n++) {
          rx = pm_MathRandNumber(&seed);
          ry = pm_MathRandNumber(&seed);
          rz = pm_MathRandNumber(&seed);

          v04 = cell_data[0] + rz*(cell_data[4] - cell_data[0]);
          v37 = cell_data[3] + rz*(cell_data[7] - cell_data[3]);
          v15 = cell_data[1] + rz*(cell_data[5] - cell_data[1]);
          v26 = cell_data[2] + rz*(cell_data[6] - cell_data[2]);

          vl = v04 + ry*(v37 - v04);
          vr = v15 + ry*(v26 - v15);
          val = vl + rx*(vr - vl);

          if ((val >= vmin) && (val <= vmax)) {
            pt[0] = x0 + rx*cell_dx;
            pt[1] = y0 + ry*cell_dy;
            pt[2] = z0 + rz*cell_dz;
            verts.push_back(pt);
            vals.push_back(val);
            }
          }
        }
      }
    }

  num_verts = verts.size();

  if (!num_verts) {
    return;
    }

  gverts = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gverts[i] = verts[i];
    }

  // map data values into a color map //

  grid_MapData(vmin, vmax, vals, colors);
  gcolors = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gcolors[i] = colors[i];
    }

  buildGeomName (geom_name);
  getGeometry (geom_name, &geom);

  color.set(1,1,1);
  atts.setColor(color);
  atts.setMarker(false);
  atts.setScale(1.0);

  geom = new PmGraphicsPoint(geom_name, num_verts, gverts);
  geom->setAttributes(atts); 
  geom->setColors(num_verts, gcolors); 
  geom->display(); 
  }

//*============================================================*
//*==========          displayIsosurface             ==========*
//*============================================================*

void 
PmGrid::displayIsosurface(const float level, bool seeded, PmVector3 point, 
                          PmDbPmInterface *db)
  {
  int dims[3]; 
  float cell_sizes[3]; 
  vector<PmVector3> verts; 
  vector<int> polys;
  PmGraphicsAttributes atts;
  PmGraphicsPolygon *poly;
  int num_tri, num_verts; 
  PmVector3 *gverts;
  PmConn *conn;
  int iseed[3];
  float slevel;

  dims[0] = idim; dims[1] = jdim; dims[2] = kdim;
  cell_sizes[0] = cell_dx; cell_sizes[1] = cell_dy; cell_sizes[2] = cell_dz; 

  if (seeded) {
    this->getIndex(point, iseed[0], iseed[1], iseed[2]);
    interpScalar(point, iseed[0], iseed[1], iseed[2], slevel); 
    grid_IsoSeededExtract (dims, cell_sizes, grid_data, slevel, iseed, verts, polys);
    }
  else {
    grid_IsoExtract (dims, cell_sizes, grid_data, level, verts, polys);
    }

  num_verts = verts.size();
  gverts = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gverts[i] = verts[i];
    } 

  num_tri = polys.size() / 3;
  conn = new PmConn(3*num_tri);
  int *pconn = conn->vals;

  for (int i = 0; i < 3*num_tri; i++) {
    pconn[i] = polys[i];
    } 

  poly = new PmGraphicsPolygon(name, num_tri, conn, 3, num_verts, gverts);

  /*
  PmGraphicsSurface *geom = new PmGraphicsSurface (name, num_verts, gverts, num_tri, 
                                                   conn, 3);
  */

  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setColor(color);
  poly->setAttributes(atts);
  poly->display();

  // write isosurface //

  if (db) {
    db->writeSurface(this->name, verts, polys);
    }
  }

//*============================================================*
//*==========          grid_IsoExtract               ==========*
//*============================================================*

static void
grid_IsoExtract (int dims[3], float cell_sizes[3], float *data, float level,
                 vector<PmVector3>& verts, vector<int>& polys) 
  {

  int i, j, k, pindex, num_polys;

  int num_verts;

  int loc, edge_num, v, v1, v2;

  PmVector3 pt1, pt2, vert;

  float d1, d2;

  int num_pts, iverts_num[100], num_edges, *edges;

  int num_ipolys, *ipolys, polyi[4];

  float t, x, y, z;

  int xdim, ydim, zdim;

  float xmin, ymin, zmin, dx, dy, dz;

  int xi, yi, zi;

  float xp[8], yp[8], zp[8], dp[8]; 

  float xf, yf, zf, xf1, yf1, zf1;

  int j_off, k_off, j_off_1, k_off_1, i_begin, i_end, j_begin, j_end, k_begin, k_end;

  int edge_pt_table_size;

  EdgePt **edge_pt_table, *edge_pt_ptr;

  int min_v, max_v, index, vert_num;

  int np[8], n1, n2;

  int num_malloc; 

  int ptr_block_size; 

  EdgePt *ptr_block;

  EdgePtBlockList *block_list, *ptr;

  int k_inc;

  int conn_size, conn_inc, verts_size, verts_inc;

  int num_pmalloc, num_vmalloc;

#define debug_grid_IsoExtract
#ifdef debug_grid_IsoExtract
  fprintf (stderr, "\n  --------  grid_IsoExtract --------  \n");
#endif

  grid_IsoInit();

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  dx = cell_sizes[0]; 
  dy = cell_sizes[1]; 
  dz = cell_sizes[2]; 

  xmin = 0.0;
  ymin = 0.0;
  zmin = 0.0;

  i_begin = 0; 
  i_end = xdim;
  j_begin = 0; 
  j_end = ydim; 
  k_begin = 0; 
  k_end = zdim; 
  k_inc = xdim * ydim; 

  num_pmalloc = num_vmalloc = 0;

  edge_pt_table_size = 1000;
  edge_pt_table_size = 10000;
  edge_pt_table = new EdgePt*[edge_pt_table_size];

  num_malloc = 0; 
  ptr_block_size = 20000; 
  ptr_block_size = 2000; 
  ptr_block = new EdgePt[ptr_block_size];

  block_list = new EdgePtBlockList;
  block_list->size = ptr_block_size; 
  block_list->ptr_block = ptr_block; 
  block_list->next_block = NULL; 

  num_polys = 0;
  num_verts = 0;
  pindex = 0;

  for (i = 0; i < edge_pt_table_size; i++) {
    edge_pt_table[i] = NULL; 
    }

  zf = zmin + (float)k_begin * dz;
  zf1 = zf + dz;
  k_off = k_inc * k_begin;

  for (zi = k_begin; zi < k_end-1 ; zi++, k_off += k_inc) {
    yf = ymin + (float)j_begin * dy;
    yf1 = yf + dy;
    k_off_1 = k_off + k_inc; 

    for (yi = j_begin; yi < j_end-1; yi++) {
      xf = xmin + dx*i_begin;
      xf1 = xf + dx;
      j_off = xdim * yi;
      j_off_1 = xdim * (yi + 1);

      for (xi = i_begin; xi < i_end-1; xi++) {
        dp[0] = data[xi + j_off + k_off_1]; 
        dp[1] = data[(xi+1) + j_off + k_off_1]; 
        dp[2] = data[(xi+1) + j_off_1 + k_off_1]; 
        dp[3] = data[(xi) + j_off_1 + k_off_1]; 
        dp[4] = data[xi + j_off + k_off]; 
        dp[5] = data[(xi+1) + j_off + k_off]; 
        dp[6] = data[(xi+1) + j_off_1 + k_off]; 
        dp[7] = data[(xi) + j_off_1 + k_off]; 

        loc = ((dp[7] >= level) << 7) |
              ((dp[6] >= level) << 6) |
              ((dp[5] >= level) << 5) |
              ((dp[4] >= level) << 4) |
              ((dp[3] >= level) << 3) |
              ((dp[2] >= level) << 2) |
              ((dp[1] >= level) << 1) |
              (dp[0] >= level); 

        if ((loc > 0) && (loc < 255)) {
          edges = iso_edge_table[loc];
          num_edges = edges[0];

          xp[0] = xf; 
          xp[1] = xf1; 
	  xp[3] = xp[4] = xp[7] = xp[0];
          xp[2] = xp[5] = xp[6] = xp[1];

          yp[0] = yf; 
          yp[2] = yf1; 
	  yp[1] = yp[4] = yp[5] = yp[0];
          yp[3] = yp[6] = yp[7] = yp[2];

          zp[0] = zf1; 
          zp[4] = zf; 
          zp[1] = zp[2] = zp[3] = zp[0];
          zp[5] = zp[6] = zp[7] = zp[4];

          np[0] = xi + j_off + k_off_1; 
          np[1] = (xi+1) + j_off + k_off_1; 
          np[2] = (xi+1) + j_off_1 + k_off_1; 
          np[3] = (xi) + j_off_1 + k_off_1; 
          np[4] = xi + j_off + k_off; 
          np[5] = (xi+1) + j_off + k_off; 
          np[6] = (xi+1) + j_off_1 + k_off; 
          np[7] = (xi) + j_off_1 + k_off; 

          for (i = 1; i <= num_edges; i++) {
            edge_num = edges[i];
            v1 = edge_list[edge_num][0];
            v2 = edge_list[edge_num][1];

            pt1[0] = xp[v1]; 
            pt1[1] = yp[v1];
            pt1[2] = zp[v1];
            pt2[0] = xp[v2]; 
            pt2[1] = yp[v2];
            pt2[2] = zp[v2];

            d1 = dp[v1];
            d2 = dp[v2];
            n1 = np[v1];
            n2 = np[v2];

	    if (n1 < n2) {
	      min_v = n1;
	      max_v = n2;
	      }
            else {
	      min_v = n2;
	      max_v = n1;
	      }

            index = min_v % edge_pt_table_size;
            edge_pt_ptr = NULL; 

            if (edge_pt_table[index] != NULL) {
	      edge_pt_ptr = edge_pt_table[index];

              while (edge_pt_ptr != NULL) {
                if (edge_pt_ptr->node == max_v) {
		  break;
		  }

		edge_pt_ptr = edge_pt_ptr->next_pt;
		}
	      }

            if (edge_pt_ptr == NULL) {
	      edge_pt_ptr = &ptr_block[num_malloc]; 
              num_malloc++; 

              if (num_malloc == ptr_block_size) { 
                num_malloc = 0; 
                ptr_block = new EdgePt[ptr_block_size];
                ptr = new EdgePtBlockList;
                ptr->size = ptr_block_size; 
                ptr->ptr_block = ptr_block; 
                ptr->next_block = block_list->next_block; 
		block_list = ptr;
		}

              t = (level - d1) / (d2 - d1);
              x = pt1[0] + t * (pt2[0] - pt1[0]);
              y = pt1[1] + t * (pt2[1] - pt1[1]);
              z = pt1[2] + t * (pt2[2] - pt1[2]);

              vert[0] = x;
              vert[1] = y;
              vert[2] = z;
              verts.push_back(vert);

              edge_pt_ptr->node = max_v;
	      edge_pt_ptr->index = num_verts;
	      vert_num = num_verts;
	      num_verts++;

              if (num_verts > (verts_size - 20)) {
                verts_size += verts_inc;
                // grow verts //
                num_vmalloc += 1;
		}

              edge_pt_ptr->next_pt = edge_pt_table[index];
              edge_pt_table[index] = edge_pt_ptr;
	      }

	    else {
              vert_num = (int)edge_pt_ptr->index; 
	      }

            iverts_num[i-1] = vert_num;
	    }

          ipolys = iso_poly_table[loc];
          num_ipolys = ipolys[0];

          for (i = 1, j = 1; i <= num_ipolys; i++) {
            num_pts = ipolys[j++];
            //polys.push_back(num_pts);
            //polys[pindex++] = num_pts;
	    num_polys++;
        
            for (k = 0; k < num_pts; k++) {
              v = ipolys[j++];
              polyi[k] = iverts_num[v];
              //polys.push_back(iverts_num[v]);
	      //polys[pindex++] = iverts_num[v];
              }

            polys.push_back(polyi[0]);
            polys.push_back(polyi[2]);
            polys.push_back(polyi[1]);
	    }

          if (pindex > (conn_size - 20)) {
            conn_size += conn_inc;
            //vs_GeometryPolyConnGrow (geom, mem_type, conn_size, DM_TRUE, &polys);
            num_pmalloc += 1;
	    }
	  }

        xf += dx;
        xf1 += dx;
        }

      yf += dy;
      yf1 += dy;
      }

    zf += dz;
    zf1 += dz;
    }

#ifdef use

  delete edge_pt_table;


  /*  free block ptrs.  */

  while (block_list != DmNil(VsEdgePtBlockList *)) {
    dm_MemFree (0, block_list->ptr_block);
    block_list = block_list->next_block;
    }

#endif

#ifdef debug_grid_IsoExtract 
  fprintf (stderr,"\n  num_polys  [%d]  \n", num_polys); 
  fprintf (stderr,"  num_verts  [%d]  \n", num_verts); 
  fprintf (stderr,"  pindex     [%d]  \n", pindex); 
  fprintf (stderr,"  num malloc [%d]  \n", num_malloc);
  fprintf (stderr,"  num poly malloc [%d]  \n", num_pmalloc);
  fprintf (stderr,"  num vert malloc [%d]  \n", num_vmalloc);
#endif
  }

//*============================================================*
//*==========          grid_IsoSeededExtract         ==========*
//*============================================================*

static void
grid_IsoSeededExtract (int dims[3], float cell_sizes[3], float *data, float level,
                       int iseed[3], vector<PmVector3>& verts, vector<int>& polys)
  {

  int i, j, k, l, xi, yi, zi, new_xi, new_yi, new_zi;

  int pindex, num_polys, num_verts, loc, edge_num, v, v1, v2;

  PmVector3 pt1, pt2;

  int num_pts, iverts_num[100], num_edges, *edges;

  int num_ipolys, *ipolys, polyi[10];

  float t, x, y, z, d1, d2;

  int xdim, ydim, zdim;

  float xmin, ymin, zmin, dx, dy, dz;

  float xp[8], yp[8], zp[8], dp[8];

  float xf, yf, zf, xf1, yf1, zf1;

  int j_off, k_off, i_begin, i_end, j_begin, j_end, k_begin, k_end;

  int edge_pt_table_size;

  EdgePt **edge_pt_table, *edge_pt_ptr, *ptr_block;

  int min_v, max_v, index, vert_num, np[8], n1, n2;

  int num_malloc, num_pmalloc, num_vmalloc, ptr_block_size; 

  EdgePtBlockList *block_list, *ptr;

  int conn_size, conn_inc, verts_size, verts_inc, list_size;

  int i_list[5000], j_list[5000], k_list[5000];

  int new_list_size, new_i_list[5000], new_j_list[5000], new_k_list[5000];

  unsigned long *blocks;

  int size, num_blocks, new_cell; 

  bool found;

  PmVector3 vec;

  fprintf (stderr, "\n ----------  grid_IsoSeedExtract ---------- \n");
  fprintf (stderr, ">>> level = %f \n", level);

  grid_IsoInit();

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  dx = cell_sizes[0];
  dy = cell_sizes[1];
  dz = cell_sizes[2];

  xmin = 0.0;
  ymin = 0.0;
  zmin = 0.0;

  i_begin = 0;
  i_end = xdim;
  j_begin = 0;
  j_end = ydim;
  k_begin = 0;
  k_end = zdim;
  k_off = xdim * ydim;
  j_off = xdim;

  num_pmalloc = num_vmalloc = 0;

  edge_pt_table_size = 1000;
  edge_pt_table_size = 10000;
  edge_pt_table = new EdgePt*[edge_pt_table_size];

  num_malloc = 0; 
  ptr_block_size = 20000; 
  ptr_block_size = 2000; 
  ptr_block = new EdgePt[ptr_block_size];

  block_list = new EdgePtBlockList;
  block_list->size = ptr_block_size; 
  block_list->ptr_block = ptr_block; 
  block_list->next_block = NULL; 

  num_polys = 0;
  num_verts = 0;
  pindex = 0;

  for (i = 0; i < edge_pt_table_size; i++) {
    edge_pt_table[i] = NULL; 
    }

  /*  initialize bit mask blocks.  */

  size = xdim * ydim * zdim;
  grid_SubsetCreate (size, &num_blocks, &blocks);


  // start with initial (seed) cell. // 

  xi = iseed[0];
  yi = iseed[1];
  zi = iseed[2];
  fprintf (stderr,">>>  seed index = %d %d %d  \n", xi, yi, zi); 

  list_size = 0;
  i_list[0] = xi;
  j_list[0] = yi;
  k_list[0] = zi;

  new_list_size = 1;
  new_i_list[0] = xi;
  new_j_list[0] = yi;
  new_k_list[0] = zi;

  while (new_list_size != 0) {
    for (l = 0; l < new_list_size; l++) {
      i_list[l] = new_i_list[l];
      j_list[l] = new_j_list[l];
      k_list[l] = new_k_list[l];
      }

    list_size = new_list_size;
    new_list_size = 0;

    for (l = 0; l < list_size; l++) {
      xi = i_list[l];
      yi = j_list[l];
      zi = k_list[l];

      dp[0] = data[xi + j_off*yi + k_off*(zi + 1)];
      dp[1] = data[(xi+1) + j_off*yi + k_off*(zi + 1)];
      dp[2] = data[(xi+1) + j_off*(yi+1) + k_off*(zi + 1)];
      dp[3] = data[(xi) + j_off*(yi+1) + k_off*(zi + 1)];
      dp[4] = data[xi + j_off*yi + k_off*zi];
      dp[5] = data[(xi+1) + j_off*yi + k_off*zi];
      dp[6] = data[(xi+1) + j_off*(yi+1) + k_off*zi];
      dp[7] = data[(xi) + j_off*(yi+1) + k_off*zi];

      loc = ((dp[7] >= level) << 7) |
             ((dp[6] >= level) << 6) |
             ((dp[5] >= level) << 5) |
             ((dp[4] >= level) << 4) |
             ((dp[3] >= level) << 3) |
             ((dp[2] >= level) << 2) |
             ((dp[1] >= level) << 1) |
             (dp[0] >= level); 

      /*
      fprintf (stderr, " >>>>>> loc [%d] \n", loc);
      */

      if ((loc > 0) && (loc < 255)) {
        edges = iso_edge_table[loc];
        num_edges = edges[0];
        xf = xmin + xi*dx;
        yf = ymin + yi*dy;
        zf = zmin + zi*dz;
        xf1 = xf + dx;
        yf1 = yf + dy;
        zf1 = zf + dz;

        xp[0] = xf; 
        xp[1] = xf1; 
	xp[3] = xp[4] = xp[7] = xp[0];
        xp[2] = xp[5] = xp[6] = xp[1];

        yp[0] = yf; 
        yp[2] = yf1; 
	yp[1] = yp[4] = yp[5] = yp[0];
        yp[3] = yp[6] = yp[7] = yp[2];

        zp[0] = zf1; 
        zp[4] = zf; 
        zp[1] = zp[2] = zp[3] = zp[0];
        zp[5] = zp[6] = zp[7] = zp[4];

        np[0] = xi     + j_off*yi     + k_off*(zi + 1);
        np[1] = (xi+1) + j_off*yi     + k_off*(zi + 1);
        np[2] = (xi+1) + j_off*(yi+1) + k_off*(zi + 1);
        np[3] = (xi)   + j_off*(yi+1) + k_off*(zi + 1);
        np[4] = xi     + j_off*yi     + k_off*zi;
        np[5] = (xi+1) + j_off*yi     + k_off*zi;
        np[6] = (xi+1) + j_off*(yi+1) + k_off*zi;
        np[7] = (xi)   + j_off*(yi+1) + k_off*zi;

        for (i = 1; i <= num_edges; i++) {
          edge_num = edges[i];
          v1 = edge_list[edge_num][0];
          v2 = edge_list[edge_num][1];
          pt1[0] = xp[v1]; pt1[1] = yp[v1]; pt1[2] = zp[v1];
          pt2[0] = xp[v2]; pt2[1] = yp[v2]; pt2[2] = zp[v2];
          d1 = dp[v1];
          d2 = dp[v2];
          n1 = np[v1];
          n2 = np[v2];

	  if (n1 < n2) {
	    min_v = n1;
	    max_v = n2;
	    }
          else {
	    min_v = n2;
	    max_v = n1;
	    }

          index = min_v % edge_pt_table_size;
          edge_pt_ptr = NULL;

          if (edge_pt_table[index] != NULL) {
	    edge_pt_ptr = edge_pt_table[index];

            while (edge_pt_ptr != NULL) {
              if (edge_pt_ptr->node == max_v) {
                break;
	        }

              edge_pt_ptr = edge_pt_ptr->next_pt;
	      }
	    }

          if (edge_pt_ptr == NULL) {
	    edge_pt_ptr = &ptr_block[num_malloc]; 
            num_malloc++; 

            if (num_malloc == ptr_block_size) { 
              num_malloc = 0; 
              ptr_block = new EdgePt[ptr_block_size];
              ptr = new EdgePtBlockList;
              ptr->size = ptr_block_size; 
              ptr->ptr_block = ptr_block; 
              ptr->next_block = block_list->next_block; 
	      block_list = ptr;
	      }

            t = (level - d1) / (d2 - d1);
            x = pt1[0] + t * (pt2[0] - pt1[0]);
            y = pt1[1] + t * (pt2[1] - pt1[1]);
            z = pt1[2] + t * (pt2[2] - pt1[2]);

            vec.set(x,y,z);
            verts.push_back(vec);

            edge_pt_ptr->node = max_v;
	    edge_pt_ptr->index = num_verts;
	    vert_num = num_verts;
	    num_verts++;

            if (num_verts > (verts_size - 20)) {
              verts_size += verts_inc;
              num_vmalloc += 1;
	      }

            edge_pt_ptr->next_pt = edge_pt_table[index];
            edge_pt_table[index] = edge_pt_ptr;
	    }

	  else {
            vert_num = (int)edge_pt_ptr->index; 
	    }

          iverts_num[i-1] = vert_num;
	  }

        ipolys = iso_poly_table[loc];
        num_ipolys = ipolys[0];

        for (i = 1, j = 1; i <= num_ipolys; i++) {
          num_pts = ipolys[j++];
	  num_polys++;
        
          for (k = 0; k < num_pts; k++) {
            v = ipolys[j++];
            polyi[k] = iverts_num[v];
            }

          polys.push_back(polyi[0]);
          polys.push_back(polyi[2]);
          polys.push_back(polyi[1]);
	  }

        if (pindex > (conn_size - 20)) {
          conn_size += conn_inc;
          num_pmalloc += 1;
	  }

        /*  update adjacent cells.  */

        for (i = 1; i <= num_edges; i++) {
          edge_num = edges[i];
          v1 = edge_list[edge_num][0];
          v2 = edge_list[edge_num][1];

          for (j = 0; j < 3; j++) {
            new_xi = xi + cell_table[v1][v2].list[j].i;
            new_yi = yi + cell_table[v1][v2].list[j].j;
            new_zi = zi + cell_table[v1][v2].list[j].k;

            if ((new_xi < i_end) && (new_xi >= i_begin) &&
                (new_yi < j_end) && (new_yi >= j_begin) &&
                (new_zi < k_end) && (new_zi >= k_begin)) {
              new_cell = new_xi + xdim*(new_yi + ydim*new_zi);
              found = true;

              if (found) {
                /*
                fprintf (stderr, " >>>>>> add cell [%d] \n", new_cell);
                */

                if (!grid_SubsetBitIsSet(num_blocks, blocks, new_cell)) {
                  new_i_list[new_list_size] = new_xi;
                  new_j_list[new_list_size] = new_yi;
                  new_k_list[new_list_size] = new_zi;
                  new_list_size++;
                  grid_SubsetBitSet (num_blocks, blocks, new_cell);
                  }
                }
              }
            }
          }

        }  /*  if ((loc > 0) && (loc < 255))  */
      }  /*  for l < list_size  */
    }  /*  while new_list_size  */


  //delete edge_pt_table;


  /*  free block ptrs.  */

  while (block_list != NULL) {
    //dm_MemFree (0, block_list->ptr_block);
    block_list = block_list->next_block;
    }

#ifdef debug_grid_IsoSeededExtract 
  fprintf (stderr,"\n  num_polys  [%d]  \n", num_polys); 
  fprintf (stderr,"  num_verts  [%d]  \n", num_verts); 
  fprintf (stderr,"  pindex     [%d]  \n", pindex); 
  fprintf (stderr,"  num malloc [%d]  \n", num_malloc);
  fprintf (stderr,"  num poly malloc [%d]  \n", num_pmalloc);
  fprintf (stderr,"  num vert malloc [%d]  \n", num_vmalloc);
#endif

  }


//*============================================================*
//*==========          grid_IsoInit                  ==========*
//*============================================================*
//  initialiaze the tables used to extract an isosurface from  *
//  a cell.                                                    *

static void
grid_IsoInit()
  {

  int i, j, k;

  if (iso_init) {
    return;
    }

  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      for (k = 0; k < 3; k++) {
        cell_table[i][j].list[k].i = 0; 
        cell_table[i][j].list[k].j = 0; 
        cell_table[i][j].list[k].k = 0; 
	}
      }
    }


  /*  (0, 1)  */

  i = 0;
  j = 1;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 1; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = -1; 
  cell_table[i][j].list[1].k = 1; 

  cell_table[i][j].list[2].i = 0; 
  cell_table[i][j].list[2].j = -1; 
  cell_table[i][j].list[2].k = 0; 


  /*  (0, 3)  */

  i = 0;
  j = 3;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 1; 

  cell_table[i][j].list[1].i = -1; 
  cell_table[i][j].list[1].j = 0; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = -1; 
  cell_table[i][j].list[2].j = 0; 
  cell_table[i][j].list[2].k = 1; 


  /*  (0, 4)  */

  i = 0;
  j = 4;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = -1; 
  cell_table[i][j].list[0].k = 0; 

  cell_table[i][j].list[1].i = -1; 
  cell_table[i][j].list[1].j = -1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = -1; 
  cell_table[i][j].list[2].j = 0; 
  cell_table[i][j].list[2].k = 0; 


  /*  (1, 2)  */

  i = 1;
  j = 2;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 1; 

  cell_table[i][j].list[1].i = 1; 
  cell_table[i][j].list[1].j = 0; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 1; 
  cell_table[i][j].list[2].j = 0; 
  cell_table[i][j].list[2].k = 1; 


  /*  (1, 5)  */

  i = 1;
  j = 5;
  cell_table[i][j].list[0].i = 1; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 0; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = -1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 1; 
  cell_table[i][j].list[2].j = -1; 
  cell_table[i][j].list[2].k = 0; 


  /*  (2, 3)  */

  i = 2;
  j = 3;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 1; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = 1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 0; 
  cell_table[i][j].list[2].j = 1; 
  cell_table[i][j].list[2].k = 1; 


  /*  (2, 6)  */

  i = 2;
  j = 6;
  cell_table[i][j].list[0].i = 1; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 0; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = 1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 1; 
  cell_table[i][j].list[2].j = 1; 
  cell_table[i][j].list[2].k = 0; 


  /*  (3, 7)  */

  i = 3;
  j = 7;
  cell_table[i][j].list[0].i = -1; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = 0; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = 1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = -1; 
  cell_table[i][j].list[2].j = 1; 
  cell_table[i][j].list[2].k = 0; 


  /*  (4, 5)  */

  i = 4;
  j = 5;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = -1; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = -1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 0; 
  cell_table[i][j].list[2].j = -1; 
  cell_table[i][j].list[2].k = -1; 


  /*  (4, 7)  */

  i = 4;
  j = 7;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = -1; 

  cell_table[i][j].list[1].i = -1; 
  cell_table[i][j].list[1].j = 0; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = -1; 
  cell_table[i][j].list[2].j = 0; 
  cell_table[i][j].list[2].k = -1; 


  /*  (5, 6)  */

  i = 5;
  j = 6;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = -1; 

  cell_table[i][j].list[1].i = 1; 
  cell_table[i][j].list[1].j = 0; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 1; 
  cell_table[i][j].list[2].j = 0; 
  cell_table[i][j].list[2].k = -1; 


  /*  (6, 7)  */

  i = 6;
  j = 7;
  cell_table[i][j].list[0].i = 0; 
  cell_table[i][j].list[0].j = 0; 
  cell_table[i][j].list[0].k = -1; 

  cell_table[i][j].list[1].i = 0; 
  cell_table[i][j].list[1].j = 1; 
  cell_table[i][j].list[1].k = 0; 

  cell_table[i][j].list[2].i = 0; 
  cell_table[i][j].list[2].j = 1; 
  cell_table[i][j].list[2].k = -1; 


  /*  define isosurface table.  */

  iso_edge_table[0] = case_0_edges; 
  iso_poly_table[0] = case_0_polys; 

  iso_edge_table[1] = case_1_edges; 
  iso_poly_table[1] = case_1_polys; 

  iso_edge_table[2] = case_2_edges; 
  iso_poly_table[2] = case_2_polys; 

  iso_edge_table[3] = case_3_edges; 
  iso_poly_table[3] = case_3_polys; 

  iso_edge_table[4] = case_4_edges; 
  iso_poly_table[4] = case_4_polys; 

  iso_edge_table[5] = case_5_edges; 
  iso_poly_table[5] = case_5_polys; 

  iso_edge_table[6] = case_6_edges; 
  iso_poly_table[6] = case_6_polys; 

  iso_edge_table[7] = case_7_edges; 
  iso_poly_table[7] = case_7_polys; 

  iso_edge_table[8] = case_8_edges; 
  iso_poly_table[8] = case_8_polys; 

  iso_edge_table[9] = case_9_edges; 
  iso_poly_table[9] = case_9_polys; 

  iso_edge_table[10] = case_10_edges; 
  iso_poly_table[10] = case_10_polys; 

  iso_edge_table[11] = case_11_edges; 
  iso_poly_table[11] = case_11_polys; 

  iso_edge_table[12] = case_12_edges; 
  iso_poly_table[12] = case_12_polys; 

  iso_edge_table[13] = case_13_edges; 
  iso_poly_table[13] = case_13_polys; 

  iso_edge_table[14] = case_14_edges; 
  iso_poly_table[14] = case_14_polys; 

  iso_edge_table[15] = case_15_edges; 
  iso_poly_table[15] = case_15_polys; 

  iso_edge_table[16] = case_16_edges; 
  iso_poly_table[16] = case_16_polys; 

  iso_edge_table[17] = case_17_edges; 
  iso_poly_table[17] = case_17_polys; 

  iso_edge_table[18] = case_18_edges; 
  iso_poly_table[18] = case_18_polys; 

  iso_edge_table[19] = case_19_edges; 
  iso_poly_table[19] = case_19_polys; 

  iso_edge_table[20] = case_20_edges; 
  iso_poly_table[20] = case_20_polys; 

  iso_edge_table[21] = case_21_edges; 
  iso_poly_table[21] = case_21_polys; 

  iso_edge_table[22] = case_22_edges; 
  iso_poly_table[22] = case_22_polys; 

  iso_edge_table[23] = case_23_edges; 
  iso_poly_table[23] = case_23_polys; 

  iso_edge_table[24] = case_24_edges; 
  iso_poly_table[24] = case_24_polys; 

  iso_edge_table[25] = case_25_edges; 
  iso_poly_table[25] = case_25_polys; 

  iso_edge_table[26] = case_26_edges; 
  iso_poly_table[26] = case_26_polys; 

  iso_edge_table[27] = case_27_edges; 
  iso_poly_table[27] = case_27_polys; 

  iso_edge_table[28] = case_28_edges; 
  iso_poly_table[28] = case_28_polys; 

  iso_edge_table[29] = case_29_edges; 
  iso_poly_table[29] = case_29_polys; 

  iso_edge_table[30] = case_30_edges; 
  iso_poly_table[30] = case_30_polys; 

  iso_edge_table[31] = case_31_edges; 
  iso_poly_table[31] = case_31_polys; 

  iso_edge_table[32] = case_32_edges; 
  iso_poly_table[32] = case_32_polys; 

  iso_edge_table[33] = case_33_edges; 
  iso_poly_table[33] = case_33_polys; 

  iso_edge_table[34] = case_34_edges; 
  iso_poly_table[34] = case_34_polys; 

  iso_edge_table[35] = case_35_edges; 
  iso_poly_table[35] = case_35_polys; 

  iso_edge_table[36] = case_36_edges; 
  iso_poly_table[36] = case_36_polys; 

  iso_edge_table[37] = case_37_edges; 
  iso_poly_table[37] = case_37_polys; 

  iso_edge_table[38] = case_38_edges; 
  iso_poly_table[38] = case_38_polys; 

  iso_edge_table[39] = case_39_edges; 
  iso_poly_table[39] = case_39_polys; 

  iso_edge_table[40] = case_40_edges; 
  iso_poly_table[40] = case_40_polys; 

  iso_edge_table[41] = case_41_edges; 
  iso_poly_table[41] = case_41_polys; 

  iso_edge_table[42] = case_42_edges; 
  iso_poly_table[42] = case_42_polys; 

  iso_edge_table[43] = case_43_edges; 
  iso_poly_table[43] = case_43_polys; 

  iso_edge_table[44] = case_44_edges; 
  iso_poly_table[44] = case_44_polys; 

  iso_edge_table[45] = case_45_edges; 
  iso_poly_table[45] = case_45_polys; 

  iso_edge_table[46] = case_46_edges; 
  iso_poly_table[46] = case_46_polys; 

  iso_edge_table[47] = case_47_edges; 
  iso_poly_table[47] = case_47_polys; 

  iso_edge_table[48] = case_48_edges; 
  iso_poly_table[48] = case_48_polys; 

  iso_edge_table[49] = case_49_edges; 
  iso_poly_table[49] = case_49_polys; 

  iso_edge_table[50] = case_50_edges; 
  iso_poly_table[50] = case_50_polys; 

  iso_edge_table[51] = case_51_edges; 
  iso_poly_table[51] = case_51_polys; 

  iso_edge_table[52] = case_52_edges; 
  iso_poly_table[52] = case_52_polys; 

  iso_edge_table[53] = case_53_edges; 
  iso_poly_table[53] = case_53_polys; 

  iso_edge_table[54] = case_54_edges; 
  iso_poly_table[54] = case_54_polys; 

  iso_edge_table[55] = case_55_edges; 
  iso_poly_table[55] = case_55_polys; 

  iso_edge_table[56] = case_56_edges; 
  iso_poly_table[56] = case_56_polys; 

  iso_edge_table[57] = case_57_edges; 
  iso_poly_table[57] = case_57_polys; 

  iso_edge_table[58] = case_58_edges; 
  iso_poly_table[58] = case_58_polys; 

  iso_edge_table[59] = case_59_edges; 
  iso_poly_table[59] = case_59_polys; 

  iso_edge_table[60] = case_60_edges; 
  iso_poly_table[60] = case_60_polys; 

  iso_edge_table[61] = case_61_edges; 
  iso_poly_table[61] = case_61_polys; 

  iso_edge_table[62] = case_62_edges; 
  iso_poly_table[62] = case_62_polys; 

  iso_edge_table[63] = case_63_edges; 
  iso_poly_table[63] = case_63_polys; 

  iso_edge_table[64] = case_64_edges; 
  iso_poly_table[64] = case_64_polys; 

  iso_edge_table[65] = case_65_edges; 
  iso_poly_table[65] = case_65_polys; 

  iso_edge_table[66] = case_66_edges; 
  iso_poly_table[66] = case_66_polys; 

  iso_edge_table[67] = case_67_edges; 
  iso_poly_table[67] = case_67_polys; 

  iso_edge_table[68] = case_68_edges; 
  iso_poly_table[68] = case_68_polys; 

  iso_edge_table[69] = case_69_edges; 
  iso_poly_table[69] = case_69_polys; 

  iso_edge_table[70] = case_70_edges; 
  iso_poly_table[70] = case_70_polys; 

  iso_edge_table[71] = case_71_edges; 
  iso_poly_table[71] = case_71_polys; 

  iso_edge_table[72] = case_72_edges; 
  iso_poly_table[72] = case_72_polys; 

  iso_edge_table[73] = case_73_edges; 
  iso_poly_table[73] = case_73_polys; 

  iso_edge_table[74] = case_74_edges; 
  iso_poly_table[74] = case_74_polys; 

  iso_edge_table[75] = case_75_edges; 
  iso_poly_table[75] = case_75_polys; 

  iso_edge_table[76] = case_76_edges; 
  iso_poly_table[76] = case_76_polys; 

  iso_edge_table[77] = case_77_edges; 
  iso_poly_table[77] = case_77_polys; 

  iso_edge_table[78] = case_78_edges; 
  iso_poly_table[78] = case_78_polys; 

  iso_edge_table[79] = case_79_edges; 
  iso_poly_table[79] = case_79_polys; 

  iso_edge_table[80] = case_80_edges; 
  iso_poly_table[80] = case_80_polys; 

  iso_edge_table[81] = case_81_edges; 
  iso_poly_table[81] = case_81_polys; 

  iso_edge_table[82] = case_82_edges; 
  iso_poly_table[82] = case_82_polys; 

  iso_edge_table[83] = case_83_edges; 
  iso_poly_table[83] = case_83_polys; 

  iso_edge_table[84] = case_84_edges; 
  iso_poly_table[84] = case_84_polys; 

  iso_edge_table[85] = case_85_edges; 
  iso_poly_table[85] = case_85_polys; 

  iso_edge_table[86] = case_86_edges; 
  iso_poly_table[86] = case_86_polys; 

  iso_edge_table[87] = case_87_edges; 
  iso_poly_table[87] = case_87_polys; 

  iso_edge_table[88] = case_88_edges; 
  iso_poly_table[88] = case_88_polys; 

  iso_edge_table[89] = case_89_edges; 
  iso_poly_table[89] = case_89_polys; 

  iso_edge_table[90] = case_90_edges; 
  iso_poly_table[90] = case_90_polys; 

  iso_edge_table[91] = case_91_edges; 
  iso_poly_table[91] = case_91_polys; 

  iso_edge_table[92] = case_92_edges; 
  iso_poly_table[92] = case_92_polys; 

  iso_edge_table[93] = case_93_edges; 
  iso_poly_table[93] = case_93_polys; 

  iso_edge_table[94] = case_94_edges; 
  iso_poly_table[94] = case_94_polys; 

  iso_edge_table[95] = case_95_edges; 
  iso_poly_table[95] = case_95_polys; 

  iso_edge_table[96] = case_96_edges; 
  iso_poly_table[96] = case_96_polys; 

  iso_edge_table[97] = case_97_edges; 
  iso_poly_table[97] = case_97_polys; 

  iso_edge_table[98] = case_98_edges; 
  iso_poly_table[98] = case_98_polys; 

  iso_edge_table[99] = case_99_edges; 
  iso_poly_table[99] = case_99_polys; 

  iso_edge_table[100] = case_100_edges; 
  iso_poly_table[100] = case_100_polys; 

  iso_edge_table[101] = case_101_edges; 
  iso_poly_table[101] = case_101_polys; 

  iso_edge_table[102] = case_102_edges; 
  iso_poly_table[102] = case_102_polys; 

  iso_edge_table[103] = case_103_edges; 
  iso_poly_table[103] = case_103_polys; 

  iso_edge_table[104] = case_104_edges; 
  iso_poly_table[104] = case_104_polys; 

  iso_edge_table[105] = case_105_edges; 
  iso_poly_table[105] = case_105_polys; 

  iso_edge_table[106] = case_106_edges; 
  iso_poly_table[106] = case_106_polys; 

  iso_edge_table[107] = case_107_edges; 
  iso_poly_table[107] = case_107_polys; 

  iso_edge_table[108] = case_108_edges; 
  iso_poly_table[108] = case_108_polys; 

  iso_edge_table[109] = case_109_edges; 
  iso_poly_table[109] = case_109_polys; 

  iso_edge_table[110] = case_110_edges; 
  iso_poly_table[110] = case_110_polys; 

  iso_edge_table[111] = case_111_edges; 
  iso_poly_table[111] = case_111_polys; 

  iso_edge_table[112] = case_112_edges; 
  iso_poly_table[112] = case_112_polys; 

  iso_edge_table[113] = case_113_edges; 
  iso_poly_table[113] = case_113_polys; 

  iso_edge_table[114] = case_114_edges; 
  iso_poly_table[114] = case_114_polys; 

  iso_edge_table[115] = case_115_edges; 
  iso_poly_table[115] = case_115_polys; 

  iso_edge_table[116] = case_116_edges; 
  iso_poly_table[116] = case_116_polys; 

  iso_edge_table[117] = case_117_edges; 
  iso_poly_table[117] = case_117_polys; 

  iso_edge_table[118] = case_118_edges; 
  iso_poly_table[118] = case_118_polys; 

  iso_edge_table[119] = case_119_edges; 
  iso_poly_table[119] = case_119_polys; 

  iso_edge_table[120] = case_120_edges; 
  iso_poly_table[120] = case_120_polys; 

  iso_edge_table[121] = case_121_edges; 
  iso_poly_table[121] = case_121_polys; 

  iso_edge_table[122] = case_122_edges; 
  iso_poly_table[122] = case_122_polys; 

  iso_edge_table[123] = case_123_edges; 
  iso_poly_table[123] = case_123_polys; 

  iso_edge_table[124] = case_124_edges; 
  iso_poly_table[124] = case_124_polys; 

  iso_edge_table[125] = case_125_edges; 
  iso_poly_table[125] = case_125_polys; 

  iso_edge_table[126] = case_126_edges; 
  iso_poly_table[126] = case_126_polys; 

  iso_edge_table[127] = case_127_edges; 
  iso_poly_table[127] = case_127_polys; 

  iso_edge_table[128] = case_128_edges; 
  iso_poly_table[128] = case_128_polys; 

  iso_edge_table[129] = case_129_edges; 
  iso_poly_table[129] = case_129_polys; 

  iso_edge_table[130] = case_130_edges; 
  iso_poly_table[130] = case_130_polys; 

  iso_edge_table[131] = case_131_edges; 
  iso_poly_table[131] = case_131_polys; 

  iso_edge_table[132] = case_132_edges; 
  iso_poly_table[132] = case_132_polys; 

  iso_edge_table[133] = case_133_edges; 
  iso_poly_table[133] = case_133_polys; 

  iso_edge_table[134] = case_134_edges; 
  iso_poly_table[134] = case_134_polys; 

  iso_edge_table[135] = case_135_edges; 
  iso_poly_table[135] = case_135_polys; 

  iso_edge_table[136] = case_136_edges; 
  iso_poly_table[136] = case_136_polys; 

  iso_edge_table[137] = case_137_edges; 
  iso_poly_table[137] = case_137_polys; 

  iso_edge_table[138] = case_138_edges; 
  iso_poly_table[138] = case_138_polys; 

  iso_edge_table[139] = case_139_edges; 
  iso_poly_table[139] = case_139_polys; 

  iso_edge_table[140] = case_140_edges; 
  iso_poly_table[140] = case_140_polys; 

  iso_edge_table[141] = case_141_edges; 
  iso_poly_table[141] = case_141_polys; 

  iso_edge_table[142] = case_142_edges; 
  iso_poly_table[142] = case_142_polys; 

  iso_edge_table[143] = case_143_edges; 
  iso_poly_table[143] = case_143_polys; 

  iso_edge_table[144] = case_144_edges; 
  iso_poly_table[144] = case_144_polys; 

  iso_edge_table[145] = case_145_edges; 
  iso_poly_table[145] = case_145_polys; 

  iso_edge_table[146] = case_146_edges; 
  iso_poly_table[146] = case_146_polys; 

  iso_edge_table[147] = case_147_edges; 
  iso_poly_table[147] = case_147_polys; 

  iso_edge_table[148] = case_148_edges; 
  iso_poly_table[148] = case_148_polys; 

  iso_edge_table[149] = case_149_edges; 
  iso_poly_table[149] = case_149_polys; 

  iso_edge_table[150] = case_150_edges; 
  iso_poly_table[150] = case_150_polys; 

  iso_edge_table[151] = case_151_edges; 
  iso_poly_table[151] = case_151_polys; 

  iso_edge_table[152] = case_152_edges; 
  iso_poly_table[152] = case_152_polys; 

  iso_edge_table[153] = case_153_edges; 
  iso_poly_table[153] = case_153_polys; 

  iso_edge_table[154] = case_154_edges; 
  iso_poly_table[154] = case_154_polys; 

  iso_edge_table[155] = case_155_edges; 
  iso_poly_table[155] = case_155_polys; 

  iso_edge_table[156] = case_156_edges; 
  iso_poly_table[156] = case_156_polys; 

  iso_edge_table[157] = case_157_edges; 
  iso_poly_table[157] = case_157_polys; 

  iso_edge_table[158] = case_158_edges; 
  iso_poly_table[158] = case_158_polys; 

  iso_edge_table[159] = case_159_edges; 
  iso_poly_table[159] = case_159_polys; 

  iso_edge_table[160] = case_160_edges; 
  iso_poly_table[160] = case_160_polys; 

  iso_edge_table[161] = case_161_edges; 
  iso_poly_table[161] = case_161_polys; 

  iso_edge_table[162] = case_162_edges; 
  iso_poly_table[162] = case_162_polys; 

  iso_edge_table[163] = case_163_edges; 
  iso_poly_table[163] = case_163_polys; 

  iso_edge_table[164] = case_164_edges; 
  iso_poly_table[164] = case_164_polys; 

  iso_edge_table[165] = case_165_edges; 
  iso_poly_table[165] = case_165_polys; 

  iso_edge_table[166] = case_166_edges; 
  iso_poly_table[166] = case_166_polys; 

  iso_edge_table[167] = case_167_edges; 
  iso_poly_table[167] = case_167_polys; 

  iso_edge_table[168] = case_168_edges; 
  iso_poly_table[168] = case_168_polys; 

  iso_edge_table[169] = case_169_edges; 
  iso_poly_table[169] = case_169_polys; 

  iso_edge_table[170] = case_170_edges; 
  iso_poly_table[170] = case_170_polys; 

  iso_edge_table[171] = case_171_edges; 
  iso_poly_table[171] = case_171_polys; 

  iso_edge_table[172] = case_172_edges; 
  iso_poly_table[172] = case_172_polys; 

  iso_edge_table[173] = case_173_edges; 
  iso_poly_table[173] = case_173_polys; 

  iso_edge_table[174] = case_174_edges; 
  iso_poly_table[174] = case_174_polys; 

  iso_edge_table[175] = case_175_edges; 
  iso_poly_table[175] = case_175_polys; 

  iso_edge_table[176] = case_176_edges; 
  iso_poly_table[176] = case_176_polys; 

  iso_edge_table[177] = case_177_edges; 
  iso_poly_table[177] = case_177_polys; 

  iso_edge_table[178] = case_178_edges; 
  iso_poly_table[178] = case_178_polys; 

  iso_edge_table[179] = case_179_edges; 
  iso_poly_table[179] = case_179_polys; 

  iso_edge_table[180] = case_180_edges; 
  iso_poly_table[180] = case_180_polys; 

  iso_edge_table[181] = case_181_edges; 
  iso_poly_table[181] = case_181_polys; 

  iso_edge_table[182] = case_182_edges; 
  iso_poly_table[182] = case_182_polys; 

  iso_edge_table[183] = case_183_edges; 
  iso_poly_table[183] = case_183_polys; 

  iso_edge_table[184] = case_184_edges; 
  iso_poly_table[184] = case_184_polys; 

  iso_edge_table[185] = case_185_edges; 
  iso_poly_table[185] = case_185_polys; 

  iso_edge_table[186] = case_186_edges; 
  iso_poly_table[186] = case_186_polys; 

  iso_edge_table[187] = case_187_edges; 
  iso_poly_table[187] = case_187_polys; 

  iso_edge_table[188] = case_188_edges; 
  iso_poly_table[188] = case_188_polys; 

  iso_edge_table[189] = case_189_edges; 
  iso_poly_table[189] = case_189_polys; 

  iso_edge_table[190] = case_190_edges; 
  iso_poly_table[190] = case_190_polys; 

  iso_edge_table[191] = case_191_edges; 
  iso_poly_table[191] = case_191_polys; 

  iso_edge_table[192] = case_192_edges; 
  iso_poly_table[192] = case_192_polys; 

  iso_edge_table[193] = case_193_edges; 
  iso_poly_table[193] = case_193_polys; 

  iso_edge_table[194] = case_194_edges; 
  iso_poly_table[194] = case_194_polys; 

  iso_edge_table[195] = case_195_edges; 
  iso_poly_table[195] = case_195_polys; 

  iso_edge_table[196] = case_196_edges; 
  iso_poly_table[196] = case_196_polys; 

  iso_edge_table[197] = case_197_edges; 
  iso_poly_table[197] = case_197_polys; 

  iso_edge_table[198] = case_198_edges; 
  iso_poly_table[198] = case_198_polys; 

  iso_edge_table[199] = case_199_edges; 
  iso_poly_table[199] = case_199_polys; 

  iso_edge_table[200] = case_200_edges; 
  iso_poly_table[200] = case_200_polys; 

  iso_edge_table[201] = case_201_edges; 
  iso_poly_table[201] = case_201_polys; 

  iso_edge_table[202] = case_202_edges; 
  iso_poly_table[202] = case_202_polys; 

  iso_edge_table[203] = case_203_edges; 
  iso_poly_table[203] = case_203_polys; 

  iso_edge_table[204] = case_204_edges; 
  iso_poly_table[204] = case_204_polys; 

  iso_edge_table[205] = case_205_edges; 
  iso_poly_table[205] = case_205_polys; 

  iso_edge_table[206] = case_206_edges; 
  iso_poly_table[206] = case_206_polys; 

  iso_edge_table[207] = case_207_edges; 
  iso_poly_table[207] = case_207_polys; 

  iso_edge_table[208] = case_208_edges; 
  iso_poly_table[208] = case_208_polys; 

  iso_edge_table[209] = case_209_edges; 
  iso_poly_table[209] = case_209_polys; 

  iso_edge_table[210] = case_210_edges; 
  iso_poly_table[210] = case_210_polys; 

  iso_edge_table[211] = case_211_edges; 
  iso_poly_table[211] = case_211_polys; 

  iso_edge_table[212] = case_212_edges; 
  iso_poly_table[212] = case_212_polys; 

  iso_edge_table[213] = case_213_edges; 
  iso_poly_table[213] = case_213_polys; 

  iso_edge_table[214] = case_214_edges; 
  iso_poly_table[214] = case_214_polys; 

  iso_edge_table[215] = case_215_edges; 
  iso_poly_table[215] = case_215_polys; 

  iso_edge_table[216] = case_216_edges; 
  iso_poly_table[216] = case_216_polys; 

  iso_edge_table[217] = case_217_edges; 
  iso_poly_table[217] = case_217_polys; 

  iso_edge_table[218] = case_218_edges; 
  iso_poly_table[218] = case_218_polys; 

  iso_edge_table[219] = case_219_edges; 
  iso_poly_table[219] = case_219_polys; 

  iso_edge_table[220] = case_220_edges; 
  iso_poly_table[220] = case_220_polys; 

  iso_edge_table[221] = case_221_edges; 
  iso_poly_table[221] = case_221_polys; 

  iso_edge_table[222] = case_222_edges; 
  iso_poly_table[222] = case_222_polys; 

  iso_edge_table[223] = case_223_edges; 
  iso_poly_table[223] = case_223_polys; 

  iso_edge_table[224] = case_224_edges; 
  iso_poly_table[224] = case_224_polys; 

  iso_edge_table[225] = case_225_edges; 
  iso_poly_table[225] = case_225_polys; 

  iso_edge_table[226] = case_226_edges; 
  iso_poly_table[226] = case_226_polys; 

  iso_edge_table[227] = case_227_edges; 
  iso_poly_table[227] = case_227_polys; 

  iso_edge_table[228] = case_228_edges; 
  iso_poly_table[228] = case_228_polys; 

  iso_edge_table[229] = case_229_edges; 
  iso_poly_table[229] = case_229_polys; 

  iso_edge_table[230] = case_230_edges; 
  iso_poly_table[230] = case_230_polys; 

  iso_edge_table[231] = case_231_edges; 
  iso_poly_table[231] = case_231_polys; 

  iso_edge_table[232] = case_232_edges; 
  iso_poly_table[232] = case_232_polys; 

  iso_edge_table[233] = case_233_edges; 
  iso_poly_table[233] = case_233_polys; 

  iso_edge_table[234] = case_234_edges; 
  iso_poly_table[234] = case_234_polys; 

  iso_edge_table[235] = case_235_edges; 
  iso_poly_table[235] = case_235_polys; 

  iso_edge_table[236] = case_236_edges; 
  iso_poly_table[236] = case_236_polys; 

  iso_edge_table[237] = case_237_edges; 
  iso_poly_table[237] = case_237_polys; 

  iso_edge_table[238] = case_238_edges; 
  iso_poly_table[238] = case_238_polys; 

  iso_edge_table[239] = case_239_edges; 
  iso_poly_table[239] = case_239_polys; 

  iso_edge_table[240] = case_240_edges; 
  iso_poly_table[240] = case_240_polys; 

  iso_edge_table[241] = case_241_edges; 
  iso_poly_table[241] = case_241_polys; 

  iso_edge_table[242] = case_242_edges; 
  iso_poly_table[242] = case_242_polys; 

  iso_edge_table[243] = case_243_edges; 
  iso_poly_table[243] = case_243_polys; 

  iso_edge_table[244] = case_244_edges; 
  iso_poly_table[244] = case_244_polys; 

  iso_edge_table[245] = case_245_edges; 
  iso_poly_table[245] = case_245_polys; 

  iso_edge_table[246] = case_246_edges; 
  iso_poly_table[246] = case_246_polys; 

  iso_edge_table[247] = case_247_edges; 
  iso_poly_table[247] = case_247_polys; 

  iso_edge_table[248] = case_248_edges; 
  iso_poly_table[248] = case_248_polys; 

  iso_edge_table[249] = case_249_edges; 
  iso_poly_table[249] = case_249_polys; 

  iso_edge_table[250] = case_250_edges; 
  iso_poly_table[250] = case_250_polys; 

  iso_edge_table[251] = case_251_edges; 
  iso_poly_table[251] = case_251_polys; 

  iso_edge_table[252] = case_252_edges; 
  iso_poly_table[252] = case_252_polys; 

  iso_edge_table[253] = case_253_edges; 
  iso_poly_table[253] = case_253_polys; 

  iso_edge_table[254] = case_254_edges; 
  iso_poly_table[254] = case_254_polys; 

  iso_init = true;
  }

//*============================================================*
//*==========          SubsetCreate                  ==========*
//*============================================================*
// create the set of blocks for bit masking.                         

static void
grid_SubsetCreate (int size, int *p_num_blocks, Subset **p_blocks)
  {

  int num_blocks;

  Subset *blocks;

  num_blocks = size / SubsetSize + 1;
  blocks = new Subset[num_blocks];

  for (int i = 0; i < num_blocks; i++) {
    blocks[i] = 0;
    }

  *p_num_blocks = num_blocks;
  *p_blocks = blocks;
  }

//*============================================================*
//*==========          grid_SubsetBitSet             ==========*
//*============================================================*
// set bits.                                     

static void
grid_SubsetBitSet (int num_blocks, Subset *blocks, int n)
  {

  Subset block, mask;

  int bit_num;

  int block_num;

  block_num = (int)(n / SubsetSize + 0.5);
  block = blocks[block_num];

  bit_num = n % SubsetSize;
  mask = 1;
  mask <<= (bit_num);
  block = block | mask;
  blocks[block_num] = block;
  }

//*============================================================*
//*==========          grid_SubsetBitIsSet           ==========*
//*============================================================*
// check is a bit is set.

static bool 
grid_SubsetBitIsSet (int num_blocks, Subset *blocks, int n)
  {

  Subset block, mask;

  int bit_num;

  int block_num;

  block_num = (int)(n / SubsetSize + 0.5);
  block = blocks[block_num];

  bit_num = n % SubsetSize;
  mask = 1;
  mask <<= (bit_num);

  return (block & mask);
  }

//*============================================================*
//*==========             displaySlice               ==========*
//*============================================================*

void
PmGrid::displaySlice(const int id, const int dim, const int index, const float vmin,
                     const float vmax)
  {
  int dims[3];
  float cell_sizes[3];
  vector<PmVector3> verts;
  vector<int> polys;
  vector<float> vals;
  PmGraphicsAttributes atts;
  PmGraphicsPolygon *poly;
  int num_tri, num_verts;
  PmVector3 *gverts;
  PmVector3 *gcolors;
  PmConn *conn;
  vector<PmVector3> colors;

  //fprintf (stderr, "\n>>>>>> PmGrid::displaySlice \n");

  dims[0] = idim; dims[1] = jdim; dims[2] = kdim;
  cell_sizes[0] = cell_dx; cell_sizes[1] = cell_dy; cell_sizes[2] = cell_dz;


  // extract slice geometry and data values //

  grid_SliceExtract (dims, cell_sizes, grid_data, dim, index, verts, polys, vals);

  num_verts = verts.size();
  gverts = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gverts[i] = verts[i];
    }

  num_tri = polys.size() / 3;
  conn = new PmConn(3*num_tri);
  int *pconn = conn->vals;

  for (int i = 0; i < 3*num_tri; i++) {
    pconn[i] = polys[i];
    }

  // map data values into a color map //

  grid_MapData(vmin, vmax, vals, colors);
  gcolors = new PmVector3[num_verts];

  for (int i = 0; i < num_verts; i++) {
    gcolors[i] = colors[i];
    }

  atts.setShadingType(shading_type);
  atts.setDisplayType(display_type);
  atts.setColor(color);
  poly = new PmGraphicsPolygon(name, num_tri, conn, 3, num_verts, gverts);
  poly->setAttributes(atts);
  poly->setColors(num_verts, gcolors);
  poly->display();

  /*
  color.set(1,0,0);
  atts.setColor(color);
  atts.setMarker(false);
  atts.setScale(2.0);

  PmGraphicsPoint *geom = new PmGraphicsPoint(name, num_verts, gverts);
  geom->setAttributes(atts);
  geom->setColors(gcolors);
  geom->display();
  */
  }

//*============================================================*
//*==========          grid_SliceExtract             ==========*
//*============================================================*

static void
grid_SliceExtract (int dims[3], float cell_sizes[3], float *data, int dim, int index,
                   vector<PmVector3>& verts, vector<int>& polys, vector<float>& vals)
  {
  int xdim, ydim, zdim;
  float dx, dy, dz; 
  int i_begin, i_end, j_begin, j_end, k_begin, k_end;
  int i, j, k;
  float ox, oy, oz, x, y, z;
  int num_rows, num_cols;
  int skip;
  float val;
  float ux, uy, uz, vx, vy, vz, wx, wy, wz;
  int i1, i2;
  PmVector3 vec;

  // set up initial indexes into the grid data //

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];
  dx = cell_sizes[0];
  dy = cell_sizes[1];
  dz = cell_sizes[2];
  ox = oy = oz = 0.0;

  skip = 1;
  i_begin = 0;
  i_end = xdim - 1;
  j_begin = 0;
  j_end = ydim - 1;
  k_begin = 0;
  k_end = zdim - 1;

  ux = dx; 
  uy = 0.0; 
  uz = 0.0; 

  vx = 0.0;
  vy = dy; 
  vz = 0.0;

  wx = 0.0;
  wy = 0.0; 
  wz = dz; 

  if (dim == 1) {
    i = index;
    num_rows = j_end - j_begin + 1;
    num_cols = k_end - k_begin + 1;

    for (j = j_begin; j <= j_end; j++) {
      for (k = k_begin; k <= k_end; k++) {
        val = data[i + xdim*(j + ydim*k)];
        x = ox + i*ux + vx*j + wx*k;
        y = oy + i*uy + vy*j + wy*k;
        z = oz + i*uz + vz*j + wz*k;
        vec.set(x,y,z);
        verts.push_back(vec); 
        vals.push_back(val); 
        }
      }
    }

  else if (dim == 2) {
    j = index;
    num_rows = k_end - k_begin + 1;
    num_cols = i_end - i_begin + 1;

    for (k = k_begin; k <= k_end; k += skip) {
      for (i = i_begin; i <= i_end; i += skip) {
        val = data[i + xdim*(j + ydim*k)];
        x = ox + i*ux + vx*j + wx*k;
        y = oy + i*uy + vy*j + wy*k;
        z = oz + i*uz + vz*j + wz*k;
        vec.set(x,y,z);
        verts.push_back(vec); 
        vals.push_back(val); 
        }
      }
    }

  else if (dim == 3) {
    k = index;
    num_rows = j_end - j_begin + 1;
    num_cols = i_end - i_begin + 1;

    for (j = j_begin; j <= j_end; j += skip) {
      for (i = i_begin; i <= i_end; i += skip) {
        val = data[i + xdim*(j + ydim*k)];
        x = ox + i*ux + vx*j + wx*k;
        y = oy + i*uy + vy*j + wy*k;
        z = oz + i*uz + vz*j + wz*k;
        vec.set(x,y,z);
        verts.push_back(vec);
        vals.push_back(val);
        }
      }
    }

  i1 = 0;
  i2 = num_cols;

  for (i = 0; i < num_cols-1; i++) {
    for (j = 0; j < num_rows-1; j++) {
      polys.push_back(i1);
      polys.push_back(i1+1);
      polys.push_back(i2+1);
    
      polys.push_back(i2);
      polys.push_back(i1);
      polys.push_back(i2+1);

      i1 += 1;
      i2 += 1;
      }

    i1 = (i+1)*num_cols;
    i2 = i1 + num_cols;
    }
  }

//*============================================================*
//*==========            grid_MapData                ==========*
//*============================================================*

static void
grid_MapData(float vmin, float vmax, vector<float>& data, vector<PmVector3>& colors)
  {
  vector<PmVector3> cmap;
  float val, f, dv;
  int ci, num_colors;

  // generate color map //

  grid_ColorMapSpectrumCreate (64, cmap);
  num_colors = cmap.size();
  dv = vmax - vmin;

  if (dv != 0.0) {
    f = (num_colors-1) / dv;
    }
  else{
    f = 0.0; 
    }

  for (unsigned int i = 0; i < data.size(); i++ ) {
    val = data[i];
    ci = (int)(f*(val - vmin));
    if (ci < 0) ci = 0;
    if (ci > num_colors-1) ci = num_colors-1;
    colors.push_back(cmap[ci]);
    }
  }

//*============================================================*
//*==========      grid_ColorMapSpectrumCreate       ==========*
//*============================================================*

static void
grid_ColorMapSpectrumCreate (int num_colors, vector<PmVector3>& colors)
  {

  float inc, r, g, b, hue;

  PmVector3 color;

  int i;

  inc = 240.0 / (float)(num_colors - 1);

  for (i = 0, hue = 240.0; i < num_colors; i++, hue -= inc) {
    grid_ColorMapHsvToRgb (hue, 1.0, 1.0, &r, &g, &b);
    if (r > 1.0) r = 1.0;
    if (r < 0.0) r = 0.0;
    if (g > 1.0) g = 1.0;
    if (g < 0.0) g = 0.0;
    if (b > 1.0) b = 1.0;
    if (b < 0.0) b = 0.0;
    color[0] = r;
    color[1] = g;
    color[2] = b;
    colors.push_back(color);
    }
  }


//*============================================================*
//*==========          grid_ColorMapHsvToRgb         ==========*
//*============================================================*

static void
grid_ColorMapHsvToRgb (float h, float s, float v, float *r, float *g, float *b)
  {

  float f, p, q, t;

  int i;

 /**************
  ***  body  ***
  **************/

  if (s == 0.0) {
    *r = v;
    *g = v;
    *b = v;
    }
  else {
    if (h == 360.0) {
      h = 0.0;
      }

    h = h / 60.0;
    i = (int)h;
    f = h - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s*f));
    t = v * (1.0 - (s * (1.0 - f)));

    switch (i) {
      case 0:
        *r = v, *g = t, *b = p;
      break;

      case 1:
        *r = q, *g = v, *b = p;
      break;

      case 2:
        *r = p, *g = v, *b = t;
      break;

      case 3:
        *r = p, *g = q, *b = v;
      break;

      case 4:
        *r = t, *g = p, *b = v;
      break;

      case 5:
        *r = v, *g = p, *b = q;
      break;
      }
    }
  }

}


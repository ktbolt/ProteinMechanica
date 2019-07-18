
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
//* grid:                    g r i d                           *
//*============================================================*

#ifndef _GRID_PM_H_
#define _GRID_PM_H_

#include "pm/pm.h"
#include "graphics.h"
#include "pm/mth.h"
#include "db.h"

namespace ProteinMechanica {

  class PmDbPmInterface;


// PmGrid     
// ------

 class PM_EXPORT PmGrid {

    public:

       PmGrid(const string name);

       PmGrid(const string name, const int idim, const int jdim, const int kdim,
              const float dx, const float dy, const float dz, float *data);
       ~PmGrid();

       void getDataRange(float& vmin, float& vmax);

       void getExtent (PmExtent& extent);
       void setExtent (const PmExtent& extent);
       void getIndex(PmVector3& pt, int& i, int& j, int& k);
       void interpScalar(PmVector3& pt, int& i, int& j, int& k, float& ival);

       void getName (string& name);

       // graphics //

       void displayIsosurface(const float level, bool seeded, PmVector3 point, 
                              PmDbPmInterface *db);
       void displayPoints(const float vmin, const float vmax, PmDbPmInterface *db);
       void displaySlice(const int id, const int dim, const int index, const float vmin,
                         const float vmax);
       void displayVolume(const float vmin, const float vmax, const int num_samples);
       //void setColor(PmVector3& color) { this->datts.color = color; }
       //void setXform (PmXform& xform);

       void setColor (PmVector3& color) { this->color = color; }
       void setDisplayType (PmGeometryDisplayType type) { this->display_type = type; }
       void setShadingType (PmGeometryShadingType type) { this->shading_type = type; }

    private:
       string name;
       int idim, jdim, kdim;
       float cell_dx, cell_dy, cell_dz;
       float data_vmin, data_vmax, *grid_data;
       PmExtent extent;
       PmXform xform;
       PmVector3 origin;
       PmVector3 basis_u, basis_v, basis_w; 
       void init();
       void getCellData(int i, int j, int k, float cell_data[8]);

       PmVector3 color;
       PmGeometryDisplayType display_type;
       PmGeometryShadingType shading_type;

       void buildGeomName (string& grname);
       void addGeometry (PmGraphicsPoint *geom);
       void getGeometry (string name, PmGraphicsPoint **geom);
       vector<PmGraphicsPoint*>  graphics_point_geometry;
   };

}

#endif

